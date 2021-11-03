/*
 * fwCoreMex.cpp - 
 *
 * The calling syntax is:
 *
 *		C = fwCoreMatlab(X, initC, lambda, iterStart, epsilon, maxIters,...
                innerVerboseInterval, objTol)
 *
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <math.h>       /* exp */
#include <algorithm>    // std::min
#include <thread>
#include <vector>
#include <iostream>

typedef Eigen::Triplet<double> T;
// template<typename T> struct TD;

void multipleMatrixVector(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 
        Eigen::RowMajor>& tXT,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& tCT,
    int start,
    int tIndex,
    std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic, 
        Eigen::RowMajor>>& tG1); 
void computeSumOverColC(Eigen::SparseMatrix<double, Eigen::RowMajor> tC, 
        Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>& tSumOverColC1,
        Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>& tSumOverColC2,
        double epsilon);
double computePhiC(Eigen::SparseMatrix<double, Eigen::RowMajor> tCT, 
        double epsilon);
void fwMain(double *X, double *sr, mwIndex *jcs, mwIndex *irs, mwSize m, 
        mwSize n, size_t bSize, size_t estimatedN, size_t maxIters,
        size_t iterStart, size_t innerVerboseInterval, double lambda, 
        double epsilon, double objTol, bool verbose, mxArray *plhs[]);



/* The gateway function 
function [C] = fwCoreMatlab(X, initC, lambda, iterStart, epsilon, maxIters,...
         innerVerboseInterval, objTol) */
void mexFunction(int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *X;               /* MxN input matrix */
    double *C;               /* MxN input matrix */
    size_t ncols;                   /* size of matrix */
    size_t nrows;                   /* size of matrix */
    double *sr;
    mwIndex *irs,*jcs,j,k;



    /* X  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    X = mxGetDoubles(prhs[0]);
    #else
    X = mxGetPr(prhs[0]);
    #endif
    /* Get the size and pointers to input data */
    nrows =mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);

    C = mxGetPr(prhs[1]);
    /* Check data type dimensions of the matrix input y. */
    if (!(mxIsDouble(prhs[1]))){
        mexErrMsgIdAndTxt( "MATLAB:DKU", "Input argument must be of type double.");
    }
    if (mxGetNumberOfDimensions(prhs[1]) != 2){
        mexErrMsgIdAndTxt( "MATLAB:DKU", "Input argument must be two dimensional\n");
    }
    // Verify the matrix C is infact sparse
    if (!(mxIsSparse(prhs[1]))){
        mexErrMsgIdAndTxt( "MATLAB: DKU", "Input matrix is not sparse\n");
    } else {    
        sr  = mxGetPr(prhs[1]);
        jcs = mxGetJc(prhs[1]);    
        irs = mxGetIr(prhs[1]);
    }

    double lambda = mxGetScalar(prhs[2]);
    size_t iterStart = mxGetScalar(prhs[3]);
    double epsilon = mxGetScalar(prhs[4]);
    size_t maxIters = mxGetScalar(prhs[5]);
    size_t innerVerboseInterval = mxGetScalar(prhs[6]);
    double objTol = mxGetScalar(prhs[7]);
    bool verbose = mxGetScalar(prhs[8]);
    int estimatedN = mxGetScalar(prhs[9]);
    size_t numThreads = mxGetScalar(prhs[10]);
    mexPrintf("Verify your input: \n");
    mexPrintf("lambda = %f \n", lambda);
    mexPrintf("iterStart = %d \n", iterStart);
    mexPrintf("epsilon = %f \n", epsilon);
    mexPrintf("maxIters = %d \n", maxIters);
    mexPrintf("innerVerboseInterval = %d \n", innerVerboseInterval);
    mexPrintf("objTol = %f \n", objTol);
    mexPrintf("estimatedN = %d \n", estimatedN);
    mexPrintf("verbose = %d \n", verbose);
    mexPrintf("Size of X: %d x %d \n", nrows, ncols);
    mexPrintf("numThreads = %d \n", numThreads);


    /* create the output matrix */
    if (estimatedN<0) {
        estimatedN = 100;
        mexPrintf("WARNING: using 100 as estimated N to pre-allocate memory\n");
    }

    mexPrintf("There are %d cpus. Mex will use %d threads \n", Eigen::nbThreads(), numThreads);
    /* call the computational routine */
    fwMain(X,sr,jcs,irs,nrows,ncols,numThreads,estimatedN,maxIters,iterStart,
            innerVerboseInterval,lambda,epsilon,objTol,verbose,plhs);
}


void multipleMatrixVector(
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 
        Eigen::RowMajor>& tXT,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& tCT,
    int start,
    int tIndex,
    std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic, 
        Eigen::RowMajor>>& tG1
    ) {
    tG1[tIndex] = (2*(tCT.row(start+tIndex)*tXT - tXT.row(start+tIndex))) *
        tXT.transpose();
}


void computeSumOverColC(Eigen::SparseMatrix<double, Eigen::RowMajor> tCT, 
        Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>& tSumOverColC1,
        Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>& tSumOverColC2,
        double epsilon) {
    size_t colIndex;
    size_t anchor=0;
    size_t i=0;
    size_t L = tCT.outerSize();
    double tmp;
    Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> eMaxValueInverse(L);

    // tSumOverColC2 = tCT.colwise().maxCoeff();
	for (size_t k=0; k<L; ++k)
    {
        anchor=0;
		for (Eigen::SparseMatrix<double, 
                Eigen::RowMajor>::InnerIterator it(tCT,k); it; ++it) {
            colIndex = it.col();
            tmp = it.value()/epsilon;
            if (tmp > tSumOverColC2(colIndex)) {
                tSumOverColC2(colIndex) = tmp;
            }
		}
    }

    // Compute e^(-maxValue) once for every column, and reuse this value for zero elements
    eMaxValueInverse = (-tSumOverColC2.array()).exp();

	for (size_t k=0; k<L; ++k)
    {
        anchor=0;
		for (Eigen::SparseMatrix<double, 
                Eigen::RowMajor>::InnerIterator it(tCT,k); it; ++it) {
            colIndex = it.col();

            for (i = anchor; i < colIndex; ++i) {
                tSumOverColC1(i) += eMaxValueInverse(i); // e^(0-maxValue)
            }
            // e^(current-value - maxValue)
            tSumOverColC1(colIndex) += exp(it.value()/epsilon - tSumOverColC2(colIndex));

            anchor = colIndex+1;
		}
        // Finish remaining element of current columns
        for (i = anchor; i < L; ++i) {
            tSumOverColC1(i) += eMaxValueInverse(i);
        }
    }
}


/* The computational routine */
void fwMain(double *X, double *sr, mwIndex *jcs, mwIndex *irs, mwSize m, 
        mwSize n, size_t bSize, size_t estimatedN, size_t maxIters,
        size_t iterStart, size_t innerVerboseInterval, double lambda, 
        double epsilon, double objTol, bool verbose, mxArray *plhs[]) {
    mexPrintf("MaxIter from fwMain: %d \n", maxIters);
    mexPrintf("Input size: %d x %d \n", m, n);
    mwSize i,j;
    size_t row, col;
    double val;

    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> tX(X,m,n);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tXT = tX.transpose();

    Eigen::SparseMatrix<double> tZ(m,n);
    Eigen::SparseVector<double,Eigen::RowMajor> tVecTmp(n);

    std::vector<T> tripletList;
    size_t cCapacity = estimatedN*n;
    tripletList.reserve(cCapacity);

    /* Construct Sparse matrix in Eigen */
    for (i = 1; i < n+1; ++i) {
        for (j = (size_t)jcs[i-1]; j < size_t(jcs[i]); ++j) {
            // row = irs[j];
            // col = i-1;
            row = i-1;
            col = irs[j];
            val = sr[j];
            tripletList.push_back(T((size_t)row, (size_t)col, val));
        }
    }
    /* End construction */


    Eigen::SparseMatrix<double, Eigen::RowMajor> tCT(n,n);
    tCT.setFromTriplets(tripletList.begin(), tripletList.end());
    
    Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> tSumOverColC1(n); 
    Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> tSumOverColC2(n); 
    double stepSize = 0;
    std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>> tG1(bSize);
    size_t colIndex, k, iii, anchor, iMin, bIter, cNnz;
    double vMin, currentObj, newObj, relObjChange, fillingRate, fittingError, regValue;
    currentObj = 100000000;
    relObjChange = 100000000;
    double tg2Current = 0.;
    bool foundIt;
    mexPrintf("Start main loop \n");
    for (size_t iter = iterStart; iter <= maxIters; ++iter) {
        stepSize = 2.0/(iter + 1);
        if (verbose) {
            mexPrintf("------------------   Iter: %d   ------------------\n", iter);
            mexPrintf("Stepsize: %f \n", stepSize);
        }

        tSumOverColC1.setZero();
        tSumOverColC2.setZero();
        computeSumOverColC(tCT, tSumOverColC1, tSumOverColC2, epsilon);
        tripletList.clear();
        for (bIter = 0; bIter < n; bIter=bIter+bSize) {
            if (verbose && (bIter % innerVerboseInterval < bSize)) {
                mexPrintf("BatchIter: %d \n", bIter);
            }
            std::thread threads[bSize];
            for (int i = 0; i < std::min(bSize, n-bIter); ++i) {
                threads[i] = std::thread(multipleMatrixVector,
                        std::ref(tXT), 
                        std::ref(tCT),
                        bIter, 
                        i,
                        std::ref(tG1)
                        );
            }

            for (int i = 0; i < std::min(bSize, n-bIter); ++i) {
              threads[i].join();
            }

            for (i = 0; i < std::min(bSize, n-bIter); ++i) {
                k = i+bIter;

                /* BEGIN Computing G2, and also evaluate G at the same time */
                anchor = 0;
                iMin = 0;
                vMin = 100000000.;
                tg2Current = 0;

                for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(tCT,k); it; ++it) {
                    colIndex = it.col();   // col index
                    for (iii = anchor; iii < colIndex; ++iii) {
                        tg2Current = tG1[i](iii) + lambda *exp(-tSumOverColC2(iii)) / tSumOverColC1(iii);
                        if (tg2Current < vMin) {
                            vMin = tg2Current;
                            iMin = iii;
                        }
                    }
                    tg2Current = tG1[i](colIndex) + 
                                 lambda*exp(it.value()/epsilon - tSumOverColC2(iii)) 
                                        / tSumOverColC1(colIndex);

                    if (tg2Current < vMin) {
                        vMin = tg2Current;
                        iMin = colIndex;
                    }
                    anchor = colIndex+1;
                }

                for (iii = anchor; iii < n; ++iii) {
                    tg2Current = tG1[i](iii) +lambda*exp(-tSumOverColC2(iii))/ tSumOverColC1(iii);
                    if (tg2Current < vMin) {
                        vMin = tg2Current;
                        iMin = iii;
                    }
                }
                /* ENDING Computing G2 */

                /* START Update column ck 
                 * This is equivalent to update column k of CNew
                 * */
                tVecTmp = tCT.row(k);
                foundIt = false;
                for (Eigen::SparseVector<double, Eigen::RowMajor>::InnerIterator it(tVecTmp); it; ++it)
                {
                    iii = it.index();
                    if (iii == iMin) {
                        foundIt = true;
                    } else {
                        tripletList.push_back(T(k, (size_t)iii,(double)((1.-stepSize)*it.value())));
                    }
                }
                if (foundIt) {
                    tripletList.push_back(T(k, iMin, (double)((1.-stepSize)*
                                    tVecTmp.coeffRef(iMin) + stepSize)));
                } else {
                    tripletList.push_back(T(k, iMin, (double)stepSize));
                }
                /* END Update column ck */
            }
        }
        tCT.setZero();
        tCT.setFromTriplets(tripletList.begin(), tripletList.end());

        fittingError = (tXT - tCT*tXT).squaredNorm();
        // mexPrintf("DEBUG: sum of tX: %f \n", tX.sum());
        // mexPrintf("DEBUG: sum of tCT: %f \n", tCT.sum());
        regValue = epsilon * ( tSumOverColC1.array().log() + tSumOverColC2.array() - log(n)).sum();
        newObj = fittingError + lambda*regValue;

        cNnz = tCT.nonZeros();
        fillingRate = ((double)cNnz)/cCapacity;

        if (fillingRate > 0.7) {
            cCapacity = (size_t)(1.5*cCapacity);
            mexPrintf("WARNING: Increase capacity of C to: %d \n", cCapacity);
            tripletList.clear();
            tripletList.reserve(cCapacity);
        }

        relObjChange = (currentObj-newObj)/currentObj;
        if (relObjChange<0)
            relObjChange = -relObjChange;
        
        if (verbose) {
            mexPrintf("Current obj: %f \n", currentObj);
            mexPrintf("New obj: %f \n", newObj);
            mexPrintf("Fitting error: %f \n", fittingError);
            mexPrintf("Reg value: %f \n", regValue);
            mexPrintf("Relative change: %f \n", relObjChange);
            mexPrintf("C filling up: %f \t (%d/%d) \n" , 
                        fillingRate, cNnz, cCapacity);
            mexPrintf("\n");
        }
        currentObj = newObj;

        if (relObjChange < objTol) {
            break;
        }
    }

    // It should be tC now, not tCT
    Eigen::SparseMatrix<double, Eigen::ColMajor> tC = tCT.transpose();

    // mexPrintf("Matrix C is: \n");
    // for (i = 0; i < n; ++i) {
    //     for (j = 0; j < n; ++j) {
    //         mexPrintf("%f \t", tC.coeffRef(i, j));
    //     }
    //     mexPrintf("\n");
    // }

    plhs[0] = mxCreateSparse((mwSize)n,(mwSize)n,cNnz,mxREAL);
    double *zSr  = mxGetPr(plhs[0]);
    mwIndex *zIrs = mxGetIr(plhs[0]);
    mwIndex *zJcs = mxGetJc(plhs[0]);
   
    /* Copy tCT to output Z*/
    k = 0;
	for (size_t i=0; i<n; ++i)
    {
        zJcs[i] = k;
		for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it(tC,i); it; ++it) {
            zIrs[k] = it.row();
            zSr[k] = it.value();
            // mexPrintf("C[%d, %d] = %f \n", zIrs[k], zJcs[i], zSr[k]);
            k++;
		}
    }
    zJcs[n] = k;
}


