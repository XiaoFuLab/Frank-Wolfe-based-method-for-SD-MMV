/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2012 The MathWorks, Inc.
 *
 *========================================================*/

#include "mex.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <math.h>       /* exp */
#include <algorithm>    // std::min
#include <thread>
#include <vector>
#include <iostream>

const size_t ESTIMATED_NNZ = 500;
const double EPSILON = 0.01;
typedef Eigen::Triplet<double> T;
// template <typename Derived>;
template<typename T> struct TD;
// void eigenArrayToC(Eigen::MatrixXd *X, double *Z)

void multipleMatrixVector(
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& tXT,
        // Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& tX,
        Eigen::SparseMatrix<double, Eigen::RowMajor>& tCT,
        int start,
        int tIndex,
        std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>& tG1
        ) {

    tG1[tIndex] = (2*(tCT.row(start+tIndex)*tXT - tXT.row(start+tIndex))) * tXT.transpose();
    // tG1[tIndex].setZero(); 
}

/* The computational routine */
void fwMain(double *X, double *sr, mwIndex *jcs, mwIndex *irs, mwSize m, mwSize n, double* Z, size_t bSize)
{
    mwSize i,j;
    size_t row, col;
    double val;

    mexPrintf("m=%d n=%d \n", m, n);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>> tX(X,m,n);
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tXT = tX.transpose();

    Eigen::SparseMatrix<double> tZ(m,n);
    Eigen::SparseVector<double,Eigen::RowMajor> tVecTmp(n);

    std::vector<T> tripletList;
    tripletList.reserve(ESTIMATED_NNZ*n);
    mexPrintf("ESTIMATED_NNZ=%u and n=%u \n", ESTIMATED_NNZ, n);

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
    // mexPrintf("Len : %d \n", tripletList.size());
    tCT.setFromTriplets(tripletList.begin(), tripletList.end());
    
    // mexPrintf("Matrix CT before is: \n");
    // for (i = 0; i < n; ++i) {
    //     for (j = 0; j < n; ++j) {
    //         mexPrintf("%f \t", tCT.coeffRef(i, j));
    //     }
    //     mexPrintf("\n");
    // }

    Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> tSumOverColC(n); 
    double stepSize = 0;

    // Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor> tG1i(n);
    std::vector<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>> tG1(bSize);

    size_t colIndex, k, iii, anchor, iMin, bIter;
    double vMin;
    double tg2Current = 0.;
    bool foundIt;
    mexPrintf("Start main loop \n");
    for (size_t iter = 1; iter <= 1; ++iter) {
        stepSize = 2.0/(iter + 1);

        tSumOverColC.setZero();
        // mexPrintf("sumC = ");
        // for (iii = 0; iii < n; ++iii) {
        //     mexPrintf("%f ", tSumOverColC(iii));
        // }
        // mexPrintf("\n");

        tripletList.clear();
        // Loop row by row of C^T
        for (bIter = 0; bIter < n; bIter=bIter+bSize) {
            mexPrintf("BatchIter: %d \n", bIter);
            // tG1 = (2*(tXT.rows(k)))*tX;
            // tG1 = (2*( tCT.rows(k:bSize) * tXT - tXT.rows(k:bSize) ))*tX;
            // mexPrintf("BatchIter: %d\n", bIter);
            // tG1 = (2*(tCT.block(bIter, 0, std::min(bSize, n-bIter), n) * tXT - 
            //             tXT.block(bIter, 0, std::min(bSize, n-bIter), n)))*tX;
            
            std::thread threads[bSize];
            for (int i = 0; i < std::min(bSize, n-bIter); ++i) {
            //std::cout << "Starting thread " << i << std::endl;
                threads[i] = std::thread(multipleMatrixVector,
                        std::ref(tXT), 
                        std::ref(tCT),
                        bIter, 
                        i,
                        std::ref(tG1)
                        );
            }

            //std::cout << "Calculating...." << std::endl;

            for (int i = 0; i < bSize; ++i) {
              //std::cout << "Joining thread " << i << std::endl;
              threads[i].join();
            }

            for (i = 0; i < bSize; ++i) {
                k = i+bIter;

                // BEGIN Computing G2
                anchor = 0;
                iMin = 0;
                vMin = 100000000.;
                tg2Current = 0;

                for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(tCT,k); it; ++it) {
                    colIndex = it.col();   // col index
                    mexPrintf("CT[%d, %d] = %f \n", it.row(), it.col(), it.value());
                    for (iii = anchor; iii < colIndex; ++iii) {
                        tg2Current = tG1[i](iii) + 1. / tSumOverColC(iii);
                        if (tg2Current < vMin) {
                            vMin = tg2Current;
                            iMin = iii;
                        }
                    }
                    tg2Current = tG1[i](colIndex) + exp(it.value()/EPSILON) / tSumOverColC(colIndex);

                    if (tg2Current < vMin) {
                        vMin = tg2Current;
                        iMin = colIndex;
                    }
                    anchor = colIndex+1;
                }

                for (iii = anchor; iii < n; ++iii) {
                    tg2Current = tG1[i](iii) +1. / tSumOverColC(iii);
                    if (tg2Current < vMin) {
                        vMin = tg2Current;
                        iMin = iii;
                    }
                }

                mexPrintf("iMin of column %d: %d\n", k, iMin);
                /* ENDING Computing G2 
                 * This is equivalent to update column k of CNew
                 * */

                /* START Update column ck */
                tVecTmp = tCT.row(k);
                // mexPrintf("sum of tVecTmp: %f \n", tVecTmp.sum());
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
            }
            // mexPrintf("Stepsize: %f \n", stepSize);
            // mexPrintf("Len of tripletList: %d \n", tripletList.size());
            /* END Update column ck */
        }
        tCT.setZero();
        tCT.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    mexPrintf("Matrix C is: \n");
    for (i = 0; i < n; ++i) {
        for (j = 0; j < n; ++j) {
            mexPrintf("%f \t", tCT.coeffRef(i, j));
        }
        mexPrintf("\n");
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *X;               /* MxN input matrix */
    double *C;               /* MxN input matrix */
    size_t ncols;                   /* size of matrix */
    size_t nrows;                   /* size of matrix */
    double *Z;              /* output matrix */
    double *sr;
    mwIndex *irs,*jcs,j,k;


    /* create a pointer to the real data in the input matrix  */
    #if MX_HAS_INTERLEAVED_COMPLEX
    X = mxGetDoubles(prhs[0]);
    #else
    X = mxGetPr(prhs[0]);
    #endif

    C = mxGetPr(prhs[1]);
    /* Check data type dimensions of the matrix input y. */
    if (!(mxIsDouble(prhs[1]))){
        mexErrMsgIdAndTxt( "MATLAB:DKU", "Input argument must be of type double.");
    }
    if (mxGetNumberOfDimensions(prhs[1]) != 2){
        mexErrMsgIdAndTxt( "MATLAB:DKU", "Input argument must be two dimensional\n");
    }
    /* Get the size and pointers to input data */
    nrows =mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    // Verify the matrix y is infact sparse
    if (!(mxIsSparse(prhs[1]))){
        mexErrMsgIdAndTxt( "MATLAB: DKU", "Input matrix is not sparse\n");
    } else {    
        sr  = mxGetPr(prhs[1]);
        jcs = mxGetJc(prhs[1]);    
        irs = mxGetIr(prhs[1]);
    }

    /* create a pointer to the real data in the input matrix  */
    // #if MX_HAS_INTERLEAVED_COMPLEX
    // y = mxGetDoubles(prhs[1]);
    // #else
    // y = mxGetPr(prhs[1]);
    // #endif

    /* get dimensions of the input matrix */
    // nrows = mxGetM(prhs[0]);
    // ncols = mxGetN(prhs[0]);
    // mexPrintf("nrows = %u \n", nrows);
    // mexPrintf("ncols = %u \n", ncols);

    /* create the output matrix */
    plhs[0] = mxCreateSparse((mwSize)nrows,(mwSize)ncols,ESTIMATED_NNZ*ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    Z = mxGetDoubles(plhs[0]);
    #else
    Z = mxGetPr(plhs[0]);
    #endif

    mexPrintf("Running with %d threads \n", Eigen::nbThreads());
    /* call the computational routine */
    fwMain(X,sr,jcs,irs,nrows,ncols,Z, 2);
}
