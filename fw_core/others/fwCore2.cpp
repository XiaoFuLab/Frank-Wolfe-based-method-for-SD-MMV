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

const size_t ESTIMATED_NNZ = 500;
const double EPSILON = 0.01;
typedef Eigen::Triplet<double> T;


template<typename T> struct TD;
// void eigenArrayToC(Eigen::MatrixXd *X, double *Z)

void computeSumC(Eigen::SparseMatrix<double> tC, Eigen::VectorXd& tSumRowC) {
    size_t rowIndex;
    size_t anchor=0;
    size_t i=0;
    size_t L = tC.outerSize();
	for (size_t k=0; k<L; ++k)
    {
        anchor=0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(tC,k); it; ++it)
		{
            rowIndex = it.row();

            // exp(0) = 1
            for (i = anchor; i < rowIndex; ++i) {
                tSumRowC(i) += 1.;
            }

            // compute exp(.)
            tSumRowC(rowIndex) += exp(it.value()/EPSILON);

            anchor = rowIndex+1;

            // mexPrintf("[%d, %d] = %f \n", rowIndex, it.col(), it.value()/EPSILON);
            // mexPrintf("Adding %f to cell(%d, %d) \n", exp(it.value()/EPSILON), rowIndex, it.col());
		}
        // Finish remaining element of current columns
        for (i = anchor; i < L; ++i) {
            tSumRowC(i) += 1.;
        }
    }
}


/* The computational routine */
void fwMain(double *X, double *sr, mwIndex *jcs, mwIndex *irs, mwSize m, mwSize n, double* Z)
{
    mwSize i,j;
    size_t row, col;
    double val;

    mexPrintf("m=%d n=%d \n", m, n);
    Eigen::Map<Eigen::MatrixXd> tX(X,m,n);
    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tXT = tX.transpose();

    Eigen::SparseMatrix<double> tZ(m,n);
    // Eigen::SparseVector<double> tVecTmp(n);

    std::vector<T> tripletList;
    tripletList.reserve(ESTIMATED_NNZ*n);
    mexPrintf("ESTIMATED_NNZ=%u and n=%u \n", ESTIMATED_NNZ, n);


    /* Construct Sparse matrix in Eigen */
    for (i = 1; i < n+1; ++i) {
        for (j = (size_t)jcs[i-1]; j < size_t(jcs[i]); ++j) {
            row = irs[j];
            col = i-1;
            val = sr[j];
            tripletList.push_back(T((size_t)row, (size_t)col, val));
        }
    }
    /* End construction */

    Eigen::SparseMatrix<double, Eigen::ColMajor> tC(n,n);
    Eigen::Vector<double, Eigen::RowMajor> tCRow(n);
    // mexPrintf("Len : %d \n", tripletList.size());
    tC.setFromTriplets(tripletList.begin(), tripletList.end());
    tCRow = tC.row(10);
    // mexPrintf("Sum of C: %f \n", tC.sum()); 
    
    Eigen::VectorXd tSumRowC(n); // todo: make sure this is zero vector
    tSumRowC.setZero();
    // for (i = 0; i < n; ++i) {
    //     tSumRowC(i) = 0.0;
    // }
    
    double stepSize = 0;
    Eigen::VectorXd tG1(n);
    Eigen::VectorXd tG2(n);
    Eigen::VectorXd tG(n);
    size_t rowIndex, k, iii, anchor, iMin;
    bool foundIt;
    mexPrintf("Start compute here\n");

    // Eigen::MatrixXd A = Eigen::MatrixXd::Ones(n,n);
    // Eigen::MatrixXd b = Eigen::MatrixXd::Ones(4,n);
    // Eigen::MatrixXd c = Eigen::MatrixXd::Ones(4,n);
    //
    // for (size_t iter = 0; iter <= 10000; iter=iter+4) {
    //     c = b*A; 
    // }
    
    Eigen::MatrixXd A = Eigen::MatrixXd::Ones(n,n);
    Eigen::MatrixXd b = Eigen::MatrixXd::Ones(1,n);
    Eigen::MatrixXd c = Eigen::MatrixXd::Ones(1,n);

    for (size_t iter = 0; iter <= 10000; iter=iter+1) {
        c = b*A; 
    }

    // Eigen::Matrix<double,Eigen::Dynamic, 1, Eigen::ColMajor> tmp1(n);
    // Eigen::Matrix<double,Eigen::Dynamic, 1, Eigen::ColMajor> tmp2(n);
    // tmp1 = tX.col(i);
    // for (size_t iter = 1; iter <= 1000; ++iter) {
    //     tmp2 = tXT*tmp1; 
    // }

    // mexPrintf("Matrix C is: \n");
    // for (i = 0; i < n; ++i) {
    //     for (j = 0; j < n; ++j) {
    //         mexPrintf("%f \t", tC.coeffRef(i, j));
    //     }
    //     mexPrintf("\n");
    // }


    // for (i = 0; i < m; ++i) {
    //     for (j = 0; j < n; ++j) {
    //
    //     }
    //     z[i] = tz(i);
    // }
    // z[0] = m(1, 1);
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
    fwMain(X,sr,jcs,irs,nrows,ncols,Z);
}
