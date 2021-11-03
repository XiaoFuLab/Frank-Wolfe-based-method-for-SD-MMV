#include <omp.h>
#include "mex.h"
#include "matrix.h"

#define ll long long

void omp_smm(double* A, double*B, double* C, ll m, ll p, ll n, ll* irs, ll* jcs)
{
    for (ll j=0; j<p; ++j)
    {
        ll istart = jcs[j];
        ll iend = jcs[j+1];
        #pragma omp parallel for
        for (ll ii=istart; ii<iend; ++ii)
        {
            ll i = irs[ii];
            double aa = A[ii];
            for (ll k=0; k<n; ++k)
            {
                C[i+k*m] += B[j+k*p]*aa;
            }
        }
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C; /* pointers to input & output matrices*/
    size_t m,n,p;      /* matrix dimensions */

    A = mxGetPr(prhs[0]); /* first sparse matrix */
    B = mxGetPr(prhs[1]); /* second full matrix */

    mwIndex * irs = mxGetIr(prhs[0]);
    mwIndex * jcs = mxGetJc(prhs[0]);

    m = mxGetM(prhs[0]);  
    p = mxGetN(prhs[0]);
    n = mxGetN(prhs[1]);

    /* create output matrix C */
    plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
    C = mxGetPr(plhs[0]);

    omp_smm(A,B,C, m, p, n, (ll*)irs, (ll*)jcs);
}
