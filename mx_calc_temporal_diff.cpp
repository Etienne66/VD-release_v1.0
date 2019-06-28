#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))

void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{    
    double* L_reference = mxGetPr(prhs[0]);
    double* L_target = mxGetPr(prhs[1]);
    double* u = mxGetPr(prhs[2]);
    double* v = mxGetPr(prhs[3]);
     
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
//     printf("%d %d\n", M, N);
    
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    
    double* L_diff = (double*)mxGetPr(out[0]);
    int i, j;
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int idx = j*M + i;
            int tmp_x = max(0, min(N-1, j + round(u[idx])));
            int tmp_y = max(0, min(M-1, i + round(v[idx])));

            int tmp_idx = tmp_x*M + tmp_y;
            L_diff[idx] = L_reference[idx] - L_target[tmp_idx];
        }
    }
}