#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))

void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{    
    double* dual_vars = mxGetPr(prhs[0]);
    double* u = mxGetPr(prhs[1]);
    double* v = mxGetPr(prhs[2]);
     
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
//     printf("%d %d\n", M, N);
    
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    out[1] = mxCreateDoubleMatrix(M, N, mxREAL);
    
    double* out1 = (double*)mxGetPr(out[0]);
    double* out2 = (double*)mxGetPr(out[1]);
    
    int i, j;
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int idx = j*M + i;
            out1[idx] = 0;
            out2[idx] = 0;
        }
    }
    
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int src_idx = j*M + i;
            int tmp_x = max(0, min(N-1, j + round(u[src_idx])));
            int tmp_y = max(0, min(M-1, i + round(v[src_idx])));

            int target_idx = tmp_x*M + tmp_y;
//             L_diff[idx] = dual_vars[idx] - L_target[tmp_idx];
//             L_diff[idx] = L_target[idx] - dual_vars[tmp_idx];
            out1[src_idx] += dual_vars[src_idx];
            out2[target_idx] -= dual_vars[src_idx];
        }
    }
}