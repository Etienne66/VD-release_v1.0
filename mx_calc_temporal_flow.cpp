#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))

void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    
    
    double* u_reference = mxGetPr(prhs[0]);
    double* v_reference = mxGetPr(prhs[1]);
    double* u_target = mxGetPr(prhs[2]);
    double* v_target = mxGetPr(prhs[3]);
    double* m = mxGetPr(prhs[4]);
    
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
    
    int i, j, idx;
    int tmp_x, tmp_y, tmp_idx;
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);//u_new
    out[1] = mxCreateDoubleMatrix(M, N, mxREAL);//v_new
    out[2] = mxCreateDoubleMatrix(M, N, mxREAL);//boundary
    
    double* u = (double*)mxGetPr(out[0]);
    double* v = (double*)mxGetPr(out[1]);
    double* is_valid = (double*)mxGetPr(out[2]);
    
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            idx = j*M + i;            
            u[idx] = round(u_reference[idx]);
            v[idx] = round(v_reference[idx]);
            
            tmp_x = j + u[idx];
            tmp_y = i + v[idx];                
            if (tmp_x < 0 || tmp_x > N - 1 || tmp_y < 0 || tmp_y > M - 1)
            {                
                is_valid[idx] = 0;
                
                tmp_idx = max(0, min(N-1, tmp_x))*M + max(0, min(M-1, tmp_y));
                u[idx] += u_target[tmp_idx];
                v[idx] += v_target[tmp_idx];
            }
            else
            {
                is_valid[idx] = m[idx];
                
                tmp_idx = tmp_x*M + tmp_y;
                u[idx] += u_target[tmp_idx];
                v[idx] += v_target[tmp_idx];
            }            
        }
    }
}