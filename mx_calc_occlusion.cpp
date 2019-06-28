#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))
double calc_dist(double* ref_patch, double* qry_patch, int dims);



void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    double* u_forward = mxGetPr(prhs[0]);
    double* v_forward = mxGetPr(prhs[1]);
    double* u_backward = mxGetPr(prhs[2]);
    double* v_backward = mxGetPr(prhs[3]);
    
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
    const int MN = M*N;
    
    int i, j, k;    
    
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    double* out_occlusion = mxGetPr(out[0]); 
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int idx = j*M + i;
            int tmp_u_forward = max(0, min(N - 1, j + round(u_forward[idx])));
            int tmp_v_forward = max(0, min(M - 1, i + round(v_forward[idx])));
            int tmp_idx = tmp_u_forward*M + tmp_v_forward;
            
            int tmp_u_backward = max(0, min(N - 1, tmp_u_forward + round(u_backward[tmp_idx])));
            int tmp_v_backward = max(0, min(M - 1, tmp_v_forward + round(v_backward[tmp_idx])));
//             
//             if(idx != (tmp_u_backward*M + tmp_v_backward))            

            
            if(tmp_u_backward == j &&  tmp_v_backward == i)                 
                out_occlusion[idx] = 0;
            else
                out_occlusion[idx] = 1;                
        }
    }    
}