#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define NUM_FRAME 70
// #define ONE_OVER_NUM_FRAME 1./(2*70+1)

// #define NUM_FRAME 140
// #define ONE_OVER_NUM_FRAME 1./140
#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))

void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    int i, j, k;
    
    
    double* L = mxGetPr(prhs[0]);
    double* L_dx = mxGetPr(prhs[1]);
    double* L_dy = mxGetPr(prhs[2]);
    double* u_forward = mxGetPr(prhs[3]);
    double* v_forward = mxGetPr(prhs[4]);
    double* u_backward = mxGetPr(prhs[5]);
    double* v_backward = mxGetPr(prhs[6]);    
    double tau = mxGetScalar(prhs[7]);
    
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
    const int MN = M*N;
    
    double* KL = new double[MN];
    double* KL_dx = new double[MN];
    double* KL_dy = new double[MN];
    
  

    float one_over_num_frame = 1./(2*NUM_FRAME + 1);
    double* tau_f_over_num_frame = new double[NUM_FRAME + 1];
    for(int f = 1; f <=NUM_FRAME; f++)
    {
        tau_f_over_num_frame[f] = tau*f/(double)(NUM_FRAME);
    }

    
    int idx, tmp_idx;
    int tmp_x, tmp_y;
    int dx, dy;
    
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            idx = j*M + i;
                      
            KL[idx] = one_over_num_frame*L[idx];
            KL_dx[idx] = one_over_num_frame*L_dx[idx];
            KL_dy[idx] = one_over_num_frame*L_dy[idx];
            for(int f = 1; f <=NUM_FRAME; f++)
            {  
                int dx = round(u_forward[idx]*tau_f_over_num_frame[f]);
                int dy = round(v_forward[idx]*tau_f_over_num_frame[f]);
                
                tmp_x = max(0, min(N-1, j - dx));
                tmp_y = max(0, min(M-1, i - dy));
                tmp_idx = tmp_x*M + tmp_y;                 
                KL[idx] += one_over_num_frame*L[tmp_idx];
                KL_dx[idx] += one_over_num_frame*L_dx[tmp_idx];
                KL_dy[idx] += one_over_num_frame*L_dy[tmp_idx];
             
                dx = round(u_backward[idx]*tau_f_over_num_frame[f]);
                dy = round(v_backward[idx]*tau_f_over_num_frame[f]);
                
                tmp_x = max(0, min(N-1, j - dx));
                tmp_y = max(0, min(M-1, i - dy));
                tmp_idx = tmp_x*M + tmp_y; 
                KL[idx] += one_over_num_frame*L[tmp_idx];
                KL_dx[idx] += one_over_num_frame*L_dx[tmp_idx];
                KL_dy[idx] += one_over_num_frame*L_dy[tmp_idx];
            }
        }
    }
    delete[] tau_f_over_num_frame;
    
    
    
    
    
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    out[1] = mxCreateDoubleMatrix(M, N, mxREAL);
    out[2] = mxCreateDoubleMatrix(M, N, mxREAL);
    double* out_KL = (double*)mxGetPr(out[0]);
    double* out_KL_dx = (double*)mxGetPr(out[1]);
    double* out_KL_dy = (double*)mxGetPr(out[2]);
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int idx = j*M + i;
            out_KL[idx] = KL[idx];
            out_KL_dx[idx] = KL_dx[idx];
            out_KL_dy[idx] = KL_dy[idx];
        }
    }
    
    
    delete []KL;
    delete []KL_dx;
    delete []KL_dy;
}