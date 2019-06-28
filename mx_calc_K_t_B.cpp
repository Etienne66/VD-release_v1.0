#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))
#define NUM_FRAME 70
void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    int i, j, k;
    
    
    double* B = mxGetPr(prhs[0]);
    double* u_forward = mxGetPr(prhs[1]);
    double* v_forward = mxGetPr(prhs[2]);
    double* u_backward = mxGetPr(prhs[3]);
    double* v_backward = mxGetPr(prhs[4]);    
    double tau = mxGetScalar(prhs[5]);
    
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
    const int MN = M*N;
    
    double* K_t_B = new double[MN];
    
//     double u_forward = 4.8751;
//     double v_forward = 3.6450;
//     double u_backward = -2.8259;
//     double v_backward = -2;
    

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
            K_t_B[idx] = 0;
        }
    }
    
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            idx = j*M + i;
                      
            K_t_B[idx] += one_over_num_frame*B[idx];
            for(int f = 1; f <=NUM_FRAME; f++)
            {
                dx = round(u_forward[idx]*tau_f_over_num_frame[f]);
                dy = round(v_forward[idx]*tau_f_over_num_frame[f]);
                
                tmp_x = max(0, min(N-1, j - dx));
                tmp_y = max(0, min(M-1, i - dy));
                tmp_idx = tmp_x*M + tmp_y;
                K_t_B[tmp_idx] += one_over_num_frame*B[idx];
             
                dx = round(u_backward[idx]*tau_f_over_num_frame[f]);
                dy = round(v_backward[idx]*tau_f_over_num_frame[f]);
                
                tmp_x = max(0, min(N-1, j - dx));
                tmp_y = max(0, min(M-1, i - dy));
                tmp_idx = tmp_x*M + tmp_y;
                K_t_B[tmp_idx] += one_over_num_frame*B[idx];
            }
        }
    }
    delete[] tau_f_over_num_frame;
    
    
    
    
    
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    double* out_K_t_B = (double*)mxGetPr(out[0]);
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int idx = j*M + i;
            out_K_t_B[idx] = K_t_B[idx];
        }
    }
    
    
    delete []K_t_B;
}