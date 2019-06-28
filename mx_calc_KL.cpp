#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))
#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))


#if 0


void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
//     int *data = (int *) mxGetData(prhs[0]);
    
//     double* K_i = mxGetPr(prhs[0]);
//     double* K_j = mxGetPr(prhs[1]);
     int* K_i = (int *) mxGetData(prhs[0]);
    int* K_j = (int *) mxGetData(prhs[1]);
    double K_s = mxGetScalar(prhs[2]);
    double* L = mxGetPr(prhs[3]);
    
    const int n_elems = mxGetM(prhs[0])*mxGetN(prhs[0]);
    
    const int M = mxGetM(prhs[3]);
    const int N = mxGetN(prhs[3]);
    
    printf("n_elems:%d, %d %d\n", n_elems, M, N);
    
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    double* KL = (double*)mxGetPr(out[0]);
//     double* KL = new double[MN]();
    for(int k = 0; k < n_elems; k++)
    {
        KL[K_i[k]] += L[K_j[k]];
//         KL[(int)K_i[k]] += L[(int)K_j[k]];
//         KL[(int)K_i[k]] += L[(int)K_j[k]];
//         KL[(int)K_i[k]] += L[(int)K_j[k]];
//         KL[(int)K_i[k]] += K_s*L[(int)K_j[k]];
//         KL[(int)K_i[k]] += K_s*L[(int)K_j[k]];
    }
        
    
    
//     for(i = 0; i < M; i++)
//         {
//             for(j = 0; j < N; j++)
//             {
//                 out_rho[i*N + j] =  L[i*N + j];
//             }
// 
//         }
    
}




#else



#define NUM_FRAME 70
void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    int i, j, k;
    
    
    double* L = mxGetPr(prhs[0]);
    double* u_forward = mxGetPr(prhs[1]);
    double* v_forward = mxGetPr(prhs[2]);
    double* u_backward = mxGetPr(prhs[3]);
    double* v_backward = mxGetPr(prhs[4]);    
    double tau = mxGetScalar(prhs[5]);
    
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
    const int MN = M*N;
    
    double* KL = new double[MN];      

    float one_over_num_frame = 1./(2*NUM_FRAME + 1);
    double* tau_f_over_num_frame = new double[NUM_FRAME + 1];
    for(int f = 1; f <=NUM_FRAME; f++)
    {
        tau_f_over_num_frame[f] = tau*f/(double)(NUM_FRAME);
    }

    
    int idx, tmp_idx;
    int tmp_x, tmp_y;
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            idx = j*M + i;                      
            KL[idx] = 0;
        }
    }
    
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            idx = j*M + i;
                      
            KL[idx] = KL[idx] + one_over_num_frame*L[idx];
            for(int f = 1; f <=NUM_FRAME; f++)
            {
//                 tmp_x = max(0, min(N-1, 0.5 + j + u_forward[idx]*tau_f_over_num_frame[f]));
//                 tmp_y = max(0, min(M-1, 0.5 + i + v_forward[idx]*tau_f_over_num_frame[f]));
//                 tmp_x = max(0, min(N-1, round(j + u_forward*tau_f_over_num_frame[f])));
//                 tmp_y = max(0, min(M-1, round(i + v_forward*tau_f_over_num_frame[f])));
//                 tmp_idx = tmp_x*M + tmp_y;
                
                int dx = round(u_forward[idx]*tau_f_over_num_frame[f]);
                int dy = round(v_forward[idx]*tau_f_over_num_frame[f]);
                
                tmp_x = max(0, min(N-1, j - dx));
                tmp_y = max(0, min(M-1, i - dy));
                tmp_idx = tmp_x*M + tmp_y; 
                KL[idx] += one_over_num_frame*L[tmp_idx];
//                 KL[tmp_idx] += one_over_num_frame*L[idx];
                
//                 tmp_x = max(0, min(N-1, 0.5 + j + u_backward[idx]*tau_f_over_num_frame[f]));
//                 tmp_y = max(0, min(M-1, 0.5 + i + v_backward[idx]*tau_f_over_num_frame[f]));
//                 tmp_x = max(0, min(N-1, round(j + u_backward*tau_f_over_num_frame[f])));
//                 tmp_y = max(0, min(M-1, round(i + v_backward*tau_f_over_num_frame[f])));
//                 tmp_idx = tmp_x*M + tmp_y;
             
                dx = round(u_backward[idx]*tau_f_over_num_frame[f]);
                dy = round(v_backward[idx]*tau_f_over_num_frame[f]);
                
                tmp_x = max(0, min(N-1, j - dx));
                tmp_y = max(0, min(M-1, i - dy));
                tmp_idx = tmp_x*M + tmp_y; 
                KL[idx] += one_over_num_frame*L[tmp_idx];
//                 KL[tmp_idx] += one_over_num_frame*L[idx];
            }
        }
    }
    delete[] tau_f_over_num_frame;
    
    
    
    
    
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    double* out_rho = (double*)mxGetPr(out[0]);
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int idx = j*M + i;
            out_rho[idx] = KL[idx];
        }
    }
    
    
    delete []KL;
}


#endif