#include "mex.h"
#include <stdio.h>
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))

void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    double* u_forward = mxGetPr(prhs[0]);
    double* v_forward = mxGetPr(prhs[1]);
    double* u_backward = mxGetPr(prhs[2]);
    double* v_backward = mxGetPr(prhs[3]);   
  
    const double tau = mxGetScalar(prhs[4]);       
    const int num_frame = mxGetScalar(prhs[5]);
    
    printf("tau:%f, num_frame:%d\n", tau, num_frame);
    
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
        
    
    float one_over_num_frame = 1./(2*num_frame + 1);
    double* tau_f_over_num_frame = new double[num_frame + 1];
    for(int f = 1; f <=num_frame; f++)
    {
        tau_f_over_num_frame[f] = tau*f/(double)(num_frame);
    }
    
    
    int i, j, ii, jj;
    int idx, tmp_idx_col, tmp_idx_row;    
    int tmp_x, tmp_y;    
    int dx, dy;
    
    const int n_elems = M*N*(2*num_frame + 1);
//     const int n_elems = M*N;
    out[0] = mxCreateDoubleMatrix(n_elems, 1, mxREAL);
    out[1] = mxCreateDoubleMatrix(n_elems, 1, mxREAL); 
    out[2] = mxCreateDoubleMatrix(n_elems, 1, mxREAL); 
    
    int count = 0;
    double* out_i = mxGetPr(out[0]);
    double* out_j = mxGetPr(out[1]);
    double* out_k = mxGetPr(out[2]);
       
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            idx = j*M + i;
               
            tmp_idx_row = i*N + j+1;
            tmp_idx_col = i*N + j+1;
            out_i[count] = tmp_idx_row;
            out_j[count] = tmp_idx_col;
            
            double tmp_dist_forward = sqrt(u_forward[idx]*u_forward[idx] + v_forward[idx]*v_forward[idx]);
            double tmp_dist_backward = sqrt(u_backward[idx]*u_backward[idx] + v_backward[idx]*v_backward[idx]);
            double tmp_coeff_of_zero = 0.5*(1./(1+tmp_dist_forward) + 1./(1+tmp_dist_backward));
            out_k[count] = tmp_coeff_of_zero;
            count++;
            
            for(int f = 1; f <=num_frame; f++)
            {
                int dx = round(u_forward[idx]*tau_f_over_num_frame[f]);
                int dy = round(v_forward[idx]*tau_f_over_num_frame[f]);
                
                tmp_x = max(0, min(N-1, j - dx));
                tmp_y = max(0, min(M-1, i - dy));
                tmp_idx_col = (tmp_y)*N + tmp_x+1; 
                
                out_i[count] = tmp_idx_row;
                out_j[count] = tmp_idx_col;
                if(tmp_x == 0 && tmp_y == 0)
                    out_k[count] = 0;
                else
                    out_k[count] = tmp_coeff_of_zero + one_over_num_frame;
                
                count++;
                
                dx = round(u_backward[idx]*tau_f_over_num_frame[f]);
                dy = round(v_backward[idx]*tau_f_over_num_frame[f]);
                
                tmp_x = max(0, min(N-1, j - dx));
                tmp_y = max(0, min(M-1, i - dy));
                tmp_idx_col = (tmp_y)*N + tmp_x+1; 
                
                out_i[count] = tmp_idx_row;
                out_j[count] = tmp_idx_col;
                if(tmp_x == 0 && tmp_y == 0)
                    out_k[count] = 0;
                else
                    out_k[count] = tmp_coeff_of_zero + one_over_num_frame;                
                count++;
            }
        }
    }
    delete[] tau_f_over_num_frame;
    
    
}