#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))
double calc_dist(double* ref_patch, double* qry_patch, int dims);



void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    double* ref_img = mxGetPr(prhs[0]);
    double* qry_img = mxGetPr(prhs[1]);
    
    double* u = mxGetPr(prhs[2]);
    double* v = mxGetPr(prhs[3]);    
    const int block_size = mxGetScalar(prhs[4]);//should be odd
    const int window_size = mxGetScalar(prhs[5]);//should be odd
    
    const int half_block_size = int(block_size/2);
    const int half_window_size = int(window_size/2);
            
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
    const int MN = M*N;
    
    int i, j, k;    
    int ii, jj;
    int iii,jjj;
    
    const int block_size_sqr = block_size*block_size;
    const int window_size_sqr = window_size*window_size;
    
    double* ref_patch = new double[block_size_sqr];
    double* qry_patch = new double[block_size_sqr];
    
    out[0] = mxCreateDoubleMatrix(window_size*window_size, M*N, mxREAL);
    out[1] = mxCreateDoubleMatrix(window_size*window_size, M*N, mxREAL);
    double* out_dist = mxGetPr(out[0]); 
    double* out_val = mxGetPr(out[1]); 
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int idx = j*M + i;
            int tmp_u = j + round(u[idx]);
            int tmp_v = i + round(v[idx]);
            
            
            
            int count = 0;
            for(iii = -half_block_size; iii <= half_block_size; iii++)
            {
                for(jjj = -half_block_size; jjj <= half_block_size; jjj++)
                {
                    int tmp_x = max(0, min(N - 1, j + jjj));
                    int tmp_y = max(0, min(M - 1, i + iii));
                    int tmp_idx = tmp_x*M + tmp_y;
                    ref_patch[count++] = ref_img[tmp_idx];
                }                
            }
            
            int count_widnow = 0;
            for(ii = -half_window_size; ii <= half_window_size; ii++)
            {
                for(jj = -half_window_size; jj <= half_window_size; jj++)
                {
                    int count = 0;
                    for(iii = -half_block_size; iii <= half_block_size; iii++)
                    {
                        for(jjj = -half_block_size; jjj <= half_block_size; jjj++)
                        {
                            int tmp_x = max(0, min(N -1, tmp_u + jj + jjj));
                            int tmp_y = max(0, min(M - 1, tmp_v + ii + iii));
                            int tmp_idx = tmp_x*M + tmp_y;
                            qry_patch[count++] = qry_img[tmp_idx];
                        }                
                    }
                    
                    int out_idx = idx*window_size_sqr + count_widnow;
                    out_dist[out_idx] = calc_dist(ref_patch, qry_patch, block_size_sqr);
                    
                    int tmp_x = max(0, min(N -1, tmp_u + jj ));
                    int tmp_y = max(0, min(M - 1, tmp_v + ii ));
                    int tmp_idx = tmp_x*M + tmp_y;
                    out_val[out_idx] = qry_img[tmp_idx];
                    count_widnow++;
                }
            }
        }
    }
    delete[] ref_patch;
    delete[] qry_patch;
    
}

double calc_dist(double* ref_patch, double* qry_patch, int dims)
{
    double sigma = 10./255;
//     double sigma = 7.5/255;
//     double sigma = 8./255;
//     double sigma = 5./255;    
    sigma = (sigma*dims);
    double var = 2*sigma*sigma;
    double ret = 0;
    for(int i = 0; i < dims; i++)
    {
        double dist = (ref_patch[i] - qry_patch[i])*(ref_patch[i] - qry_patch[i]);
        ret += dist;
    }
    ret = ret;
    ret = exp(-ret/var);
    return ret;
}