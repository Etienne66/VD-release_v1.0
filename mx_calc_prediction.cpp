#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))

void quantize_4_bins(double* A, int M, int N, int* bins, double bin_max[4]);


void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    int i, j;   
    
    double* Px = mxGetPr(prhs[0]);
    double* Py = mxGetPr(prhs[1]);    
    
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
    const int MN = M*N;   
    
    int* Px_bins = new int[MN]();
    int* Py_bins = new int[MN]();
    
    double Px_bin_max[4] = {0, 0, 0, 0};
    double Py_bin_max[4] = {0, 0, 0, 0};
    
    quantize_4_bins(Px, M, N, Px_bins, Px_bin_max);
    quantize_4_bins(Py, M, N, Py_bins, Py_bin_max);

    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);
    out[1] = mxCreateDoubleMatrix(M, N, mxREAL);
    
    double* out_Px = (double*)mxGetPr(out[0]);
    double* out_Py = (double*)mxGetPr(out[1]);
    
    double thr = 0.05;
    for(i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
            int idx = j*M + i;
            
            if(fabs(Px[idx]) < thr*Px_bin_max[Px_bins[idx]])
                out_Px[idx] = 0;
            else
                out_Px[idx] = Px[idx];
            
            if(fabs(Py[idx]) < thr*Py_bin_max[Py_bins[idx]])
                out_Py[idx] = 0;                
            else
                out_Py[idx] = Py[idx];
        }
    }
    
    delete[] Px_bins;
    delete[] Py_bins;
}


void quantize_4_bins(double* A, int M, int N, int* bins, double bin_max[4])
{    
    int i, j, k;    
    double tmp_magnitude[8];
    
    for(i = 1; i < M-1; i++)
    {
        for(j = 1; j < N-1; j++)
        {
            int idx = j*M + i;
            
            tmp_magnitude[0] = fabs(A[idx] - A[(j+1)*M + (i)]);
            tmp_magnitude[1] = fabs(A[idx] - A[(j+1)*M + (i-1)]);
            tmp_magnitude[2] = fabs(A[idx] - A[(j)*M + (i-1)]);
            tmp_magnitude[3] = fabs(A[idx] - A[(j-1)*M + (i-1)]);
            tmp_magnitude[4] = fabs(A[idx] - A[(j-1)*M + (i)]);
            tmp_magnitude[5] = fabs(A[idx] - A[(j-1)*M + (i+1)]);
            tmp_magnitude[6] = fabs(A[idx] - A[(j)*M + (i+1)]);
            tmp_magnitude[7] = fabs(A[idx] - A[(j+1)*M + (i+1)]);
            
            int tmp_max_val = 0;
            int tmp_max_bin = 0;
            for(k = 0; k < 8; k++)
            {
                if(tmp_magnitude[k] > tmp_max_val)
                {
                    tmp_max_bin = k;
                    tmp_max_val = tmp_magnitude[k];
                }
            }
            
            if(tmp_max_bin == 0 || tmp_max_bin == 4)
            {
                bins[idx] = 0;
            }
            else if(tmp_max_bin == 1 || tmp_max_bin == 5)
            {
                bins[idx] = 1;
            }
            else if(tmp_max_bin == 2 || tmp_max_bin == 6)
            {
                bins[idx] = 2;
            }
            else
            {
                bins[idx] = 3;
            }
            
            if (bin_max[bins[idx]] < fabs(A[idx]))
                bin_max[bins[idx]] = fabs(A[idx]);            
        }
    }
}
