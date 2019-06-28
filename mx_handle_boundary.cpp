#include "mex.h"
#include <math.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

#define round(x) (x >= 0 ? int(x + 0.5) : int(x-0.5))

void mexFunction(int nlhs,mxArray *out[],int nrhs,const mxArray *prhs[])
{
    double* in = mxGetPr(prhs[0]);
    double* m_boundary = mxGetPr(prhs[1]);    
    
    const int M =mxGetM(prhs[0]);
    const int N =mxGetN(prhs[0]);
//     printf("%d %d\n", M, N);
    
    out[0] = mxCreateDoubleMatrix(M, N, mxREAL);    
    double* boundary_handled = (double*)mxGetPr(out[0]);
    
    int i, j;
    for( i = 0; i < M; i++)
    {
        for(j = 0; j < N; j++)
        {
//             printf("%d %d\n", j, i);
            int idx = j*M + i;
            if(m_boundary[idx] == 0)
                boundary_handled[idx] = 0;
            else
                boundary_handled[idx] = in[idx];            
        }
    }
    
}