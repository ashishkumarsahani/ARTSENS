
#include "mex.h"
#include "math.h"
#include "stdlib.h"

/* The computational routine */
void peak_detect(double *inData, double delta, double in_Data_Rows, double *peak_Locs, double *peak_Values, double *valley_Locs, double *valley_Values, double *peak_Arr_Len, double *valley_Arr_Len)
{
    double mn, mx, mnpos, mxpos, max_Val, min_Val;
    double this=0.0;
    unsigned int lookformax;
    int i =0;
    //peak_Locs = (double *)mxMalloc((unsigned int)(in_Data_Rows/2)*sizeof(double));
    //valley_Locs = (double *)mxMalloc((unsigned int)(in_Data_Rows/2)*sizeof(double));
    //peak_Values = (double *)mxMalloc((unsigned int)(in_Data_Rows/2)*sizeof(double));
    //valley_Values = (double *)mxMalloc((unsigned int)(in_Data_Rows/2)*sizeof(double));
    *peak_Arr_Len = 0;
    *valley_Arr_Len = 0;
    
    mn = inData[0];
    mx= inData[0];
    mnpos = 0; mxpos = 0;

    lookformax = 1;
    for(i=0; i<in_Data_Rows;i++)
    {
        this = inData[i];
        if (lookformax ==1)
        {
            if (this > mx) 
            {
                mx = this; 
                mxpos =i;
            }
            else if (this < mx-delta)
            {
              peak_Locs[(unsigned int)*peak_Arr_Len] = mxpos+1; //for matlab compatibility added 1
              peak_Values[(unsigned int)*peak_Arr_Len] = mx;
              *peak_Arr_Len = *peak_Arr_Len+1;
              //mexPrintf("%lf %lf %lf\n",(*peak_Arr_Len)-1,peak_Locs[((unsigned int)*peak_Arr_Len)-1], peak_Values[((unsigned int)*peak_Arr_Len)-1]);
              mn = this; mnpos = i;
              lookformax = 0;
            }
        }
        else
        {
            if (this < mn)
            {
                mn = this; 
                mnpos =i;
            }
            else if (this > mn+delta)
            {
              valley_Locs[(unsigned int)*valley_Arr_Len] = mnpos+1; //for matlab compatibility added 1
              valley_Values[(unsigned int)*valley_Arr_Len] = mn;
              *valley_Arr_Len = *valley_Arr_Len+1;
              mx = this; mxpos = i;
              lookformax = 1;
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double delta;
    double *inData;
    mwSize nrows;
    double *valley_Locs_Ptr;
    double *peak_Locs_Ptr;
    double *valley_Vals_Ptr;
    double *peak_Vals_Ptr;
    double *peak_No_Ptr;
    double *valley_No_Ptr;
    double pointer_Init_Val = 0;
    int i=0;
    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Six outputs required.");
    }
    /* make sure the second input argument is scalar */
    if( !((mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1]) ||
         mxGetNumberOfElements(prhs[1])!=1) && (mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2]) ||
         mxGetNumberOfElements(prhs[2])!=1) )) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","delta and length must be a scalar.");
    }
    
    /* check that number of column in first input argument is 1 */
    if(mxGetN(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notColumnVector","Input must be a column vector.");
    }
    
    /* get the value of the scalar input  */
    delta = mxGetScalar(prhs[1]);

    /* create a pointer to the real data in the input matrix  */
    inData = mxGetPr(prhs[0]);

    /* get dimensions of the input matrix */
    nrows = mxGetM(prhs[0]);
    //mexPrintf("%u\n", nrows);
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((unsigned int)(nrows/2),1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((unsigned int)(nrows/2),1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix((unsigned int)(nrows/2),1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix((unsigned int)(nrows/2),1,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,1,mxREAL);
    peak_Locs_Ptr = mxGetPr(plhs[0]);
    peak_Vals_Ptr = mxGetPr(plhs[1]);
    valley_Locs_Ptr = mxGetPr(plhs[2]);
    valley_Vals_Ptr = mxGetPr(plhs[3]);
    peak_No_Ptr = mxGetPr(plhs[4]);
    valley_No_Ptr = mxGetPr(plhs[5]);
    peak_detect(inData, delta, nrows, peak_Locs_Ptr,peak_Vals_Ptr, valley_Locs_Ptr, valley_Vals_Ptr, peak_No_Ptr, valley_No_Ptr);
}
