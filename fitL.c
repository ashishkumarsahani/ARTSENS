/*==========================================================
 * arrayProduct.c - example in MATLAB External Interfaces
 *
 * Multiplies an input scalar (multiplier) 
 * times a 1xN matrix (inMatrix)
 * and outputs a 1xN matrix (outMatrix)
 *
 * The calling syntax is:
 *
 *		outMatrix = arrayProduct(multiplier, inMatrix)
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2007-2008 The MathWorks, Inc.
 *
 *========================================================*/
/* $Revision: 1.1.10.2 $ */

#include "mex.h"
#include <math.h>
#include "stdlib.h"
#define round(x) floor(x + 0.5);

void fitL(double *well_Half, int Tot_Len, double *inflection_Ptr, double *error_Ptr, double *curve_Ptr, double *trial_Count)
{
    
    int Jump = 1;
    int x = 0;
    int xb, xm, xe;
    double yb,ym,ye;
    double *L1;
    double *L2;
    double Err1=0.0,Err2=0.0;
    int i=0;
    *error_Ptr = 100.0;
    *trial_Count = 0.0;
    L1 = (double *)malloc(Tot_Len*sizeof(double));
    L2 = (double *)malloc(Tot_Len*sizeof(double));
    x = (int)round(2*Tot_Len/3);
   
    
    while(1)
    {
        *trial_Count = *trial_Count +1.0;
        x = x- Jump;
        if(x<=1  || x>Tot_Len-1) break;     
        xb =1; yb = well_Half[0]; xm = x; ym = well_Half[x-1]; ye = well_Half[Tot_Len-1]; xe = Tot_Len;
        for(i = 1;i<= xm;i++)
            L1[i-1] = ((ym-yb)/(xm-xb))*(i-xb)+yb;
        for(i = 1;i<=xe-xm;i++)
            L1[xm+i-1] = ((ye-ym)/(xe-xm))*(i)+ym;

        Err1=0.0;
        for(i = 0;i< xe;i++)
            Err1 = Err1+(L1[i]-well_Half[i])*(L1[i]-well_Half[i]);

        if (Err1 < *error_Ptr)
        {
            *inflection_Ptr = (double)x;
            memcpy(curve_Ptr, L1, Tot_Len*sizeof(double));
            *error_Ptr = Err1;
            if(*error_Ptr < 0.1){
                break;
            }
        }
        xb =1; yb = well_Half[0]; xm = x+1; ym = well_Half[x+1-1]; ye = well_Half[Tot_Len-1]; xe = Tot_Len;
        for(i = 1;i<= xm;i++)
            L2[i-1] = ((ym-yb)/(xm-xb))*(i-xb)+yb;
        for(i = 1;i<=xe-xm;i++)
            L2[xm+i-1] = ((ye-ym)/(xe-xm))*(i)+ym;
        
        Err2=0.0;
        for(i = 0;i< xe;i++)
            Err2 = Err2+(L2[i]-well_Half[i])*(L2[i]-well_Half[i]);
        
        if(*trial_Count>=20.0) {
            break;
        }
        if(Err2 ==Err1) continue;
        Jump = (int)round(Err1/(Err2-Err1));
    }
    free(L1);
    free(L2);
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double length;
    double *well_Half;
    mwSize nrows;
    double *inflection_Ptr;
    double *error_Ptr;
    double *curve_Ptr;
    double *trial_Count;

    double pointer_Init_Val = 0;
    int i=0;
    /* check for proper number of arguments */
    if(nrhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","One input required - Well half");
    }
    if(nlhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","Four output required - Inflection point, Error, Curve, trial_Count");
    }
    
    /* check that number of column in first input argument is 1 */
    if(mxGetN(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notColumnVector","Well half must be a column vector.");
    }

    /* create a pointer to the real data in the input matrix  */
    well_Half = mxGetPr(prhs[0]);

    /* get dimensions of the input matrix */
    length = mxGetM(prhs[0]);
    //mexPrintf("%u\n", nrows);
    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix((unsigned int)length,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);
    
    inflection_Ptr = mxGetPr(plhs[0]);
    error_Ptr = mxGetPr(plhs[1]);
    curve_Ptr = mxGetPr(plhs[2]);
    trial_Count =  mxGetPr(plhs[3]);
    fitL(well_Half, (int)length, inflection_Ptr, error_Ptr, curve_Ptr, trial_Count);
}
