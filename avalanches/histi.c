#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int32_T *X;                                   	/*input*/
    double *C, *F, *xlim;                           /*output*/
    int i, n, xmin, xmax, d;
    
    /* check number of variables */
    if(nrhs!=1)
        mexErrMsgTxt("There is one input variable.");
    if(nlhs>3)
        mexErrMsgTxt("There are three output variables.");
    
    /* initialize and verify input */
    if(X = (int32_T *)mxGetData(prhs[0]), !mxIsInt32(prhs[0]))
        mexErrMsgTxt("X must be in int32 format.");
    if(!(n = mxGetNumberOfElements(prhs[0])))
        mexErrMsgTxt("X must contain at least one element.");
    
    /* get xmin and xmax */
    for(xmin=xmax=*X, i=1; i<n; i++)
        if(*(X+i)>xmax)
            xmax=*(X+i);
        else if(*(X+i)<xmin)
            xmin=*(X+i);
    d=xmax-xmin;
    
    /* inialize output */
    plhs[0] = mxCreateDoubleMatrix(1,d+1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,d+1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,2,mxREAL);
    F = mxGetPr(plhs[0]);                           /* frequency distribution */
    C = mxGetPr(plhs[1]);                           /* cumulative frequency distribution */
    xlim = mxGetPr(plhs[2]);                        /* xmin and xmax */
    
    /* get frequency distribution */
    for(i=0; i<n; i++)
        *(F+*(X+i)-xmin)+=1;
    
    /* get cumulative frequency distribution */
    *(C+d)=*(F+d);
    for(i=d-1; i>=0; i--)
        *(C+i)=*(C+i+1)+*(F+i);
    
    /* get xlim */
    *(xlim+0)=xmin;
    *(xlim+1)=xmax;
}
