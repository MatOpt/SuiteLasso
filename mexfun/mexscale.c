/***********************************************************************
* mexscale: compute
* 
* mex -O -largeArrayDims mexscale.c
* [sigma,bscale2,cscale2,sbc,sboc,bscale,cscale] 
* = mexscale(sigma,normx,normuxi,bscale,cscale)
************************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

#if !defined(max)
#define max(a,b) (a>b?a:b)
#endif

/**********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )

{    double   sigma, normx, normuxi;
     double   bscale2, cscale2, cst = 1.0;
     double   bscale, cscale, sbc, sboc;
     int      i;
/* CHECK THE DIMENSIONS */

   if (nrhs !=5) {
       mexErrMsgTxt(" mexscale_FL: must have 5 inputs"); }

   /***** assign pointers *****/
       
       sigma = (double) *mxGetPr(prhs[0]);
       normx = (double) *mxGetPr(prhs[1]);
       normuxi = (double) *mxGetPr(prhs[2]);
       bscale = (double) *mxGetPr(prhs[3]);
       cscale = (double) *mxGetPr(prhs[4]);
       
   /****** main body ****************/ 
       if (normx < 1e-7) { normx = 1; normuxi = 1; }
       bscale2 = normx*cst;
       cscale2 = normuxi*cst;
       sbc  = sqrt(bscale2*cscale2);       
       sboc = sqrt(bscale2/cscale2);
       sigma = sigma*(cscale2/bscale2);
       cscale = cscale*cscale2;
       bscale = bscale*bscale2;
       
       for(i=0; i<=7; i++){
       plhs[i] = mxCreateDoubleMatrix(1,1,mxREAL);}
       
       *mxGetPr(plhs[0]) = sigma;
       *mxGetPr(plhs[1]) = bscale2;
       *mxGetPr(plhs[2]) = cscale2;
       *mxGetPr(plhs[3]) = sbc;
       *mxGetPr(plhs[4]) = sboc;
       *mxGetPr(plhs[5]) = bscale;
       *mxGetPr(plhs[6]) = cscale;
 return;
}
/**********************************************************/

