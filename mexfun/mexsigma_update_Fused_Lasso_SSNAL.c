/***********************************************************************
* mexsigma_update: compute
* 
* mex -O -largeArrayDims mexsigma_update_Fused_Lasso_SSNAL.c
* 
************************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

#if !defined(max)
#define max(a,b) (a>b?a:b)
#endif

#if !defined(min)
#define min(a,b) (a<b?a:b)
#endif
/**********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[], 
                 int nrhs,   const mxArray  *prhs[] )

{    double   sigma, sigmamax, sigmamin, prim_win, dual_win;
     int      sigma_update_iter, iter, innerop, i;
     double  sigmascale =  5.0;/*10;*/
/* CHECK THE DIMENSIONS */

   if (nrhs <7) {
       mexErrMsgTxt(" mexsigma_update: must have 7 inputs"); }

   /***** assign pointers *****/
       
       sigma = (double) *mxGetPr(prhs[0]);
       sigmamax = (double) *mxGetPr(prhs[1]);
       sigmamin = (double) *mxGetPr(prhs[2]);
       prim_win = (double) *mxGetPr(prhs[3]);
       dual_win = (double) *mxGetPr(prhs[4]);
       iter = (int) *mxGetPr(prhs[5]);
       innerop = (int) *mxGetPr(prhs[6]);
       
   /****** main body ****************/ 

       if ( iter < 10)
           sigma_update_iter = 1;
       else if(iter < 20)
           sigma_update_iter = 2;
       else if(iter < 200)
           sigma_update_iter = 2;
       else if(iter < 500)
           sigma_update_iter = 3;
       else
           sigma_update_iter = 3;
       
       if (innerop == 1) {
          sigmascale = sqrt(5.0);
       }
              
       if ((iter%sigma_update_iter == 0)) {
          if (prim_win > max(1, 1.2*dual_win)){
             prim_win = 0.0;
             sigma = min(sigmamax,sigma*sigmascale);
          }else if (dual_win > max(1, 3*prim_win)){
             dual_win = 0.0;
             sigma = max(sigmamin,2*sigma/sigmascale);
          } 
       }
              
       for(i=0; i<= 2; i++){
       plhs[i] = mxCreateDoubleMatrix(1,1,mxREAL);}
       
       *mxGetPr(plhs[0]) = sigma;
       *mxGetPr(plhs[1]) = prim_win;
       *mxGetPr(plhs[2]) = dual_win;
 return;
}
/**********************************************************/

