/**********************************************************
* mexFusedLassoJacobian.c
*            
* [h,U] = mexFusedLassoJacobian(r); 
* Jacobian = spdiags(h,0,nr+1,nr+1)+U*U'; 
*
* mex  -O -largeArrayDims  mexFusedLassoJacobian.c
* Copyright (c) 2016 by
* Xudong Li, Defeng Sun, and Kim-Chuan Toh 
**********************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif
        
#if !defined(MIN)
#define  MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif
        
#if !defined(SQR)
#define SQR(x) ((x)*(x))
#endif

/*#include "mymexheader.h"*/

/********************************************************************
  PROCEDURE mexFunction - Entry for Matlab
*********************************************************************/
void mexFunction(const int nlhs, mxArray *plhs[],
                 const int nrhs, const mxArray *prhs[])
{
  double   *rr, *vv, *hh;
  ptrdiff_t  *blklen, *blkend, *ii, *jj; 
  int      nr, len, numblk, NZ, idxstart, idxend;
  int      j, k, cnt; 
  double   tmp; 

  if(nrhs != 1)
    mexErrMsgTxt("mexFusedLassoJacobian: requires 1 input arguments.");
  if(nlhs != 2)
    mexErrMsgTxt("mexFusedLassoJacobian: requires 2 output argument.");
  if (mxIsLogical(prhs[0])) {
     mexErrMsgTxt("mexFusedLassoJacobian: input cannot be a logical array."); 
  }
  nr = MAX(mxGetM(prhs[0]),mxGetN(prhs[0])); 
  if (mxIsSparse(prhs[0])) {
     mexErrMsgTxt("mexFusedLassoJacobian: input cannot be sparse"); 
  } else {
     rr = mxGetPr(prhs[0]); 
  }
  /*****************************/
   blklen = mxCalloc(nr,sizeof(mwSize));
   blkend = mxCalloc(nr,sizeof(mwSize));   
   len = 0; numblk = 0;  NZ = 0; 
   for (k=0; k<nr; k++) {
      if (rr[k]==1) { len++; 
      } else { 
        if (len > 0) {
           blklen[numblk] = len;
           blkend[numblk] = k+1;
           NZ += len;
           numblk++;           
           len = 0; }         
      }
   }
   if (len > 0) {
      blklen[numblk] = len;
      blkend[numblk] = nr+1; 
      NZ += len; 
      numblk++;            
   }
   NZ += numblk+1; 
 /************************************************
  ************************************************/
   plhs[0] = mxCreateDoubleMatrix(nr+1,1,mxREAL);  
   plhs[1] = mxCreateSparse(nr+1,numblk,NZ,mxREAL);           
   hh = mxGetPr(plhs[0]);   
   ii = mxGetIr(plhs[1]);  
   jj = mxGetJc(plhs[1]);
   vv = mxGetPr(plhs[1]);  
   cnt = 0;    
   for (k=0; k<numblk; k++) {
      len = blklen[k]+1; 
      idxend = blkend[k];
      idxstart = idxend-len; 
      tmp = 1/sqrt(len);
      for (j=idxstart; j<idxend; j++) { 
         hh[j] = 1; 
         ii[cnt] = j;      
         vv[cnt] = tmp;
         cnt++; 
      }
      jj[k+1] = cnt;
   }
   for (k=0; k<=nr; k++) { 
      hh[k] = 1-hh[k]; 
   }
return;
}
/************************************************************************/

