/**********************************************************************************

  DELTAGAMMAINC Fast and Accurate Evaluation of a Generalized Incomplete Gamma
  Function. Copyright (C) 2016 Remy Abergel (remy.abergel AT parisdescartes.fr), 
  Lionel Moisan (Lionel.Moisan AT parisdescartes.fr).

  This file is a part of the DELTAGAMMAINC software, dedicated to the
  computation of a generalized incomplete gammafunction. See the Companion paper
  for a complete description of the algorithm.

  ``Fast and accurate evaluation of a generalized incomplete gamma function''
  (Rémy Abergel, Lionel Moisan), preprint MAP5 nº2016-14, revision 1.

  This program is free software: you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/

/*
 * Description : MEX interface for deltagammainc.
 * 
 * compilation (from the Matlab console) : mex -R2018a -silent -lm CFLAGS="\$CFLAGS -std=c99" kernel.c deltagammainc_mexinterface.c -output deltagammainc_mexinterface
 * usage : [rho,sigma] = deltagammainc_mexinterface(x,y,mu,p);
 *
 */

#include <mex.h>
#include <math.h>
#include "matrix.h"

extern void deltagammainc(double*,double*,char*,double,double,double,double); // see source file 'kernel.c'

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  double *x,*y,*mu,*p,*rho,*sigma;
  size_t size_x,size_y,size_mu,size_p;
  char method;
  int adr;

  /*******************/
  /* retrieve inputs */
  /*******************/
  
  // number of inputs //
  if(nrhs != 4) {
    //display_usage();
    mexErrMsgIdAndTxt("deltagammainc:deltagammainc_mexinterface:badinputnumber","Incorrect number of inputs");
  }

  // number of outputs //
  if(nlhs > 2) {
    //display_usage();
    mexErrMsgIdAndTxt("deltagammainc:deltagammainc_mexinterface:badoutputnumber","Incorrect number of outputs");
  }

  // first input (x) //
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
    //display_usage();
    mexErrMsgIdAndTxt("deltagammainc:deltagammainc_mexinterface:badinputdatatype","the first input (x) must be of type 'double'.");
  }
  size_x = mxGetNumberOfElements(prhs[0]);
  x = mxGetDoubles(prhs[0]);
  
  // second input (y) //
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
    //display_usage();
    mexErrMsgIdAndTxt("deltagammainc:deltagammainc_mexinterface:badinputdatatype","the second input (y) must be of type 'double'.");
  }
  size_y = mxGetNumberOfElements(prhs[1]);
  y = mxGetDoubles(prhs[1]);

  // third input (mu) //
  if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
    //display_usage();
    mexErrMsgIdAndTxt("deltagammainc:deltagammainc_mexinterface:badinputdatatype","the second input (2) must be of type 'double'.");
  }
  size_mu = mxGetNumberOfElements(prhs[2]);
  mu = mxGetDoubles(prhs[2]);

  // fourth input (p) //
  if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
    //display_usage();
    mexErrMsgIdAndTxt("deltagammainc:deltagammainc_mexinterface:badinputdatatype","the second input (2) must be of type 'double'.");
  }
  size_p = mxGetNumberOfElements(prhs[3]);
  p = mxGetDoubles(prhs[3]);

  /*************************************/
  /* test input retrivial (debug only) */
  /*************************************/
  /*
  mexPrintf("numel(x) = %d\n",size_x);
  for(adr=0;adr<size_x;adr++)
   mexPrintf("x(%d) = %g\n",adr+1,x[adr]);
  mexPrintf("numel(y) = %d\n",size_y);
  for(adr=0;adr<size_y;adr++)
   mexPrintf("y(%d) = %g\n",adr+1,y[adr]);
  mexPrintf("numel(mu) = %d\n",size_mu);
  for(adr=0;adr<size_mu;adr++)
   mexPrintf("mu(%d) = %g\n",adr+1,mu[adr]);
  mexPrintf("numel(p) = %d\n",size_p);
  for(adr=0;adr<size_p;adr++)
   mexPrintf("p(%d) = %g\n",adr+1,p[adr]);
  */
    
  /*********************/
  /* memory allocation */
  /*********************/
  plhs[0] = mxCreateDoubleMatrix((mwSize)size_x,(mwSize)1,mxREAL); // rho  
  plhs[1] = mxCreateDoubleMatrix((mwSize)size_x,(mwSize)1,mxREAL); // sigma
  rho = mxGetData(plhs[0]);
  sigma = mxGetData(plhs[1]);
    
  /***************** kernel *****************/
  for(adr=0;adr<size_x;adr++) {
    deltagammainc(rho+adr,sigma+adr,&method,x[adr],y[adr],mu[adr],p[adr]);
  }
  /**************** end kernel ***************/
  
}
