/***************************************************************************** 
POLCOF: This is the POLCOF algorithm taken from Numerical
Recipes. Here's what NR has to say about it:

"Given arrays xa[0..n] and ya[0..n] containing a tabulated function
yai = f(xai), this routine returns an array of coefficients cof[0..n]
such that yai = j cofjxaj i."

I have made the following changes:
0) Double precision inputs and outputs and returns an interger
1) Own error handling routines.
2) Have made consistent with zero-offsetting of this routine with dpolint()
*****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "fit.h"
#include "memory.h"
#include "error.h"

int dpolcof(double *xa, double *ya, int n, double *cof) {

  double xmin=0.0,dy=0.0;
  double *x=NULL,*y=NULL;
  int    i=0,j=0,k=0;

  /* Allocate memory for temporary arrays and fill them */
  if ((x=darray(n))==NULL) {
    errormsg("dpolcof(): Cannot allocate memory to x array\n\to size %s",n);
    return 0;
  }
  if ((y=darray(n))==NULL) {
    errormsg("dpolcof(): Cannot allocate memory to y array\n\to size %s",n);
    return 0;
  }
  for (i=0; i<n; i++) { x[i]=xa[i]; y[i]=ya[i]; }
  for (j=0; j<n; j++) {
    if (!dpolint(x,y,n-j,0.0,&cof[j],&dy)) {
      nferrormsg("dpolcof(): Error returned from dpolint()"); return 0;
    }
    xmin=FIT_INFIN; k=-1;
    for (i=0; i<n-j; i++) {
      if (abs(x[i])<xmin) { xmin=fabs(x[i]); k=i; }
      if (x[i]) y[i]=(y[i]-cof[j])/x[i];
    }
    for (i=k+1; i<n-j; i++) { y[i-1]=y[i]; x[i-1]=x[i]; }
  }

  /* Clean up */
  free(x); free(y);

  return 1;

}
