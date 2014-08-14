/*************************************************************************** 
KSPROB: This is the Kolmogorov-Smirnov probability routine, probks, given
in Numerical Recipes. Here's what NR has to say about it:

Kolmogorov-Smirnov probability function.

I have made the routine double precision.
****************************************************************************/

#include <math.h>
#include "stats.h"
#include "error.h"

double ksprob(double alam) {

  int     i=0,j=0;
  double  a2,fac=2.0,sum=0.0,term=0.0,termbf=0.0;

  a2 = -2.0*alam*alam;
  for (i=0, j=1; i<=KS_ITER; i++, j++) {
    term=fac*exp(a2*((double)j*j));
    sum += term;
    if (fabs(term) < KS_EPS1*termbf || fabs(term) < KS_EPS2*sum) return sum;
    fac = -fac;
    termbf=fabs(term);
  }

  warnmsg("ksprob(): Failed to converge in %d iterations",KS_ITER);
  return 1.0;
}
