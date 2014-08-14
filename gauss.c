/****************************************************************************
* Unnormalized Gaussian function
****************************************************************************/

#include <math.h>

double gauss(double x) {

  double xsq=0.0;

  if ((xsq=x*x)>84.0) return 0.0;
  else return exp(-xsq);
}
