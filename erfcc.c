/*************************************************************************** 
ERFCC: This is the erfcc algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Returns the complementary error function erfc(x) with fractional error
everywhere less than 1.2x10^{-7}"

NOTE: A more exacting version of the same function is given as in
erffc.c. This is the proper definition of the complimentary error function.

I have changed the raw routine to include the following features:
0) Double precision is used.

****************************************************************************/

#include <math.h>

double erfcc(double x) {

  double  t=0.0,z=0.0,ans=0.0;

  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  ans=t*exp(-z*z-1.26551223+t*
	    (1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*
            (0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+
	    t*0.17087277)))))))));

  return x >= 0.0 ? ans : 2.0-ans;

}
