#include <stdlib.h>
#include <math.h>
#include "error.h"

double res_func(double wl, int spectro){

  double res=0;
  if(spectro == 1)
    res=1.27464 - .07071 * wl/1000;
  else if(spectro == 2)
    res =( wl > 6150 ?( 5.342 - 1.225*wl/1000 + .08333 * pow(wl/1000,2) ) :
                       ( 1.32 - .1503 * wl/1000 +.0114 * pow(wl/1000,2) ));

  else
    errormsg("Invalid Spectrograph Entry");

  return res;
}
