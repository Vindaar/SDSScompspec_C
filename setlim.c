/* Allows the User to set the limit which will be used as
Wlimit in SDSS_bounds */


#include <stdlib.h>
#include "SDSS-spec.h"
#include "error.h"

int setlim(void){

float newlim;
int set;

  while(!set){
    printf("Set the PE threshhold:   ")
    scanf("%f", &newlim);
    if(newlim < .001 || newlim > .5)
      printf("PE threshhold should be in range [.001,.5]\n");
    else
      Wlimit = newlim;
      set++;
  }

return 1;
}
