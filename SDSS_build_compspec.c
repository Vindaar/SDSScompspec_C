#include "SDSS-spec.h"
#include "const.h"
#include "error.h"


int FindBin(double wlarray[], double wlx, int nzpart) {
  /* This program sorts through the bins looking for which bin the bottom
     pixel lies in    */
  int i=0;
  do{
    i++;
  } while(wlarray[i]<=wlx && i<201);
  i--;
  
  // TODO: Check what to return? I typed in i just so it has a return.
  return i;

}

int SDSS_build_compspec(spectrum *spec, double wlcut[], int nzpart, int forest, compspec *cspec) {

  int i=0;
  return 0;

}
  
