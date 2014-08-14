#include <stdio.h>
#include <fitsio.h>

void fitserrmsg(int errstatus) {
  char buffer[64]="\0";
  if(!errstatus) 
    return;
  else {
    fits_get_errstatus(errstatus,buffer);
    fprintf(stderr,"FITS ERROR #%i\n\t%s\n\n",errstatus,buffer);
    exit(1);
  }
}
