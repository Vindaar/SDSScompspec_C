/***************************************************************************
* Send a fatal error message to the terminal and exit
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
//#include <fitsio.h>

void errormsg(char *fmt, ...) {

  va_list       args;
  extern char   *progname;

  va_start(args, fmt);
  fprintf(stderr, "\a%s: FATAL ERROR: ", progname);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
  exit(1);
}


/***************************************************************************
* Send a non-fatal error message to the terminal
***************************************************************************/


void nferrormsg(char *fmt, ...) {

  va_list       args;
  extern char   *progname;

  va_start(args, fmt);
  fprintf(stderr, "\a%s: ERROR: ", progname);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");

}
/***************************************************************************
* Send a warning message to the terminal
***************************************************************************/

void warnmsg(char *fmt, ...) {

  va_list       args;
  extern char   *progname;

  va_start(args, fmt);
  fprintf(stderr, "%s: WARNING: ", progname);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n");
}
/*
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
*/
