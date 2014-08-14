/*************************************************************************** 
PG_GET_NEAROBJ: Find the nearest object to the input co-ords
(typically returned from CPGBAND). The routine returns -1 if an error
occurs, otherwise it returns the index of the object in the
user-supplied object-list which is closest to the coords xpos,
ypos. Note that in order to compare distances in x and y fairly, the
user must supply a normalisation constant for both the x and y
directions. Typically this will be the size of the plotting surface in
the x and y directions in the same units as used for xpos and ypos
(usually world coordinates). With the opt argument, the user can
control whether a match is searched for among point, horizontal lines
and/or veritcal lines in different combinations:

Opt  :   Objects searched
0    :   Points only
1    :   Horizontal lines only
2    :   Vertical lines only
3    :   Points, horizontal and vertical lines
4    :   Points and horizontal lines
5    :   Points and vertical lines
6    :   Horizontal and vertical lines

Note that dummy values for xnorm and ynorm (e.g. 1.0, NOT 0.0) can be used
when options 1 and 2.

****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "pg_plot.h"
#include "error.h"

int pg_get_nearobj(float xpos, float ypos, double xnorm, double ynorm,
		   plobj *obj, int nobj, int opt) {

  double x=0.0,y=0.0,d=0.0,min=0.0;
  int    i=0,j=0,n=0;

  /* Basic initial checks */
  if (nobj<1) {
    nferrormsg("pg_get_nearobj(): nobj = %d, should be >= 1",nobj); return -1;
  }
  if (obj==NULL) {
    nferrormsg("pg_get_nearobj(): Pointer to object array is NULL"); return -1;
  }
  if (xnorm<=0.0) {
    nferrormsg("pg_get_nearobj(): X-normalisation factor (=%lf) must be >0.0",
	       xnorm); return -1;
  }
  if (ynorm<=0.0) {
    nferrormsg("pg_get_nearobj(): Y-normalisation factor (=%lf) must be >0.0",
	       ynorm); return -1;
  }

  if (opt==0) {
    for (i=0,j=0,n=0; i<nobj; i++) {
      if (obj[i].m==POINT) {
	n++; x=(obj[i].x-(double)xpos)/xnorm; y=(obj[i].y-(double)ypos)/ynorm; 
	if ((d=x*x+y*y)<min || n==1) { min=d; j=i; }
      }
    }
  }
  else if (opt==1) {
    for (i=0,j=0,n=0; i<nobj; i++) {
      if (obj[i].m==HLINE) {
	n++; if ((d=abs(obj[i].y-(double)ypos))<min || n==1) { min=d; j=i; }
      }
    }
  }
  else if (opt==2) {
    for (i=0,j=0,n=0; i<nobj; i++) {
      if (obj[i].m==VLINE) {
	n++; if ((d=abs(obj[i].x-(double)xpos))<min || n==1) { min=d; j=i; }
      }
    }
  }
  else if (opt==3) {
    for (i=0,j=0,n=0; i<nobj; i++) {
      if (obj[i].m==POINT || obj[i].m==HLINE || obj[i].m==VLINE) {
	n++;
	switch (obj[i].m) {
	case POINT:
	  x=(obj[i].x-(double)xpos)/xnorm; y=(obj[i].y-(double)ypos)/ynorm;
	  d=x*x+y*y; break;
	case HLINE:
	  y=abs(obj[i].y-(double)ypos)/ynorm; d=y*y; break;
	case VLINE:
	  x=abs(obj[i].x-(double)xpos)/xnorm; d=x*x; break;
	}
	if (d<min || n==1) { min=d; j=i; }
      }
    }
  }
  else if (opt==4) {
    for (i=0,j=0,n=0; i<nobj; i++) {
      if (obj[i].m==POINT || obj[i].m==HLINE) {
	n++;
	switch (obj[i].m) {
	case POINT:
	  x=(obj[i].x-(double)xpos)/xnorm; y=(obj[i].y-(double)ypos)/ynorm;
	  d=x*x+y*y; break;
	case HLINE:
	  y=abs(obj[i].y-(double)ypos)/ynorm; d=y*y; break;
	}
	if (d<min || n==1) { min=d; j=i; }
      }
    }
  }
  else if (opt==5) {
    for (i=0,j=0,n=0; i<nobj; i++) {
      if (obj[i].m==POINT || obj[i].m==VLINE) {
	n++;
	switch (obj[i].m) {
	case POINT:
	  x=(obj[i].x-(double)xpos)/xnorm; y=(obj[i].y-(double)ypos)/ynorm;
	  d=x*x+y*y; break;
	case HLINE:
	  x=abs(obj[i].x-(double)xpos)/xnorm; d=x*x; break;
	}
	if (d<min || n==1) { min=d; j=i; }
      }
    }
  }
  else if (opt==6) {
    for (i=0,j=0,n=0; i<nobj; i++) {
      if (obj[i].m==HLINE || obj[i].m==VLINE) {
	n++;
	switch (obj[i].m) {
	case HLINE:
	  y=abs(obj[i].y-(double)ypos)/ynorm; d=y*y; break;
	case VLINE:
	  x=abs(obj[i].x-(double)xpos)/xnorm; d=x*x; break;
	}
	if (d<min || n==1) { min=d; j=i; }
      }
    }
  }
  else { nferrormsg("pg_get_nearobj(): Invalid option flag, %d",opt); return -1; }

  return j;

}
