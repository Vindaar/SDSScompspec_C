/*************************************************************************** 
PG_ADDOBJ: Initialize and/or add to an array of plot objects. 
****************************************************************************/

#include <stdlib.h>
#include "pg_plot.h"
#include "error.h"

int pg_addobj(double x, double y, int m, int s, int c, plobj **obj, int *nobj) {

  /* Check number of objects is reasonable */
  if (*nobj<0) {
    nferrormsg("pg_addobj(): Number of objects (=%d) must be >=0",*nobj); return 0;
  }

  /* Reallocate memory for object array */
  if (!((*obj)=(plobj *)realloc((*obj),(size_t)((*nobj+1)*sizeof(plobj))))) {
    nferrormsg("pg_addobj(): Cannot reallocate memory for obj array\n\
\tof size %d",*nobj+1); return 0;
  }

  /* Fill in last array element */
  (*obj)[*nobj].x=x; (*obj)[*nobj].y=y; (*obj)[*nobj].m=m; (*obj)[*nobj].s=s;
  (*obj)[*nobj].c=c;

  /* Increment number of objects */
  (*nobj)++;

  return 1;

}
