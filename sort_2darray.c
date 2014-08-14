/****************************************************************************
* Sort 2 double precision arrays based on the sort of the first array
****************************************************************************/

#include <stdlib.h>
#include "sort.h"
#include "error.h"

int sort_2darray(int n, double *data1, double *data2) {

  int     i=0;
  twodble *sarray;

  /* Allocate memory for the sort array */
  if (!(sarray=(twodble *)malloc((size_t)(n*sizeof(twodble))))) {
    nferrormsg("sort_2darray(): Cannot allocate memory for array");
    return 0;
  }

  /* Fill sort array */
  for (i=0; i<n; i++) {
    sarray[i].a=data1[i]; sarray[i].b=data2[i]; 
  }

  /* Sort the sort array */
  qsort(sarray,n,sizeof(twodble),qsort_twodarray);

  /* Fill sorted data arrays */
  for (i=0; i<n; i++) {
    data1[i]=sarray[i].a; data2[i]=sarray[i].b;
  }

  /* Clean up */
  free(sarray);

  return 1;

}
