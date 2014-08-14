/****************************************************************************
* Qsort comparison routine for a array of double precision pairs.
****************************************************************************/

#include "sort.h"

int qsort_twodarray(const void *dat1, const void *dat2) {

  if (((twodble *)dat1)->a > ((twodble *)dat2)->a) return 1;
  else if (((twodble *)dat1)->a == ((twodble *)dat2)->a) return 0;
  else return -1;

}
