/****************************************************************************
* Qsort comparison routine for a single float array
****************************************************************************/

int qsort_farray(const void *x1, const void *x2) {

  if ((*(float *)x1) > (*(float *)x2)) return 1;
  else if ((*(float *)x1) == (*(float *)x2)) return 0;
  else return -1;

}
