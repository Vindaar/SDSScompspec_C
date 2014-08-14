/****************************************************************************
* Find the index of the maximum number in an array of n floats. Returns -1
* if an error occurs, usually because an n<0 has been passed.
****************************************************************************/

int idxfmax(float *array, int n) {

  int          i=0,imax=0;
  float        max;

  if (n<0) return -1;
  max=array[0];
  for (i=1; i<n; i++) {
    if (array[i]>max) {
      max=array[i];
      imax=i;
    }
  }

  return imax;

}
