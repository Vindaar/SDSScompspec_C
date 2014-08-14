/****************************************************************************
FHIST.C: Routine to form the float precision histogram of some float
precision data.

Note that the histogram is not necessarily cleared first, meaning that you
can add to an old histogram from different arrays by calling this routine
many times.

- data contains the data to be binned. This is an array of length ndat.
- hist contains the input and output histrogram of length nhist.
- x is an array of length nhist+1 which contains the left edges of each
  bin. The final point in the x array specifies the right edge of the
  final bin.
- If opt1 = 0 then the histogram is not cleared before filling.
  If opt1 = 1 then the histogram is cleared before filling
- If opt2 = 0 then the x array is internally determined: the left edges of
  the bins are taken such that the bins are of equal size and the full
  range of data is used.
  If opt2 = 1 then the current values in x are used and data values outside
  of [x[0],x[nhist]] are excluded. If this option is selected and the x
  array has not been externally set, a warning is issued and internal
  determination continues as if opt2=0 was entered
****************************************************************************/

#include <stddef.h>
#include "error.h"

int fhist(float *data, int ndat, float *hist, float *x, int nhist, int opt1,
	  int opt2) {

  float  binsize=0.0;
  int    i=0,j=0,determine_x=0;

  /* Clear histogram if required by user */
  if (opt1) for (i=0; i<nhist; i++) hist[i]=0.0;

  /* Test x array if required */
  if (opt2 && (x==NULL || x[0]==x[nhist])) {
    warnmsg("fhist(): x array is either NULL or does not contain\n\
\tsensible binning information. Determining binning scale interanlly.");
    determine_x=1;
  }

  /* Determine x array if required */
  if (!opt2 || determine_x) {
    x[0]=x[nhist]=data[0];
    for (i=1; i<ndat; i++) {
      if (data[i]<x[0]) x[0]=data[i]; if (data[i]>x[nhist]) x[nhist]=data[i];
    }
    binsize=(x[nhist]-x[0])/(float)nhist;
    for (i=1; i<nhist; i++) x[i]=x[0]+binsize*(float)i;
  }

  /* Creat histogram */
  for (i=0; i<ndat; i++) {
    if (data[i]>=x[0] && data[i]<=x[nhist]) {
      j=1;
      while (data[i]>=x[j]) {
	j++; if (j==nhist) break;
      }
      hist[j-1]++;
    }
  }

  return 1;

}
