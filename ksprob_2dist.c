/*************************************************************************** 
KSPROB_2DIST: Calculates the KS probabilitythat two given distributions are
drawn from the same distribution. This is the kstwo routine given in
Numerical Recipes. Here's what NR has to say about it:

Given an array data1[1..n1], and an array data2[1..n2], this routine
returns the KS statistic d, and the significance level prob for the null
hypothesis that the data sets are drawn from the same distribution. Small
values of prob show that the cumulative distribution function of data1 is
significantly different from that of data2. The arrays data1 and data2 are
modified by being sorted into ascending order.

I have made routine an integer function. Have made the routine double
precision. Made some arguments pointer rather than passed parameters or
lists. Have replaced sort with qsort.

NOTE: Data array are returned have been sorted
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "stats.h"
#include "sort.h"

int ksprob_2dist(double *data1, unsigned long n1, double *data2,
		 unsigned long n2, double *d, double *prob) {

  unsigned long j1=0,j2=0;
  double        d1=0.0,d2=0.0,dt=0.0,en1=0.0,en2=0.0,en=0.0,fn1=0.0,fn2=0.0;

  /* Sort data into ascending order */
  qsort(data1,n1,sizeof(double),qsort_darray);
  qsort(data2,n2,sizeof(double),qsort_darray);
  
  en1=(double)n1; en2=(double)n2; *d=0.0;
  
  while (j1<n1 && j2<n2) {
    if ((d1=data1[j1]) <= (d2=data2[j2])) fn1=((double)++j1)/en1;
    if (d2 <= d1) fn2=((double)++j2)/en2;
    if ((dt=fabs(fn2-fn1)) > *d) *d=dt;
  }
  en=sqrt(en1*en2/(en1+en2));
  *prob=ksprob((en+0.12+0.11/en)*(*d));

  return 1;

}
