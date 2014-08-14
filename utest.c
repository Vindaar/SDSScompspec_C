/*************************************************************************** 
UTEST: Calculates the probability that two distributions are drawn from the
same parent distribution using the Mann-Whitney U-test. The two input data
arrays are combined and indexed in a way that allows easy
de-entanglement. The combined data and index arrays are simultaneously
sorted and the data array is replaced by the rank array. The index array is
used to calculate the sum of ranks for each data sample and the U is
calculated for each. The smaller value of U is taken using U_2 = n_1 * n_2
- U_1. The distribution of z(U) is approximately normal for sample sizes
>~20 and so the erfcc routine is used to calculate the probabiliy P(U).

NOTE: The U-test only works well for samples sizes >~20.
NOTE: The input data arrays are corrupted during this routine.

****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "gamm.h"
#include "sort.h"
#include "stats.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int utest(double *data1, unsigned long n1, double *data2, unsigned long n2,
	  double *prob) {

  double        U=0.0,z=0.0;
  double        ddum=0.0;
  double        *wksp=NULL,*idx=NULL;
  int           i=0;
  unsigned long N=0,n1n2=0;

  N=n1+n2; n1n2=n1*n2;

  /* Allocate memory for workspace arrays */
  if ((wksp=darray(N))==NULL)
    errormsg("utest(): Cannot allocate memory for wksp array\n\
\tof size %d",N);
  if ((idx=darray(N))==NULL)
    errormsg("utest(): Cannot allocate memory for wksp array\n\
\tof size %d",N);

  /* Fill combined data and index array */
  for (i=0; i<n1; i++) { wksp[i]=data1[i]; idx[i]=1.0; }
  for (i=0; i<n2; i++) { wksp[n1+i]=data2[i]; idx[n1+i]=-1.0; }

  /* Simultaneously sort data arrays */
  if (!sort_2darray(N,wksp,idx)) {
    nferrormsg("utest(): Problem sorting workspace arrays"); return 0;
  };

  /* Rank data array */
  crank(N,wksp,&ddum);

  /* Sum ranks for smallest sample and calculate U */
  if (n1<n2) {
    for (i=0; i<N; i++) if (idx[i]>0.0) U-=wksp[i];
    U+=(double)(n1n2+n1*(n1+1)/2);
  }
  else {
    for (i=0; i<N; i++) if (idx[i]<0.0) U-=wksp[i];
    U+=(double)(n1n2+n2*(n2+1)/2);
  }

  /* Take smallest U using normalization relation */
  if ((ddum=(double)n1n2-U)<U) U=ddum;

  /* Calculate z statistic */
  z=sqrt(12.0)*(U-(double)n1n2/2.0)/sqrt((double)(n1n2)*(double)(N+1));

  /* Calculate the probability of z */
  *prob=erfcc(C_1_SQRT2*fabs(z));

  /* Clean up */
  free(wksp); free(idx);

  return 1;

}
