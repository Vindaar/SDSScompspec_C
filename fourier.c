/*************************************************************************** 
FOURIER: This is the four1 algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Replaces data[1..2n] by its discrete Fourier tranform, if isign is input
as 1; or replaces data[1..2n] by n times its inverse discrete Fourier
transform, if isign is input as -1. data is a complex array of length n or,
equivalently, a real array of length 2n. n MUST be an integer power of 2
(this is not checked for)."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.
2) Arrays start with index 0

****************************************************************************/

#include "fft.h"
#include "const.h"
#include "error.h"

int fourier(double *data, unsigned long nn, int isign) {

  double             tempr,tempi;
  double             wtemp,wr,wpr,wpi,wi,theta;
  unsigned long      n,mmax,m,j,istep,i;

  n=nn<<1;
  j=0;
  for (i=0; i<n-1; i+=2) {
    if (j>i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n>>1;
    while (m>=2 && j>m-1) {
      j-=m;
      m>>=1;
    }
    j+=m;
  }

  mmax=2;
  while (n>mmax) {
    istep=mmax<<1;
    theta=isign*(2.0*C_PI/mmax);
    wtemp=sin(0.5*theta);
    wpr=-2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=0; m<mmax-1; m+=2) {
      for (i=m; i<n; i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i]+=tempr;
	data[i+1]+=tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wtemp*wpi+wi*wpr+wi;
    }
    mmax = istep;
  }

  return 1;

}
