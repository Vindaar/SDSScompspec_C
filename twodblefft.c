/*************************************************************************** 
TWODBLEFFT: This is the twofft algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Given two real input arrays data1[1..n] and data2[1..n], this routine
calls four1 and returns two complex output arrays, fft1[1..2n] and
fft2[1..2n], each of complex length n (i.e. real length 2n), which contain
the discrete Fourier transforms of the respective data arrays. n MUST be an
integer power of 2."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.

****************************************************************************/

#include "fft.h"
#include "error.h"

int twodblefft(double *data1, double *data2, double *fft1, double *fft2,
	       unsigned long n) {

  double          rep,rem,aip,aim;
  unsigned long   nn1,nn,j=0,jj=0;

  nn1=1+(nn=n+n);
  for (j=0,jj=1; j<n; j++,jj+=2) {
    fft1[jj-1]=data1[j];
    fft1[jj]=data2[j];
  }

  /* Fourier transform the data */
  fourier(fft1,n,1);
  /* fourier(fft1-1,n,1); */
  
  fft2[0]=fft1[1];
  fft1[1]=fft2[1]=0.0;
  for (j=2; j<n+1; j+=2) {
    rep=0.5*(fft1[j]+fft1[nn-j]);
    rem=0.5*(fft1[j]-fft1[nn-j]);
    aip=0.5*(fft1[j+1]+fft1[nn1-j]);
    aim=0.5*(fft1[j+1]-fft1[nn1-j]);
    fft1[j]=rep;
    fft1[j+1]=aim;
    fft1[nn-j]=rep;
    fft1[nn1-j]=-aim;
    fft2[j]=aip;
    fft2[j+1]=-rem;
    fft2[nn-j]=aip;
    fft2[nn1-j]=rem;
  }

  return 1;

}
