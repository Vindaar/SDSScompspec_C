/*************************************************************************** 
DBLEFFT: This is the realfft algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Calculates the Fourier transform of a set of n real-valued data
points. Replaces this data (which is stored in array data[1..n]) by the
positive frequency half of its complex Fourier transform. The real-valued
first and last components of the complex transform are returned as elements
data[1] and data[2], respectively. n MUST be a power of 2. This routine
also calculates the inverse transform of a complex data array if it is the
tranform of real data (result in this case must be multiplied by 2/n)."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.

****************************************************************************/

#include "fft.h"
#include "const.h"
#include "error.h"

int dblefft(double *data, unsigned long n, int isign) {

  double             c1=0.5,c2,h1r,h1i,h2r,h2i;
  double             wr,wi,wpr,wpi,wtemp,theta;
  unsigned long      i,i1,i2,i3,i4,np1;
  
  theta=C_PI/(double)(n>>1);
  if (isign==1) {
    c2=-0.5;
    fourier(data,n>>1,1);
  }
  else {
    c2=0.5;
    theta=-theta;
  }

  wtemp=sin(0.5*theta);
  wpr=-2.0*wtemp*wtemp;
  wpi=sin(theta);
  wr=1.0+wpr;
  wi=wpi;
  np1=n+1;

  for (i=1; i<(n>>2)+2; i++) {
    i4=1+(i3=np1-(i2=1+(i1=i+i)));
    h1r=c1*(data[i1]+data[i3]);
    h1i=c1*(data[i2]-data[i4]);
    h2r=-c2*(data[i2]+data[i4]);
    h2i=c2*(data[i1]-data[i3]);

    data[i1]=h1r+wr*h2r-wi*h2i;
    data[i2]=h1i+wr*h2i+wi*h2r;
    data[i3]=h1r-wr*h2r+wi*h2i;
    data[i4]=-h1i+wr*h2i+wi*h2r;

    wr=(wtemp=wr)*wpr-wi*wpi+wr;
    wi=wtemp*wpi+wi*wpr+wi;
  }
  if (isign==1) {
    data[0]=(h1r=data[0])+data[1];
    data[1]=h1r-data[1];
  } else {
    data[0]=c1*((h1r=data[0])+data[1]);
    data[1]=c1*(h1r-data[1]);
    fourier(data,n>>1,-1);
  }

  return 1;

}
