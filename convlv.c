/*************************************************************************** 
CONVLV: This is the convlv algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Convolves or deconvolves a real data set data[0..n-1] (including any
user-supplied zero-padding) with a response function respns[0..n-1]. The
response function must be stored in wrap-around order in the first m
elements of respns, where m is an odd integer <= n. Wrap-around order means
that the first half of the array contains the impulse response function at
negative times, counting down from the highest element respns[m-1]. On
input isign is +1 for convolution, -1 for deconvolution. The answer is
returned in the first n components of ans. However, ans must be supplied in
the calling program with dimensions [0..2*n-1], for consistency with
twofft. n MUST be an integer power of two."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.
2) Error messages are only warnings and control is given back to the
   calling program which can then act on the "0" that is returned.
3) The user does not need to zero pad the data array and it can be of any
   length.
4) If the data array is very long then the routine will automatically, and
   without warning, break the data down into chunks and use the
   "overlap-add" method as outlined in NR.
5) The response function is entered as an array without wrap-around
   order. The additional input integer, cntr, is the array index of the
   element of respns to be used as the center when contructing a
   wrap-around order version of respns.
6) The ans array need not be of dimension [0..2n-1]. It should, however, be
   of dimension [0..n-1].
7) m can be even if the user so wishes.

The user still has to be careful that the data and respns arrays do not
contain array elements which are of very different magnitudes so as to
avoid round-off errors.

****************************************************************************/

#include <stdlib.h>
#include "fft.h"
#include "memory.h"
#include "error.h"

int convlv(double *data, unsigned long n, double *respns, unsigned long m,
	   unsigned long cntr, int isign, double *ans) {

  double          *dat,*rpns,*fft1,*fft2;
  double          dum,mag2;
  unsigned long   nfft=NFFT_MIN,nffto2=0,nffto4=0,nchunks=1;
  unsigned long   i=0,j=0,k=0,ndat=0,sdat=0,edat=0;

  /* Make sure m is not too long */
  if (m>(NFFT_MAX>>1)) {
    nferrormsg("convlv(): Response array length m(=%d) is too long!",m);
    return 0;
  }

  /* Initialize answer array */
  for (i=0; i<n; i++) ans[i]=0.0;

  /* Find length of fft array required */
  while ((nfft<<=1)<n+m); nfft<<=1;

  /* If data is too long, figure out how many smaller chunks we need */
  while (nfft>NFFT_MAX) {
    nfft>>=1;
    if (++nchunks>NCHKS)
      nferrormsg("convlv(): Data is very long and more than NCHKS=%d\n\
\tchunks are required! Try increasing NCHKS in fft.h.",NCHKS);
  }

  /* Make sure that once response array buffer is factored in, we still
     have enough space to work with */
  while ((nfft>>1)<(n/nchunks+nchunks-1+m)) {
    if (++nchunks>NCHKS)
      nferrormsg("convlv(): Data is very long and more than NCHKS=%d\n\
\tchunks are required! Try increasing NCHKS in fft.h.",NCHKS);
  }

  nffto2=nfft>>1;
  nffto4=nffto2>>1;
  /* Allocate memory for rpns array */
  if (!(rpns=darray(nffto2))) {
    nferrormsg("convlv(): Required length(=%d) of response array too high!",
	       nffto2); return 0;
  }

  /* Copy response array from user input to utility array in wrap-around
     order */
  for (i=cntr,j=0; i<m; i++,j++) *(rpns+j)=respns[i];
  for (i=0,j=nffto2-cntr; i<cntr; i++,j++) *(rpns+j)=respns[i];
  /* Pad with zeros */
  for (i=m-cntr; i<nffto2-cntr; i++) *(rpns+i)=0.0;

  /** Loop around chunks of data **/
  edat=0;
  for (i=0; i<nchunks; i++) {

    /* Allocate memory for dat, rpns and fft arrays */
    if (!(dat=darray(nffto2))) {
      nferrormsg("convlv(): Required length(=%d) of dat array too high!",
		 nffto2); return 0;
    }
    if (!(fft1=darray(nfft))) {
      nferrormsg("convlv(): Required length(=%d) of FFT arrays too high.\n\
\tTry reducing NFFT_MAX(=%d) in fft.h.",nfft,NFFT_MAX); return 0;
    }
    if (!(fft2=darray(nfft))) {
      nferrormsg("convlv(): Required length(=%d) of FFT arrays too high.\n\
\tTry reducing NFFT_MAX(=%d) in fft.h.",nfft,NFFT_MAX); return 0;
    }

    /* Find start and end of present chunk */
    sdat=edat;
    if (i==nchunks-1) edat=n;
    else edat=sdat+n/nchunks;
    ndat=edat-sdat;

    /* Break up data array from user input to individual chunk */
    for (j=sdat,k=cntr; j<edat; j++,k++) *(dat+k)=data[j];
    /* Pad on both sides with zeros */
    for (j=0; j<cntr; j++) *(dat+j)=0.0;
    for (j=cntr+ndat; j<nffto2; j++) *(dat+j)=0.0;
    
    /* FFT both the data and response arrays */
    twodblefft(dat,rpns,fft1,fft2,nffto2);

    /* Multiply FFTs together to convole */
    if (isign) {
      for (j=1; j<nffto2+3; j+=2) {
	fft2[j-1]=(fft1[j-1]*(dum=fft2[j-1])-fft1[j]*fft2[j])/nffto4;
	fft2[j]=(fft1[j]*dum+fft1[j-1]*fft2[j])/nffto4;
      }
    }
    /* Or deconvolve */
    else {
      for (j=1; j<nffto2+3; j+=2) {
	if ((mag2=fft2[j-1]*fft2[j-1]+fft2[j]*fft2[j])==0.0) {
	  nferrormsg("convlv(): Deconvolving at response zero"); return 0;
	}
	fft2[j-1]=(fft1[j-1]*(dum=fft2[j-1])+fft1[j]*fft2[j])/nffto4/mag2;
	fft2[j]=(fft1[j]*dum-fft1[j-1]*fft2[j])/nffto4/mag2;
      }
    }

    /* Pack last element with first for dblefft */
    fft2[1]=fft2[nffto2];
    
    /* FFT back to real space */
    dblefft(fft2,nffto2,-1);

    /* Add convolution to the ans array */
    if (!(j=sdat)) {
      while (j+cntr<nffto2 && j<n) {
	ans[j]+=fft2[j+cntr]; j++;
      }
    }
    else {
      j-=cntr; k=0; while (k<nffto2 && j<n) {
	ans[j]+=fft2[k]; j++; k++;
      }
    }

    /* Free the allocated space */
    free(dat); free(fft1); free(fft2);
    
  }

  /* Free the allocated space */
  free(rpns);
  
  return 1;

}
