/***************************************************************************
FFT.H: Include file for operations concerning the Fast Fourier Transform
***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
#define NFFT_MIN    32
#define NFFT_MAX 16384
#define NCHKS       40

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* STRUCTURES */

/* PROTOTYPES */

int convlv(double *data, unsigned long n, double *respns, unsigned long m,
	   unsigned long cntr, int isign, double *ans);
int dblefft(double *data, unsigned long n, int isign);
int fourier(double *data, unsigned long n, int isign);
int twodblefft(double *data1, double *data2, double *fft1, double *fft2,
               unsigned long n);
