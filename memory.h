/***************************************************************************
MEMORY.H: Include file for allocation of memory routines.
***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
#define MEM_H 1 /* This include file has been called */

/* STRUCTURES */
typedef struct DblePair{
  double x1;
  double x2;
} dblepair;

/* PROTOTYPES */
char *carray(unsigned long n);
char **cmatrix(unsigned long n, unsigned long m);
double *darray(unsigned long n);
double **dmatrix(unsigned long n, unsigned long m);
double ***dcube(unsigned long n, unsigned long m, unsigned long p);
float *farray(unsigned long n);
float **fmatrix(unsigned long n, unsigned long m);
float ***fcube(unsigned long n, unsigned long m, unsigned long p);
int *iarray(unsigned long n);
int **imatrix(unsigned long n, unsigned long m);
int ***icube(unsigned long n, unsigned long m, unsigned long p);
