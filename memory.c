/****************************************************************************
* Allocate a matrix structure of double precision with size row x col.
* Function returns a pointer to the matrix structure.
****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "error.h"

  
/****************************************************************************
* Allocate an array of double precision with length n. Function returns a
* pointer to first element
****************************************************************************/

double *darray(unsigned long n) {
  
  double   *da=NULL;
  int      i=0;

  if (!(da=(double *)malloc((size_t)(n*sizeof(double))))) {
    nferrormsg("darray(): Cannot allocate memory for array");
    return NULL;
  }

  /* Initialize array */
  for (i=0; i<n; i++) *(da+i)=0.0;

  return da;

}



/****************************************************************************
* Allocate a cube of doubles with dimensions n x m x q. Function
* returns a pointer to (0,0,0) element
****************************************************************************/

double ***dcube(unsigned long n, unsigned long m, unsigned long p) {
  
  double ***dc=NULL;
  int    i=0,j=0;

  /* Allocate array of pointers to pointers */
  if (!(dc=(double ***)malloc((size_t)(n*sizeof(double **))))) {
    nferrormsg("dcube(): Cannot allocate memory for pointer\n\
\tto pointer array");
    return NULL;
  }

  /* Allocate array of pointers */
  if (!(*dc=(double **)malloc((size_t)(n*m*sizeof(double *))))) {
    nferrormsg("dcube(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize block of memory */
  if ((**dc=darray(n*m*p))==NULL) {
    nferrormsg("dcube(): Cannot allocate memory for array\n\
\tof size %dx%dx%d",n,m,p);
    return NULL;
  }

  /* Set pointers to pointers and pointers to rows */
  for (i=0; i<n; i++) {
    if (i) { dc[i]=dc[i-1]+m; dc[i][0]=dc[i-1][m-1]+p; }
    for (j=1; j<m; j++) dc[i][j]=dc[i][j-1]+p;
  }

  return dc;

}


/****************************************************************************
* Allocate a matrix of double precision with dimensions n x m. Function
* returns a pointer to (0,0) element
****************************************************************************/

double **dmatrix(unsigned long n, unsigned long m) {
  
  double   **dm=NULL;
  int      i=0;

  /* Allocate array of pointers */
  if (!(dm=(double **)malloc((size_t)(n*sizeof(double *))))) {
    nferrormsg("dmatrix(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize rows */
  if ((*dm=darray(n*m))==NULL) {
    nferrormsg("dmatrix(): Cannot allocate memory for array of size %dx%d",n,m);
    return NULL;
  }

  /* Set pointers to rows */
  for (i=1; i<n; i++) dm[i]=dm[i-1]+m;

  return dm;

}


/****************************************************************************
* Allocate an array of double precision with length n. Function returns a
* pointer to first element
****************************************************************************/

float *farray(unsigned long n) {
  
  float    *fa=NULL;
  int      i=0;

  if (!(fa=(float *)malloc((size_t)(n*sizeof(float))))) {
    nferrormsg("farray(): Cannot allocate memory for array");
    return NULL;
  }

  /* Initialize array */
  for (i=0; i<n; i++) *(fa+i)=0.0;

  return fa;

}


/****************************************************************************
* Allocate a cube of floating point with dimensions n x m x q. Function
* returns a pointer to (0,0,0) element
****************************************************************************/

float ***fcube(unsigned long n, unsigned long m, unsigned long p) {
  
  float   ***fc=NULL;
  int     i=0,j=0;

  /* Allocate array of pointers to pointers */
  if (!(fc=(float ***)malloc((size_t)(n*sizeof(float **))))) {
    nferrormsg("fcube(): Cannot allocate memory for pointer\n\
\tto pointer array");
    return NULL;
  }

  /* Allocate array of pointers */
  if (!(*fc=(float **)malloc((size_t)(n*m*sizeof(float *))))) {
    nferrormsg("fcube(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize block of memory */
  if ((**fc=farray(n*m*p))==NULL) {
    nferrormsg("fcube(): Cannot allocate memory for array\n\
\tof size %dx%dx%d",n,m,p);
    return NULL;
  }

  /* Set pointers to pointers and pointers to rows */
  for (i=0; i<n; i++) {
    if (i) { fc[i]=fc[i-1]+m; fc[i][0]=fc[i-1][m-1]+p; }
    for (j=1; j<m; j++) fc[i][j]=fc[i][j-1]+p;
  }

  return fc;

}


/****************************************************************************
* Allocate a matrix of floating point with dimensions n x m. Function
* returns a pointer to (0,0) element
****************************************************************************/

float **fmatrix(unsigned long n, unsigned long m) {
  
  float   **fm=NULL;
  int      i=0;

  /* Allocate array of pointers */
  if (!(fm=(float **)malloc((size_t)(n*sizeof(float *))))) {
    nferrormsg("fmatrix(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize rows */
  if ((*fm=farray(n*m))==NULL) {
    nferrormsg("fmatrix(): Cannot allocate memory for array of size %dx%d",n,m);
    return NULL;
  }

  /* Set pointers to rows */
  for (i=1; i<n; i++) fm[i]=fm[i-1]+m;

  return fm;

}


/****************************************************************************
* Allocate an array of integers with length n. Function returns a pointer
* to first element
****************************************************************************/

int *iarray(unsigned long n) {
  
  int   *ia=NULL;
  int   i=0;

  if (!(ia=(int *)malloc((size_t)(n*sizeof(int))))) {
    nferrormsg("iarray(): Cannot allocate memory for array");
    return NULL;
  }

  /* Initialize array */
  for (i=0; i<n; i++) *(ia+i)=0;

  return ia;

}


/****************************************************************************
* Allocate a cube of integers with dimensions n x m x q. Function
* returns a pointer to (0,0,0) element
****************************************************************************/

int ***icube(unsigned long n, unsigned long m, unsigned long p) {
  
  int i=0,j=0;
  int ***ic=NULL;

  /* Allocate array of pointers to pointers */
  if (!(ic=(int ***)malloc((size_t)(n*sizeof(int **))))) {
    nferrormsg("icube(): Cannot allocate memory for pointer\n\
\tto pointer array");
    return NULL;
  }

  /* Allocate array of pointers */
  if (!(*ic=(int **)malloc((size_t)(n*m*sizeof(int *))))) {
    nferrormsg("icube(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize block of memory */
  if ((**ic=iarray(n*m*p))==NULL) {
    nferrormsg("icube(): Cannot allocate memory for array\n\
\tof size %dx%dx%d",n,m,p);
    return NULL;
  }

  /* Set pointers to pointers and pointers to rows */
  for (i=0; i<n; i++) {
    if (i) { ic[i]=ic[i-1]+m; ic[i][0]=ic[i-1][m-1]+p; }
    for (j=1; j<m; j++) ic[i][j]=ic[i][j-1]+p;
  }

  return ic;

}


/****************************************************************************
* Allocate a matrix of integers with dimensions n x m. Function
* returns a pointer to (0,0) element
****************************************************************************/

int **imatrix(unsigned long n, unsigned long m) {
  
  int   **im=NULL;
  int   i=0;

  /* Allocate array of pointers */
  if (!(im=(int **)malloc((size_t)(n*sizeof(int *))))) {
    nferrormsg("imatrix(): Cannot allocate memory for pointer array");
    return NULL;
  }

  /* Allocate and initialize rows */
  if ((*im=iarray(n*m))==NULL) {
    nferrormsg("imatrix(): Cannot allocate memory for array of size %dx%d",n,m);
    return NULL;
  }

  /* Set pointers to rows */
  for (i=1; i<n; i++) im[i]=im[i-1]+m;

  return im;

}
