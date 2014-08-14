#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SDSS-restspec.h"


#define SPECSWAP(a,b) {tspec=a;a=b;b=tspec;}

typedef struct Pair {
  int idx;
  double value;
} pair;

int comp_pair(const void *x1, const void *x2) {
  if( ((pair *)x1)->value > ((pair *)x2)->value) return 1;
  else if (((pair *)x1)->value == ((pair *)x2)->value) return 0;
  else return -1;
}

int comp_red(const void *x1, const void *x2) {
  if( ((spectrum *)x1)->zem > ((spectrum *)x2)->zem) return 1;
  else if (((spectrum *)x1)->zem == ((spectrum *)x2)->zem) return 0;
  else return -1;
}

int SDSS_sortredshift(spectrum *spec,int nspec) {

  int i=0,*nqso=NULL;
  int comp_pair(const void *x1, const void *x2);
  char comment[32]="\0";
  spectrum *pspec=NULL,tspec;
  pair *dpair=NULL;
  extern FILE *stream;
  extern setblock set;

  /* FITS DECLERATIONS */
  int fitstatus=0;
  fitsfile *infits=NULL;

 
  /* Read in Redshifts */
  for(i=0;i<nspec;i++) {

    /* Just the completion counter */
    if( rint(100.0*i/nspec) > rint(100.0*(i-1)/nspec))
      fprintf(stream,"\b\b\b%2i%%",(int) 100*i/nspec);

    fits_open_file(&infits,spec[i].file,READONLY,&fitstatus);
    fits_read_key(infits,TDOUBLE,"Z",&(spec[i].zem),comment,&fitstatus);
    fits_read_key(infits,TDOUBLE,"CRVAL1",&(spec[i].beginwl),comment,&fitstatus);
    fits_read_key(infits,TDOUBLE,"CD1_1",&(spec[i].deltwl),comment,&fitstatus);
    fits_read_key(infits,TINT,"NAXIS1",&(spec[i].np),comment,&fitstatus);
    fits_close_file(infits,&fitstatus);
    if(fitstatus) {
      fits_get_errstatus(fitstatus,comment);
      errormsg("%s",comment);
    }
  }

  fprintf(stream,"\b\b\b");
  /* Then sort in redshift */
  qsort(spec,nspec,sizeof(spectrum),comp_red);

  if(set.optsort==0)
    return 0;
  else if(set.optsort<0) {
    for(i=0;i<nspec;i++) {
      SPECSWAP(spec[i],spec[nspec-1-i]);
    }
  }
  else {
    /* Allocate Memory */
    dpair=(pair *) malloc(nspec*sizeof(pair));
    pspec=(spectrum *)malloc(nspec*sizeof(spectrum));
    /* Initialize data pairs for sorting */
    for(i=0;i<nspec;i++) {
      dpair[i].idx=i;
      if(set.optsort>0){
	pspec[i]=spec[i];
	dpair[i].value=fabs(set.optsort-spec[i].zem);
      }
    }
    
    qsort(dpair,nspec,sizeof(pair),comp_pair);
    for(i=0;i<nspec;i++){
      spec[i]=pspec[dpair[i].idx];
    }
    free(pspec);free(dpair);
  }

  return 0;
}
