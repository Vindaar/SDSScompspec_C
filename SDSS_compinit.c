#include <stdlib.h>
#include <stdio.h>
#include "SDSS-spec.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int SDSS_compinit(compspec **cspec, setblock *set) {


  int i=0,j=0,N=0;

  // fixed value of 3000? Should at least be as long as in each spectrum?!
  //N=3000;
  // TODO: Find proper way to change N to value of nspec
  N = 6000;

  //printf("Allocating memory for array of set->nforest composite spectra\n");
  /* Allocate Memory for array of set->nforest composite spectra */
  if(((*cspec)=(compspec *)malloc((size_t)(set->nforest)*sizeof(compspec)))
     ==NULL)
    errormsg("Could not allocate Memory for %i composite Spectra",
	     set->nforest);

  //printf("Did the array for the forest, now for the structures...\n");
  //printf("  the array is %d elements big\n", set->nforest);
  
  
  /* Allocate memory for the arrays in these structures */
  for(j=0;j<set->nforest;j++) {
    
    (*cspec)[j].wl=darray(N);			// alternatively: (*cspec)[j].wl=(double *) malloc(N*sizeof(double));
    if((*cspec)[j].wl==NULL)
    {
       printf("Could not build comp wl array in SDSS_compspec");
       exit(0);
    }
    (*cspec)[j].sum=darray(N);
    if((*cspec)[j].sum==NULL)
    {
       printf("Could not build comp sum array in SDSS_compspec");
       exit(0);
    }
    (*cspec)[j].sum2=darray(N);
    if((*cspec)[j].sum2==NULL)
    {
       printf("Could not build comp sum2 array in SDSS_compspec");
       exit(0);
    }
    (*cspec)[j].wsum=darray(N);
    if((*cspec)[j].wsum==NULL)
    {
       printf("Could not build comp wsum array in SDSS_compspec");
       exit(0);
    }
    (*cspec)[j].wsum2=darray(N);
    if((*cspec)[j].wsum2==NULL)
    {
       printf("Could not build comp wsum2 array in SDSS_compspec");
       exit(0);
    }
    (*cspec)[j].nhist=iarray(N); //was darray 'assignment from incompatible pointer type [enabled by default]'. makes sense, nhist is *int !
    if((*cspec)[j].nhist==NULL)
    {
       printf("Could not build comp nhist array in SDSS_compspec");
       exit(0);
    }
    (*cspec)[j].fl=darray(N);
    if((*cspec)[j].fl==NULL)
    {
       printf("Could not build comp fl array in SDSS_compspec");
       exit(0);
    }
    (*cspec)[j].er=darray(N);
    if((*cspec)[j].er==NULL)
    {
       printf("Could not build comp er array in SDSS_compspec");
       exit(0);
    }

    
    // printf(" Allocated memory for cspec, now let's put in some wavelengths\n");
    /* Initialize the wavelength spectra to a standard wl array */
    for(i=0;i<N;i++)
    {
      //(*cspec)[j].wl[i]=pow(10,3.57520+.0001*i);
    	(*cspec)[j].wl[i]=pow(10,3.58020+.0001*i);
    }
  }
  
  return 1;
}
