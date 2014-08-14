#include <stdlib.h>
#include "SDSS-spec.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int SDSS_compspec(compspec *cspec, spectrum *spec, setblock *set) {

// TODO: PROBLEM: If cspec.wl[0] < spec[i].wl[0] the program causes big problems!
	// Therefore cspec.wl redefined in SDSS_compinit such that it's always
	// bigger than all the beginning wavelengths of the spectra

  int i=0,bin=0;
  int goodpix=0;
  int forest=set->forest;
  int cindent=0,indent=0;
  static int count=0;
  double y=0,w=0,e=0;
  FILE *params_iter = NULL;
  params_iter = fopen("params_iter.txt", "a");

  //  double weight[9]=
  //    {1.21,5.30,5.90,5.43,3.84,2.64,2.53,1.78,1.13};
    //    {1.21,3.30,5.90,7.43,6.84,3.64,1.53,1.18,1.13};
    //    {1.61,3.30,5.90,5.43,3.84,2.94,3.53,4.18,2.73}; // derived
  double Restrange[][2]=RESTRANGE;

  if (set->debug == 1){
    printf("Address of cspec wl: %p\n\n", cspec[0].wl);
  }

  printf("... made it into SDSS_compspec\n");

  for(forest=0;forest<set->nforest;forest++) {

    /* Set the bounds on this QSO - what range we will use */
    spec->lambmin =
      MAX((spec->zem+1)*(Restrange[forest][0]),spec->wl[0]);
    spec->lambmax= (spec->zem+1)*(Restrange[forest][1]);
    /* Put in an artificial bound */
  //     spec->lambmin = 1140*(spec->zem+1);


    /* Make sure some pixels will be used, or else skip to the next QSO */
    if(spec->lambmax < spec->wl[0]){
      if(!forest) return 9;
      else return 1;
    }
    i=0;
    printf("Set bounds on range - using MAX()\n");

    if(log10(cspec[forest].wl[0])<spec->beginwl) {
      cindent=rint((spec->beginwl-log10(cspec[forest].wl[0]))/.0001);
      if (cindent > 0){
    	  fprintf(params_iter, "Address of cspec wl and cindent and forest Before: %p\t%i\t%i\n\n", cspec[0].wl, cindent, forest);
      }
      cspec[forest].wl+=cindent;
      cspec[forest].fl+=cindent;
      cspec[forest].er+=cindent;
      cspec[forest].nhist+=cindent;
      cspec[forest].sum+=cindent;
      cspec[forest].sum2+=cindent;
      cspec[forest].wsum+=cindent;
      cspec[forest].wsum2+=cindent;
    }
    else {
      indent=rint((log10(cspec[forest].wl[0])-spec->beginwl)/.0001);
      spec->wl+=indent;
      spec->fl+=indent;
      spec->er+=indent;
      spec->pc+=indent;
    }
    printf("Set log scale (called rint)\n");
    fprintf(params_iter, "%i\t%i\t%s\n", indent, cindent, spec->file);

    for(i=1;i<(spec->np - 1);i++) {
      if(abs(rint(log10(cspec[forest].wl[i])/.0001)
	     - rint(log10(spec->wl[i])/.0001))>=1) {
	     errormsg("Misalignment Error. Incompatible Wavelength Scale in file\n\t%s",spec->file);
         return 9;
         }
         
      goodpix = 
        (((spec->wl[i] !=0) &&
	  (spec->wl[i] + spec->wl[i-1] > 2*spec->lambmin) && 
	  (spec->wl[i] + spec->wl[i+1] < 2*spec->lambmax) &&
	  /* To select parts of the forest region */
	  //	  (spec->wl[i]>1110*(spec->zem+1) && spec->wl[i]<1135*(spec->zem+1)) &&
	  //	  (spec->wl[i]>1122.5*(spec->zem+1)) &&
	  (!isnan(spec->fl[i]) && !isnan(spec->er[i])) &&
	  (spec->er[i]!=0)) ? 1 : 0);    

      if(goodpix) {


        /* UNweighted */
	cspec[forest].sum[i] += spec->fl[i]/spec->pc[i]; 
	cspec[forest].sum2[i] += pow(spec->fl[i]/spec->pc[i],2);

        /* Weighted */
	//y=spec->fl[i]/spec->pc[i];
        //e=spec->er[i]/spec->pc[i]+.1;
        //w=1/pow(e,2);

        //cspec[forest].sum[i] += w*y;
        //cspec[forest].sum2[i] += w*pow(y,2);
        //cspec[forest].wsum[i] += w;
        //cspec[forest].wsum2[i] += pow(w,2);

	/* Funny Weighting 
	y=spec->fl[i]/spec->pc[i];
	if((spec->mag[0]-spec->mag[1]) < 1.5 &&
	   (spec->mag[0]-spec->mag[1]) > .6 &&
	   (spec->mag[1]-spec->mag[2]) > 0 &&
	   (spec->mag[1]-spec->mag[2]) < .2) {
	  bin=(int) ((spec->mag[0]-spec->mag[1] - .6)/.1);
	  if(bin>=9) errormsg("bintoohigh");
	  w=weight[bin];
	}
	else { 
	  w=1;
	}
	cspec[forest].sum[i]+=w*y;
	cspec[forest].wsum[i]+=w;
	cspec[forest].sum2[i]+=y*y;
	*/

        /* Count that point */
        cspec[forest].nhist[i]++;
      }

      if(spec->wl[i]+spec->wl[i+1] > 2*spec->lambmax) break;
    }
  }
  forest = set->forest;

 // printf("Address of cspec wl and cindent and forest2: %p\t%i\t%i\n\n", cspec[0].wl, cindent, forest);



  if (cindent > 0){
	  fprintf(params_iter, "Address of cspec wl and cindent and forest Between: %p\t%i\t%i\n\n", cspec[0].wl, cindent, forest);
  }


  cspec[forest].wl-=cindent;
 // printf("Address of cspec wl4: %p\n\n", cspec[0].wl);

  cspec[forest].fl-=cindent;

  if (cindent > 0){
	  fprintf(params_iter, "Address of cspec wl and cindent and forest After: %p\t%i\t%i\n\n", cspec[0].wl, cindent, forest);
  }

  cspec[forest].er-=cindent;
  cspec[forest].sum-=cindent;
  cspec[forest].sum2-=cindent;
  cspec[forest].wsum-=cindent;
  cspec[forest].wsum2-=cindent;
  cspec[forest].nhist-=cindent;
  spec->wl-=indent;
  spec->fl-=indent;
  spec->er-=indent;
  spec->pc-=indent;

  fclose(params_iter);
  return 1;
}
  






