/****************************************************************************
* Fit a power-law continuum to the red part (i.e. longer than Lya
* emission) of a SDSS QSO spectrum. The routine derives the median
* flux within a series of small portions of spectrum and then uses
* those portions to determine the spectral index.
****************************************************************************/

#include <math.h>
#include "SDSS-spec.h"
#include "fit.h"
#include "stats.h"
#include "const.h"
#include "memory.h"
#include "error.h"

#define FREE_SPEC free(coeff); free(w); free(f); free(e);
#define FREE_WORK free(work_a1); free(*work_m1); free(work_m1); free(*work_m2); \
                  free(work_m2); free(*work_m3); free(work_m3);

int SDSS_red_powerlaw(spectrum *spec) {

  static double emfree[4][2] =
    {{1280.0,1292.0},{1312.0,1328.0},{1345.0,1365.0},{1440.0,1475.0}};
    //     {1680.0,1700.0}};//,{1968.0,1982.0},{2020.0,2040.0},{2150.0,2170.0},
		     // {2190.0,2250.0}}; 
  //{{1280.0,1292.0},{1320.0,1328.0},{1345.0,1358.0},{1440.0,1475.0},{1683.0,1702.0}};

  double  chisq=0.0;
  double  *w=NULL,*f=NULL,*e=NULL;
  double  *coeff=NULL,*work_a1=NULL;
  double  **work_m1=NULL,**work_m2=NULL,**work_m3=NULL;
  int     idxswl=0,idxewl=0,nmed=0;

  //int 	  itemp=0;
  int     i=0,j=0,k=0;
  //static int miss=0;
  statset stat;

  /* Initialize spectral indices and errors */
  spec->beta=spec->alpha=spec->ealpha=spec->delta=-999.0;

  /* Allocate memory for temporary spectrum arrays */
  if ((coeff=darray(2))==NULL)
    errormsg("SDSS_red_powerlaw(): Cannot allocate memory to coeff\n\
\tarray of size %d",2);
  if ((w=darray(spec->np))==NULL)
    errormsg("SDSS_red_powerlaw(): Could not allocate memory\n\
\tfor w array of size %d",spec->np);
  if ((f=darray(spec->np))==NULL)
    errormsg("SDSS_red_powerlaw(): Could not allocate memory\n\
\tfor f array of size %d",spec->np);
  if ((e=darray(spec->np))==NULL)
    errormsg("SDSS_red_powerlaw(): Could not allocate memory\n\
\tfor e array of size %d",spec->np);

  /* Initialize spectrum's emfree matrix */
  for (i=0; i<NEMFREE; i++) {
    spec->emfree[0][i]=(1.0+spec->zem)*0.5*(emfree[i][0]+emfree[i][1]);
    spec->emfree[1][i]=0.0;
    spec->emfree[2][i]=-1.0;
  }
  /** Find the medians within the emission-free regions **/
  for (i=0,j=0; i<NEMFREE; i++) {
    /* Find first pixel with rest-frame wavelength longer than start
       of emission-free region */
    if ((idxswl=idxdval(&(spec->wl[j]),spec->np-j,(1.0+spec->zem)*emfree[i][0]))
	==-1) break;
    idxswl+=j;
    /* Find last pixel with rest-frame wavelength shorter than end
       of emission-free region */
    if ((idxewl=idxdval(&(spec->wl[j]),spec->np-j,(1.0+spec->zem)*emfree[i][1]))
	==-1) idxewl=spec->np-1;
    else idxewl+=j-1;
    j=idxewl;
    // printf("idxs: %i\t%i\n", idxswl, idxewl);
    /* Determine how many pixels can be used for median */
    for (k=idxswl,nmed=0; k<=idxewl; k++)
      if (spec->er[k]>0.0) f[nmed++]=spec->fl[k];
    if (nmed>2) {
        //printf("idx %i\t%i\n", idxswl, idxewl);
        //int n = 0;
        //for(n = idxswl; n < idxewl; n++){
//        	if (i==3){
        		//printf("%f\tflux\t%f\n", spec->fl[n], spec->er[n]);
        	//}
        //}
        if (!median(f,NULL,nmed,&stat,0)) {
	nferrormsg("SDSS_red_powerlaw(): Error calculating median of\n\
\t%d pixels between pixel %d and %d (wavelengths %lf and %lf)\n\
\tin file %s",nmed,idxswl+1,idxewl+1,spec->wl[idxswl],spec->wl[idxewl],spec->file);
	FREE_SPEC;
	return 0;
      }
      spec->emfree[1][i]=stat.med;
      spec->emfree[2][i]=stat.siqr/sqrt((double)nmed);
      //printf("med and siqr and nmed: %.18lf\t%lf\t%i\n", stat.med, stat.siqr, nmed);
    }
  }

  /** Calculate the spectral index from the medians in the
      emission-free chunks **/
  /* First make sure there are 2 or more medians to work with and
     record the logarithms of wavelengths, fluxes and errors */
  for (i=0,nmed=0; i<NEMFREE; i++) {
    if (spec->emfree[1][i]>0.0 && spec->emfree[2][i]>0.0) {
      w[nmed]=log10(spec->emfree[0][i]);
      f[nmed]=log10(spec->emfree[1][i]);
      e[nmed++]=log10(1.0+spec->emfree[2][i]/spec->emfree[1][i]);
    }
  }

  //for(j=0; j<NEMFREE; j++){
//  	printf("%f\t%f\n", w[j],f[j]);
  	//printf("%f\n", spec->emfree[1][j]);
  //}

  /* If there's enough points then do a linear fit to log-log data */

  if (nmed==NEMFREE) {
    if (!svdfit(w,f,e,nmed,coeff,2,&work_m1,&work_m2,&work_a1,&chisq,
		svdfit_poly)){
    	nferrormsg("SDSS_red_powerlaw(): Error returned from svdfit() when\n\
\tfitting powerlaw to %d median values from spectrum in file\n\t%s",nmed,
	       spec->file);
    	return 0;
    }

    if (!svdvar(work_m2,2,work_a1,&work_m3))
      errormsg("SDSS_red_powerlaw(): Error returned from svdvar()\n\
\tafter fitting powerlaw to %d median values from spectrum in file\n\t%s",nmed,
	       spec->file);


    /* Record spectral index and error */
    spec->beta=coeff[1];
    spec->alpha=-spec->beta-2.0;
    //spec->alpha = spec->beta;
    spec->delta=coeff[0];
    spec->ealpha=sqrt(work_m3[1][1]);
    spec->chisq = chisq;
    // freeing was done here before
    free(*work_m1);
    free(work_m1);
    free(*work_m2);
    free(work_m2);
    free(work_a1);
    free(*work_m3);
    free(work_m3);
  }
  //  else 
  //    fprintf(stderr,"%i\n",++miss);

  /* Copy power-law to continuum */
  for (i=0; i<spec->np; i++){
    spec->pc[i]=pow(10.0,coeff[0]+coeff[1]*log10(spec->wl[i]));
  }

  //printf("alpha, delta: %lf\t%lf\t%lf\n", spec->alpha, spec->delta, work_m3[1][1]);

  /*REMOVE THIS *
  double xtemp=0,a=-1.71;
  for (i=0,xtemp=0; i<nmed;i++) {
    xtemp+=(f[i]-a*w[i])/nmed;
  }

  for (i=0; i<spec->np; i++)
    spec->pc[i]=pow(10.0,xtemp+a*log10(spec->wl[i]));
  spec->beta=-1.56; spec->alpha=-spec->beta-2.0;
  spec->delta=xtemp;spec->ealpha=0;
  *Remove till here */


  /* Clean up and return */

  // TODO: Try to free memory from here!
  free(coeff);
  free(w);
  free(f);
  free(e);

  if (spec->alpha == -999){
	  printf("\n\n\n hey, hey\n\n\n");
  }

  return 1;

}
