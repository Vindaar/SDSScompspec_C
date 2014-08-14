/****************************************************************************
* Form goemetric mean composite spectrum from input spectra. The
* spectra to contribute to the composite are given by their indices in
* the comb array and their redshifts are specified in zcomb. Note that
* because SDSS spectra are always on the same wavelength scale as each
* other, this wavelength scale is simply shifted instead of a new (in
* this case equally arbitrary) wavelenegth scale being
* created. Moreover, in this routine we shift the spectra by the
* redshift with an accuracy of only 1 pixel so that we do not have to
* rebin the spectra. Finally, many of the operations in this routine
* could be done quicker, but here we are concerned with required
* memory rather than speed. If opt=1 then the routine also makes
* composite reddening spectra from the reddening spectra contained in
* each contributing spectrum.
*
* Note that each spectrum is normalized by the value of the power-law
* continuum [derived in SDSS_red_powerlaw()] at 1360A in the restframe
* (as defined by the redshift given in the zcomb array).
****************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "SDSS-spec.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

#define FREEALL free(swlr); free(ewlr); free(norm); free(dat); \
                free(*ldat); free(ldat); free(sidx); free(eidx);

int SDSS_composite(spectrum *spec, double *zcomb, int *comb, int ncomb,
		   spectrum *cmp, int opt) {

  double  cpix=0.0,deltwl=0.0,normwl=1360.0;
  double  maxwlr=0.0;
  double  *swlr=NULL,*ewlr=NULL,*norm=NULL,*dat=NULL,**ldat=NULL;
  int     ndat=0,ncomp=0;
  int     i=0,j=0,k=0,l=0;
  int     *sidx=NULL,*eidx=NULL;
  statset stat;

  /* Assume that deltwl and cpix are the same for all spectra */
  deltwl=spec[0].deltwl; cpix=spec[0].cpix;
  cmp->deltwl=deltwl; cmp->cpix=cpix;

  /* Allocate memory for rest frame start and end wavelength arrays
     and the normalization array */
  if ((swlr=darray(ncomb))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for swlr\n\
\tarray of size %d",ncomb);
  if ((ewlr=darray(ncomb))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for ewlr\n\
\tarray of size %d",ncomb);
  if ((norm=darray(ncomb))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for norm\n\
\tarray of size %d",ncomb);

  /* Allocate memory for start and end index arrays */
  if ((sidx=iarray(ncomb))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for sidx\n\
\tarray of size %d",ncomb);
  if ((eidx=iarray(ncomb))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for eidx\n\
\tarray of size %d",ncomb);

  /* Determine the quantized redshift for each spectrum and use to
     determine the minimum and maximum rest wavelengths covered */
  for (i=0,cmp->beginwl=INFIN,maxwlr=0.0; i<ncomb; i++) {
    swlr[i]=spec[comb[i]].beginwl-
      deltwl*(double)((int)(0.5+log10(1.0+zcomb[i])/deltwl));
    ewlr[i]=swlr[i]+deltwl*((double)spec[comb[i]].np-cpix);
    cmp->beginwl=MIN(cmp->beginwl,swlr[i]); maxwlr=MAX(maxwlr,ewlr[i]);
  }
  /* Determine number of pixels in composite spectrum */
  cmp->np=(int)(0.1+(maxwlr-cmp->beginwl)/deltwl);

  /* Allocate memory for wavelength, flux and error array of composite spectrum */
  if ((cmp->wl=darray(cmp->np))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for cmp->wl\n\
\tarray of size %d",cmp->np);
  if ((cmp->fl=darray(cmp->np))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for cmp->fl\n\
\tarray of size %d",cmp->np);
  if ((cmp->er=darray(cmp->np))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for cmp->er\n\
\tarray of size %d",cmp->np);
  /* If we also want to calculate reddening spectra, allocate
     composite reddening spectra */
  if (opt==1) {
    if ((cmp->red=dmatrix(DNRED,cmp->np))==NULL)
      errormsg("SDSS_composite(): Cannot allocate memory for cmp->er\n\
\tmatrix of size %dx%d",DNRED,cmp->np);
  }

  /* Find the starting and ending indices for each spectrum once it's
     shifted back to the rest frame */
  for (i=0; i<ncomb; i++) {
    sidx[i]=(int)(0.1+(swlr[i]-cmp->beginwl)/deltwl);
    eidx[i]=(int)(0.1+(ewlr[i]-cmp->beginwl)/deltwl);
    /* Determine the normalization for each spectrum */
    if (spec[comb[i]].ealpha)
      norm[i]=pow(10.0,spec[comb[i]].delta+spec[comb[i]].beta*
		  log10(normwl*(1.0+spec[comb[i]].zem)));
    else errormsg("SDSS_composite(): Error on spectral index is 0.0 for file\n\
\t%s.\n\tCannot determine normalization for generating composite",
		  spec[comb[i]].file);
    if (norm[i]<=0.0) errormsg("SDSS_composite(): Normalization factor for file\n\
\t%s.\n\tis <=0.0. Cannot use this normalization for generating composite",
		  spec[comb[i]].file);
  }

  /** Generate the composite by forming the geometric mean of the
      contributing pixels for each composite pixel **/
  /* If opt=1 then we must allocate more logarithmic data arrays so
     that we can generate composites reddening spectra */
  ncomp=(opt==1) ? 1+DNRED: 1;
  /* Allocate memory for temporary data array */
  if ((dat=darray(ncomb))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for dat\n\
\tarray of size %d",ncomb);
  if ((ldat=dmatrix(ncomp,ncomb))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for ldat\n\
\tmatrix of size %dx%d",ncomp,ncomb);

  for (i=0; i<cmp->np; i++) {
    /* Determine wavelength for this composite pixel */
    cmp->wl[i]=pow(10.0,cmp->beginwl+deltwl*((double)i-cpix));
    /* Loop over all potentially contributing spectra to see which
       ones really contribute to this composite pixel */
    for (j=0,ndat=0; j<ncomb; j++) {
      /* Does this spectrum contain any pixels in the right rest
	 wavelength range? */
      if ((k=i-sidx[j])>=0 && i<=eidx[j]) {
	/* Is this a valid pixel? If so, add it to data array */
	if (spec[comb[j]].er[k]>0.0 && spec[comb[j]].fl[k]>0.0) {
	  dat[ndat]=spec[comb[j]].fl[k]/norm[j];
	  ldat[0][ndat]=log10(dat[ndat]);
	  if (opt==1)
	    for (l=0; l<DNRED; l++) ldat[l+1][ndat]=log10(spec[comb[j]].red[l][k]);
	  ndat++;
	}
      }
    }
    /* Calculate the geometric mean (i.e. mean logarithm) */
    if (ndat>1) {
      if (!stats(ldat[0],ldat[0],NULL,NULL,NULL,ndat,0,&stat)) {
	nferrormsg("SDSS_composite(): Error returned from stats()"); FREEALL;
	return 0;
      }
      cmp->fl[i]=pow(10.0,stat.mean);
      if (!median(dat,NULL,ndat,&stat,0)) {
	nferrormsg("SDSS_composite(): Error returned from median()"); FREEALL;
	return 0;
      }
      cmp->er[i]=stat.siqr/sqrt((double)ndat);
      if (opt==1) {
	for (j=0; j<DNRED; j++) {
	  if (!stats(ldat[j+1],ldat[j+1],NULL,NULL,NULL,ndat,0,&stat)) {
	    nferrormsg("SDSS_composite(): Error returned from stats()"); FREEALL;
	    return 0;
	  }
	  cmp->red[j][i]=pow(10.0,stat.mean);
	}
      }
    } else { cmp->fl[i]=0.0; cmp->er[i]=-INFIN; }
  }

  /* Clean up */
  FREEALL;

  return 1;

}
