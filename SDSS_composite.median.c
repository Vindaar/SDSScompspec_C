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
* memory rather than speed.
*
* Note that each spectrum is normalized by the median flux in the 51
* pixels around 1360A in the restframe (as defined by the redshift
* given in the zcomb array). If this isn't possible, a warning is
* issued and we use the power-law fit parameters derived in
* SDSS_red_powerlaw() to determine the flux level at 1360A. If the
* power-law fit parameters are not valid, an error message is
* returned.
****************************************************************************/

#include <stdio.h>

#include <stdlib.h>
#include <math.h>
#include "SDSS-spec.h"
#include "stats.h"
#include "memory.h"
#include "error.h"

#define FREEALL free(swlr); free(ewlr); free(norm); free(dat); free(ldat); \
                free(sidx); free(eidx);

int SDSS_composite(spectrum *spec, double *zcomb, int *comb, int ncomb,
		   spectrum *cmp) {

  double  cpix=0.0,deltwl=0.0,normwl=1360.0;
  double  maxwlr=0.0;
  double  *swlr=NULL,*ewlr=NULL,*norm=NULL,*dat=NULL,*med=NULL;
  int     medpix=51,ndat=0,nmed=0;
  int     i=0,j=0,k=0;
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

  /* Allocate memory for wavelength and flux array of composite spectrum */
  if ((cmp->wl=darray(cmp->np))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for cmp->wl\n\
\tarray of size %d",cmp->np);
  if ((cmp->fl=darray(cmp->np))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for cmp->fl\n\
\tarray of size %d",cmp->np);

  /* Find the starting and ending indices for each spectrum once it's
     shifted back to the rest frame */
  /* Also, determine the normalization for each spectrum using 51
     pixels around 1360A in the restframe defined by the redshifts in
     zcomb. At the moment, just return an error if this is
     unsuccessful */
  /* Allocate memory for median filter array */
  if ((med=darray(medpix))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for med\n\
\tarray of size %d",medpix);
  for (i=0; i<ncomb; i++) {
    /* Find starting and ending pixels */
    sidx[i]=(int)(0.1+(swlr[i]-cmp->beginwl)/deltwl);
    eidx[i]=(int)(0.1+(ewlr[i]-cmp->beginwl)/deltwl);
    /* Determine starting and ending pixels for median filter */
    j=(int)(0.5+(log10(normwl*(1.0+zcomb[i]))-spec[comb[i]].beginwl)/deltwl)-
      medpix/2;
    if (j<0)
      errormsg("SDSS_composite(): Normalization wavelength (%lf)\n\
\tis below starting wavelength of spectrum in file\n\t%s",normwl*(1.0+zcomb[i]),
	       spec[comb[i]].file);
    if (j+medpix-1>=spec[comb[i]].np)
      errormsg("SDSS_composite(): Normalization wavelength (%lf)\n\
\tis above ending wavelength of spectrum in file\n\t%s",
	       normwl*(1.0+zcomb[i]),spec[comb[i]].file);
    /* Load median array with data */
    for (k=j,nmed=0; k<=j+medpix-1; k++)
      if (spec[comb[i]].er[k]>0.0) med[nmed++]=spec[comb[i]].fl[k];
    /* Take median */
    if (nmed) {
      if (!median(med,NULL,nmed,&stat,0))
	errormsg("SDSS_composite(): Error returned from median()");
      norm[i]=stat.med;
    } else {
      warnmsg("SDSS_composite(): Cannot normalize spectrum in file\n\t%s.\n\
\tNo valid pixels in median region. Resort to power-law fit.",spec[comb[i]].file);
      if (spec[comb[i]].ealpha)
	norm[i]=pow(10.0,spec[comb[i]].delta+spec[comb[i]].beta*
		    log10(normwl*(1.0+zcomb[i])));
      else errormsg("SDSS_composite(): Error on spectral index is 0.0 for file\n\
\t%s.\n\tCannot determine normalization for generating composite",
		    spec[comb[i]].file);
      if (norm[i]<=0.0)
	errormsg("SDSS_composite(): Normalization factor for file\n\
\t%s.\n\tis <=0.0. Cannot use this normalization for generating composite",
		 spec[comb[i]].file);
    }
  }

  /** Generate the composite by forming the geometric mean of the
      contributing pixels for each composite pixel **/
  /* Allocate memory for temporary data array */
  if ((dat=darray(ncomb))==NULL)
    errormsg("SDSS_composite(): Cannot allocate memory for dat\n\
\tarray of size %d",ncomb);
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
	if (spec[comb[j]].er[k]>0.0 && spec[comb[j]].fl[k]>0.0)
	  dat[ndat++]=log10(spec[comb[j]].fl[k]/norm[j]);
      }
    }
    /* Calculate the geometric mean (i.e. mean logarithm) */
    if (ndat>1) {
      if (!stats(ldat,ldat,NULL,NULL,NULL,ndat,0,&stat)) {
	nferrormsg("SDSS_composite(): Error returned from stats()"); FREEALL;
	return 0;
      }
      cmp->fl[i]=pow(10.0,stat.mean);
      if (!median(dat,NULL,ndat,&stat,0)) {
	nferrormsg("SDSS_composite(): Error returned from median()"); FREEALL;
	return 0;
      }
      cmp->er[i]=stat.siqr/sqrt((double)ndat);
    } else { cmp->fl[i]=0.0; cmp->er[i]=-INFIN; }
  }

  /* Clean up */
  FREEALL;

  return 1;

}
