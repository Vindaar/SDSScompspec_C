/****************************************************************************
* Calculate the CIV balnicity using a similar, but not the same,
* algorithm of Weymann et al., 1991, ApJ, 373, 23. We use their
* definition, except for in the following ways:
*
* 1) We do not use 1549.0 as the laboratory wavelength of the CIV
* doublet, we use the SDSS estimate which contains some 600km/s
* blueshift. Therefore, we need to start the balnicity integral
* somewhat closer to the CIV line than in the Weymann
* definition. Also, we aim to be somewhat conservative as to which
* QSOs are called BALs and so we begin the balnicity integral not at
* 3000-600 km/s below CIV but 1000km/s below the CIV emission
* redshift.
*
* 2) Since we are dealing with typically much noisier data than the
* spectra Weymann et al. considered, we cannot rely on the flux
* remaining below 75% of the continuum flux when there is a genuine
* BAL. Therefore, we smooth the data within the balnicity window
* (i.e. -25000 -- -1000 km/s w.r.t. CIV), not considering any data
* outside the window (important for the edges of the window) with a
* median filter. The median value is then used as the normalised flux
* for that pixel.
*
* 3) When significance of median filtered data is <snthres(=3.0) then
* the spectrum is judged to be too poor S/N for BAL detection
*
****************************************************************************/

#include <stdlib.h>
#include "SDSS-spec.h"
#include "stats.h"
#include "memory.h"
#include "const.h"
#include "error.h"

int SDSS_CIVbalnicity(spectrum *spec) {

  double  swl=0.0,ewl=0.0,hdisp=0.0;
  double  co=0.0,cwl=0.0,lewl=0.0,brac=0.0,frac=0.0,C=0.0,BI=0.0;
  double  snthres=3.0;
  double  *dat=NULL,*med=NULL;
  int     sidx=0,eidx=0,ndat=0,nmed=0,abs=0;
  int     i=0,j=0;
  statset stat;

  /* Initialize cbal flag and spectrum's balnicity */
  spec->cbal=BAL_REJ; spec->cbalnic=0.0;
  
  /* Half-dispersion */
  hdisp=0.5*spec->disp/C_C_K;

  /* Determine whether this spectrum has enough of the CIV line
     available to determine the balnicity */
  swl=CIV*(1.0+spec->zem)*(1.0-1000.0/C_C_K);
  ewl=CIV*(1.0+spec->zem)*(1.0-25000.0/C_C_K);
  if (swl<spec->wl[0]*(1.0-hdisp) || ewl>spec->wl[spec->np-1]*(1.0+hdisp))
    return 1;
  eidx=idxdval(spec->wl,spec->np,ewl);
  if (eidx && ewl<spec->wl[eidx-1]*(1.0+hdisp)) eidx--;
  sidx=idxdval(&(spec->wl[eidx]),spec->np-eidx,swl); sidx+=eidx;
  if (swl<spec->wl[sidx-1]*(1.0+hdisp)) sidx--;

  /* Check for 2000km/s regions of bad data */
  i=eidx; while (i<sidx) {
    if (spec->er[i]<=0.0) {
      j=i+1; while (j<=sidx && spec->er[j]<=0.0) j++;
      if ((spec->wl[j]-spec->wl[i])/spec->wl[j]*C_C_K>2000.0) {
	spec->cbal=BAL_BAD; return 1;
      } else i=j;
    } else i++;
  }
  spec->cbal=BAL_ACC;

  ndat=sidx-eidx+1;
  /* Allocate memory for median, data and status arrays */
  if ((dat=darray(ndat))==NULL)
    errormsg("SDSS_CIVbalnicity(): Cannot allocate memory to dat\n\
\tarray of size %d",ndat);
  if ((med=darray(ndat))==NULL)
    errormsg("SDSS_CIVbalnicity(): Cannot allocate memory to med\n\
\tarray of size %d",ndat);

  /* Create median filtered array */
  for (i=eidx; i<=sidx; i++) {
    co=spec->pc[i];
    if (spec->er[i]>0.0 && co) dat[i-eidx]=spec->fl[i]/co;
    else dat[i-eidx]=-1.0;
  }
  nmed=13; /* Set median filter to 13 pixels */
  if (!medianrun(dat,med,&(spec->st[eidx]),ndat,nmed)) 
    errormsg("SDSS_CIVbalnicity(): Error returned from medianrun() when\n\
\tcalculating median filtered flux");
  /* Calculate mean S/N of BAL region */
  if (!stats(&(spec->sn[eidx]),NULL,NULL,NULL,&(spec->st[eidx]),ndat,0,&stat))
    errormsg("SDSS_CIVbalnicity(): Unknown error returned from stat() when\n\
\tcalculating mean snr");
  if (stat.mean<snthres) {
    spec->cbal=BAL_SNR; spec->cbalnic=0.0; return 1;
  }
  /* Balnicity calculation */
  for (i=sidx,cwl=ewl; i>=eidx; i--) { 
    frac=1.0; brac=(spec->er[i]>0.0) ? 1.0-med[i-eidx]/0.75 : -0.1;
    if (brac>0.0 && !abs) { abs=1; cwl=spec->wl[i]*(1.0-2000.0/C_C_K); }
    else if (brac>0.0 && abs) {
      lewl=(i) ? spec->wl[i-1]*(1.0+hdisp) : spec->wl[i]*(1.0-hdisp);
      if (C==0.0 && lewl<cwl) { C=1.0; frac=(cwl-lewl)/lewl*C_C_K/spec->disp; }
      else if (C==1.0 && lewl<cwl) frac=1.0;
    }
    else { abs=0; C=0.0; cwl=ewl; }
    BI+=brac*frac*C*spec->disp;
  }
  spec->cbalnic=BI;

  /* Clean up */
  free(med); free(dat);

  return 1;

}
