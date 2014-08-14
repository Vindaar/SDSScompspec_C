/****************************************************************************
* Correct the spectra and the magnitudes for Galactic reddening. This
* relies on having an E(B-V) value from SDSS_Ebv (i.e. the Schlegel et
* al. dust maps) and the Pei (1992) extinction curve fitting formulas in
* dust_extinc.c
****************************************************************************/

#include <math.h>
#include "SDSS-spec.h"
#include "error.h"

int SDSS_Gal_extinct_correct(spectrum *spec) {

  double   R_V=3.08,A_V=0.0,A_lam=0.0,redfac=0.0;
  int      i=0;

  /* First correct the magnitudes using the relations between E(B-V) and
     A_j for each of the j=u',g',r',i' and z' SDSS bands from Schneider et
     al., (2003, AJ, 126, 2579) */
  /* Only do this if no additional info was provided in the initial
     input file by the user. Thus, we have to use only the information
     in the header, which contains no info on dust corrections */
  spec->mag[0]-=(spec->ext[0]=5.155*spec->Ebv);
  spec->mag[1]-=(spec->ext[1]=3.793*spec->Ebv);
  spec->mag[2]-=(spec->ext[2]=2.751*spec->Ebv);
  spec->mag[3]-=(spec->ext[3]=2.086*spec->Ebv);
  spec->mag[4]-=(spec->ext[4]=1/479*spec->Ebv);


  /* For each pixel, find the extinction and correct the spectrum */
  A_V=R_V*spec->Ebv;
  printf("flux element: %lf\n", spec->wl[1]);
  for (i=0; i<spec->np; i++) {
    if ((A_lam=dust_extinct(spec->wl[i],R_V,0))==-1.0)
      errormsg("SDSS_Gal_extinct_correct(): Error calculating dust extinction\n\
\tfor wavelength %lf in file\n\t%s",spec->wl[i],spec->file);
    	A_lam*=A_V;
    	redfac=pow(10.0,0.4*A_lam);
    	spec->fl[i]*=redfac;
    	spec->er[i]*=redfac;
  }

  return 1;

}
