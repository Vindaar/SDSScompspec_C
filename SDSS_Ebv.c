/****************************************************************************
* Calculate the Galacitc dust extinction E(B-V) from the Schlegel et
* al. dust maps. This routine links to the lambert_getval() routine the
* package that comes with the dust maps.
****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "SDSS-spec.h"
#include <memory.h>
#include <error.h>
/* The following include files are from the dust maps package */
#include <interface.h>
#include <subs_fits.h>
#include <subs_lambert.h>


int SDSS_Ebv(spectrum *spec) {

  float    *l=NULL,*b=NULL;
  float    *Ebv;
  char     Nmap[LNGSTRLEN]="\0",Smap[LNGSTRLEN]="\0";
//  char     *Nmap={"../dust_maps/maps/SFD_dust_4096_ngp.fits"};
//  char     *Smap={"../dust_maps/maps/SFD_dust_4096_sgp.fits"};
  char     *cptr=NULL;

  /* First find out whether the dust maps can be found on the system */

/*  if ((cptr=getenv("DUST_MAPS"))==NULL) sprintf(Nmap,"%s",DUSTMAPS);
  else sprintf(Nmap,"%s",cptr);
  //sprintf(Smap,"%s",Nmap); strcat(Nmap,DUSTNMAP); strcat(Smap,DUSTSMAP);
  if (access(Nmap,R_OK))
    errormsg("SDSS_Ebv(): Cannot read northern dust map at\n\t%s",Nmap);
  if (access(Smap,R_OK))
  errormsg("SDSS_Ebv(): Cannot read southern dust map at\n\t%s",Smap);*/


  /* First find out whether the dust maps can be found on the system */
  if ((cptr=getenv("DUST_MAPS"))==NULL) sprintf(Nmap,"%s",DUSTMAPS);
  else sprintf(Nmap,"%s",cptr);
  sprintf(Smap,"%s",Nmap); strcat(Nmap,DUSTNMAP); strcat(Smap,DUSTSMAP);
  if (access(Nmap,R_OK))
    errormsg("SDSS_Ebv(): Cannot read northern dust map at\n\t%s",Nmap);
  if (access(Smap,R_OK))
    errormsg("SDSS_Ebv(): Cannot read southern dust map at\n\t%s",Smap);

  /* Allocate memory for galactic longitude and latitude arrays */
  if ((l=farray(1))==NULL)
    errormsg("SDSS_Ebv(): Cannot allocate memory for l\n\tarray of size %d",1);
  if ((b=farray(1))==NULL)
    errormsg("SDSS_Ebv(): Cannot allocate memory for b\n\tarray of size %d",1);

  /* Fill l and b arrays with l and bs from spectra */
  l[0]=spec->l; b[0]=spec->b;

  /* Determine E(B-V) using the Schlegel et al. package */
  Ebv=lambert_getval(Nmap,Smap,1,l,b,1,0,0);
  printf("Ebv: %f\t%f\t%f\n", Ebv[0], l[0], b[0]);

  /* Put the E(B-V) values into the structure of each spectrum */
  spec->Ebv=Ebv[0];

  /* Clean up */
  free(l);
  free(b);
  free(Ebv);
  
  return 1;

}
