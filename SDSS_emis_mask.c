/****************************************************************************
* Generate a mask array which goes to zero around emission lines . Use
* the emission lines found by SDSS for each spectrum. If the number of
* lines is equal to zero, then use the reference lines of the SDSS,
* redshifted to the redshift of the each spectrum, to define the
* positions of the lines and the regions around them to mask out. This
* is the behaviour for opt=0. If opt=1 then go straight to the second
* option above.
****************************************************************************/

#include <stdio.h>

#include "SDSS-spec.h"
#include "const.h"
#include "error.h"

int SDSS_emis_mask(spectrum *spec, emisline *emis, int nemis, int opt) {

  double  ewl=0.0;
  int     eidx=0;
  int     i=0,j=0;

  /* Initialize mask */
  for (i=0; i<spec->np; i++) spec->me[i]=1;

  /* Decide which method to use for this spectrum */
  if (spec->nem && !opt) {
    /* Use the emission lines found by SDSS in this spectrum. Note
       that these are always listed in increasing observed wavelength
       order */
    for (i=0,eidx=0,j=0; i<spec->nem; i++) {
      if (j<spec->np-1 && spec->em[i].wc>0.0 && spec->em[i].ws>0.0 &&
	  spec->em[i].we>0.0 && spec->em[i].hgt>0.0) {
	if ((j=idxdval(&(spec->wl[eidx]),spec->np-eidx,spec->em[i].ws))==-1)
	  return 1;
	j+=eidx; while (spec->wl[j]<=spec->em[i].we && j<spec->np) spec->me[j++]=0;
	if ((eidx=j)==spec->np) return 1;
      }
    }
  } else {
    /* Use the emission lines listed in the "reflines.dat" file
       provided by SDSS to define the positions of the emission
       lines. Then just mask out regions +/- vemis km/s away from the
       line centre */
    if (!nemis)
      errormsg("SDSS_emis_mask(): No emission lines have been read from\n\
\treflines.dat file from SDSS. An emission line mask cannot be\n\
\tcreated for file\n\t%s ",spec->file);
    for (i=0,eidx=0,j=0; i<nemis; i++) {
      if (j<spec->np-1 && emis[i].wc>0.0 && emis[i].ws>0.0 && emis[i].we>0.0) {
	if ((j=idxdval(&(spec->wl[eidx]),spec->np-eidx,
		       (1.0+spec->zem)*emis[i].ws))==-1)
	  return 1;
	j+=eidx; ewl=(1.0+spec->zem)*emis[i].we;
	while (spec->wl[j]<=ewl && j<spec->np) spec->me[j++]=0;
	if ((eidx=j)==spec->np) return 1;
      }
    }
  }

  return 1;

}
