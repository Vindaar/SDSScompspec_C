/******************************************************************************
AST_EQ2GAL: Convert equatorial coordinates to galactic coordinates. Input
 coordinates and epoch are changed to epoch of galactic coordinates.

double	ra		# Right ascension (hours)
double	dec		# Declination (degrees)
double	epoch		# Epoch of coordinates
double	l		# Galactic longitude (degrees)
double	b		# Galactic latitude (degrees)

******************************************************************************/

#include <math.h>
#include "astron.h"
#include "error.h"

int ast_eq2gal(double ra, double dec, double epoch, double *l, double *b) {

  double rar=0.0,decr=0.0,drar=0.0,cosdecg=0.0,sindecg=0.0,cosdecr=0.0;
  double x=0.0,y=0.0,z=0.0,r=0.0,temp=0.0;

  /* Precess the coordinates to 1950.0 */
  if (!ast_precess(ra,dec,epoch,C_GEPOCH,&rar,&decr))
    errormsg("ast_eq2gal(): Error returned from ast_precess()");

  /* Precompute the necessary constants */
  drar=(15.0*rar-C_RAGPOLE)*C_RPDEG;
  cosdecg=cos(C_DECGPOLE*C_RPDEG);
  sindecg=sin(C_DECGPOLE*C_RPDEG);
  cosdecr=cos(decr*C_RPDEG);

  /* Compute the tansformation equations */
  x=cosdecr*cos(drar); y=cosdecr*sin(drar); z=sin(decr*C_RPDEG);
  temp=z*cosdecg-x*sindecg; z=z*sindecg+x*cosdecg; x=temp;
  r=sqrt(x*x+y*y);

  /* Compute lii and bii and convert to degrees */
  if (r<2.22e-16) *l=0.0;
  else *l=C_LONGNCP*C_RPDEG+atan2(-y,x);
  if (*l<0.0) *l+=C_2PI;
  *b=atan2(z,r); *l/=C_RPDEG; *b/=C_RPDEG;

  return 1;

}
