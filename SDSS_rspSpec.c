/****************************************************************************
* Read in relevant SDSS FITS file header information and spectrum
* flag adinfo = 1 if some information is available in file adinfo and is thus not
* required from the FITS file header.
****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SDSS-spec.h"
#include "astron.h"
#include "memory.h"
#include "error.h"

int SDSS_rspSpec(spectrum *spec, setblock *set) {

  double   ra_sf=0.0,dec_sf=0.0;
  double   nulval=0.0;
  long     naxes[9] = {0,0,0,0,0,0,0,0,0};
  int      ra_h=0,ra_m=0,ra_s=0;
  int      dec_d=0,dec_m=0,dec_s=0,dec_sign=1;
  int      hdutype=0,hdunum=0,status=0,bitpix=0,first=1,naxis=0,anynul=0;
  int      i=0;
  char     beginwl_key[FLEN_KEYWORD]="CRVAL1";
  char     cpix_key[FLEN_KEYWORD]="CRPIX1";
  char     deltwl_key[FLEN_KEYWORD]="CD1_1";
  char     comment[FLEN_COMMENT]="\0";
  char     magstr[FLEN_VALUE]="\0";
  fitsfile *infits=NULL;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,spec->file,READONLY,&status))
    errormsg("SDSS_rspSpec(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  fits_get_hdu_type(infits,&hdutype,&status);
  if (hdutype!=IMAGE_HDU)
    errormsg("SDSS_rspSpec(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("SDSS_rspSpec(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum!=6 && hdunum!=7)
    errormsg("SDSS_rspSpec(): Number of HDUs is %d instead of %d\n\
\tin file %s",hdunum,7,spec->file);

  /* Resolving power for spectrum */
  spec->R=RESPOW;
  /* Find plate ID */
  if(fits_read_key(infits,TINT,"PLATEID",&(spec->plate),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PLATEID",spec->file);
  /* Find Fiber ID */
  if(fits_read_key(infits,TINT,"FIBERID",&(spec->fiberID),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","FIBERID",spec->file);
 

  /* Formulate object name from object RA and DEC */
  if (fits_read_key(infits,TDOUBLE,"RAOBJ",&(spec->ra),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","RAOBJ",spec->file);
  if (fits_read_key(infits,TDOUBLE,"DECOBJ",&(spec->dec),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","DECOBJ",spec->file);

  spec->ra/=15.0;
  ra_h=(int)spec->ra;
  ra_m=(int)((spec->ra-(double)ra_h)*60.0);
  ra_s=(int)(((spec->ra-(double)ra_h)*60.0-(double)ra_m)*60.0);
  ra_sf=((spec->ra-(double)ra_h)*60.0-(double)ra_m)*60.0;
  sprintf(spec->ra_str,"%2.2d:%2.2d:%06.3lf",ra_h,ra_m,ra_sf);
  spec->ra_str[strlen(spec->ra_str)-1]='\0';
  if (spec->dec<0.0) dec_sign=0;
  dec_d=(int)(fabs(spec->dec));
  dec_m=(int)((fabs(spec->dec)-(double)dec_d)*60.0);
  dec_s=(int)(((fabs(spec->dec)-(double)dec_d)*60.0-(double)dec_m)*60.0);
  dec_sf=((fabs(spec->dec)-(double)dec_d)*60.0-(double)dec_m)*60.0;
  if (dec_sign) sprintf(spec->dec_str,"+%2.2d:%2.2d:%05.2lf",dec_d,dec_m,dec_sf);
  else sprintf(spec->dec_str,"-%2.2d:%2.2d:%05.2lf",dec_d,dec_m,dec_sf);
  spec->dec_str[strlen(spec->dec_str)-1]='\0';
  sprintf(spec->obj,"SDSSJ%2.2d%2.2d%06.3lf",ra_h,ra_m,ra_sf);
  if (dec_sign) sprintf(&(spec->obj[strlen(spec->obj)-1]),"+%2.2d%2.2d%05.2lf",
	  dec_d,dec_m,dec_sf);
  else sprintf(&(spec->obj[strlen(spec->obj)-1]),"-%2.2d%2.2d%05.2lf",
	 dec_d,dec_m,dec_sf);
  spec->obj[strlen(spec->obj)-1]='\0';

  /* Calculate the galactic co-ords of object */
  if (!ast_eq2gal(spec->ra,spec->dec,2000.0,&(spec->l),&(spec->b)))
    errormsg("SDSS_rspSpec(): Error returned from ast_eq2gal()");

  /* Find the fiber magnitudes */
  if (fits_read_key(infits,TSTRING,"MAG",magstr,comment,&status))
	  errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","MAG",spec->file);
  if (sscanf(magstr,"%lf %lf %lf %lf %lf",&(spec->mag[0]),&(spec->mag[1]),
		  &(spec->mag[2]),&(spec->mag[3]),&(spec->mag[4]))!=5)
		  errormsg("SDSS_rspSpec(): Cannot read fibre mags. from MAG string\n\
\tin header of file\n\t%s",spec->file);


  /* Get primary target flag info and process it */
  if (fits_read_key(infits,TINT,"PRIMTARG",&(spec->tsf.primtarg),comment,
	      &status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PRIMTARG",spec->file);
  /* Convert the primtarg decimal to 1/0 flags via hex bitwise
     comparisons */
  if (!SDSS_primtarg(&(spec->tsf)))
    errormsg("SDSS_rspSpec(): Error returned from SDSS_primtarg()");


  /* Overwrite these with PSF magnitudes */
  /*if (fits_read_key(infits,TDOUBLE,"PSFUMAG",&(spec->u),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFUMAG",spec->file);
  if (fits_read_key(infits,TDOUBLE,"PSFGMAG",&(spec->g),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFGMAG",spec->file);
  if (fits_read_key(infits,TDOUBLE,"PSFRMAG",&(spec->r),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFRMAG",spec->file);
  if (fits_read_key(infits,TDOUBLE,"PSFIMAG",&(spec->i),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFIMAG",spec->file);
  if (fits_read_key(infits,TDOUBLE,"PSFZMAG",&(spec->z),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFZMAG",spec->file);*/



  /* Find SDSS collaborations estimate of redshift and redshift error  */
  if (fits_read_key(infits,TDOUBLE,"Z",&(spec->zem),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","Z",spec->file);

  /* Read in MJD of observation */
    if (fits_read_key(infits,TINT,"MJD",&(spec->mjd),comment,&status))
      errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
        \tfrom FITS file %s.","MJD",spec->file);


  /* Get image dimensions */
  if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    errormsg("SDSS_rspSpec(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis) errormsg("SDSS_rspSpec(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  else if (naxis!=2) errormsg("SDSS_rspSpec(): Abnormal number of dimensions\n\
\tin image file %s",spec->file);

  /* Find number of spectral pixels */
  if (fits_read_key(infits,TINT,"NAXIS1",&(spec->np),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","NAXIS1",spec->file);

  /* Allocate memory for wavelength, flux, error etc. arrays */
  if ((spec->wl=darray(spec->np))==NULL)
    errormsg("SDSS_rspSpec(): Cannot allocate memory for\n\
\tspec->wl array of size %d",spec->np);
  if ((spec->fl=darray(spec->np))==NULL)
    errormsg("SDSS_rspSpec(): Cannot allocate memory for\n\
\tspec->fl array of size %d",spec->np);
  if ((spec->er=darray(spec->np))==NULL)
    errormsg("SDSS_rspSpec(): Cannot allocate memory for\n\
\tspec->er array of size %d",spec->np);
  if ((spec->co=darray(spec->np))==NULL)
    errormsg("SDSS_rspSpec(): Cannot allocate memory for\n\
\tspec->co array of size %d",spec->np);
  if ((spec->pc=darray(spec->np))==NULL)
    errormsg("SDSS_rspSpec(): Cannot allocate memory for\n\
\tspec->pc array of size %d",spec->np);
  if ((spec->sn=darray(spec->np))==NULL)
    errormsg("SDSS_rspSpec(): Cannot allocate memory for\n\
\tspec->sn array of size %d",spec->np);
  if ((spec->st=iarray(spec->np))==NULL)
    errormsg("SDSS_rspSpec(): Cannot allocate memory for\n\
\tspec->st array of size %d",spec->np);

  /* Set default continuum level to Joe Liske's favourite value */
  for (i=0; i<spec->np; i++) { spec->co[i]=-10.0; spec->pc[i]=0.0; }

  /* Read in wavelength information */
  if (fits_read_key(infits,TDOUBLE,deltwl_key,&(spec->deltwl),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s!",deltwl_key,spec->file);
  if (fits_read_key(infits,TDOUBLE,cpix_key,&(spec->cpix),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s!",cpix_key,spec->file);
  if (fits_read_key(infits,TDOUBLE,beginwl_key,&(spec->beginwl),comment,&status))
    errormsg("SDSS_rspSpec(): Cannot read value of header card %s\n\
\tfrom FITS file %s!",beginwl_key,spec->file);

  /* Create wavelength scale from header cards */
  if (spec->beginwl<5.0)
    /* Logarithmic wavelength scale */
    for (i=0; i<spec->np; i++)
      spec->wl[i]=pow(10.0,
		      spec->beginwl+((double)(i+1)-spec->cpix)*spec->deltwl);
  else
    /* Linear wavelength scale */
    for (i=0; i<spec->np; i++)
      spec->wl[i]=spec->beginwl+((double)(i+1)-spec->cpix)*spec->deltwl;

  /* Record dispersion in km/s */
  spec->disp=(spec->wl[1]-spec->wl[0])/spec->wl[0]*C_C_K;

  /* Read in flux information */
  first=1;
  if (fits_read_img(infits,TDOUBLE,first,spec->np,&nulval,spec->fl,
		    &anynul,&status))
    errormsg("SDSS_rspSpec(): Cannot read flux array in FITS file\n\t%s",spec->file);
  /* Read in error information */
  first+=2*spec->np;
  if (fits_read_img(infits,TDOUBLE,first,spec->np,&nulval,spec->er,
		    &anynul,&status))
    errormsg("SDSS_rspSpec(): Cannot read error array in FITS file\n\t%s",
	     spec->file);

  /* Define signal-to-noise and status arrays */
  for (i=0; i<spec->np; i++) {
    spec->st[i]=(spec->er[i]>0.0) ? 1 : 0;
    spec->sn[i]=(spec->st[i]) ? spec->fl[i]/spec->er[i] : 0.0;
  }

  /* Close input FITS file */
  fits_close_file(infits,&status);


  return 1;

}
