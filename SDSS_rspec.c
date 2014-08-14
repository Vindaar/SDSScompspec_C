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

int SDSS_rspec(spectrum *spec, setblock *set) {

  double   ra_sf=0.0,dec_sf=0.0;
  double   nulval=0.0;
  //long     naxes[9] = {0,0,0,0,0,0,0,0,0};
  int      ra_h=0,ra_m=0,ra_s=0;
  int      dec_d=0,dec_m=0,dec_s=0,dec_sign=1;
  int      hdutype=0,hdunum=0,status=0,anynul=0;
  //int naxis=0;
  //int bitpix=0;
  volatile int      i=0; // volatile, because it was optimized out by compiler causing problems in writing arrays
  char     beginwl_key[FLEN_KEYWORD]="COEFF0";
  // char     cpix_key[FLEN_KEYWORD]="CRPIX1";
  char     deltwl_key[FLEN_KEYWORD]="COEFF1";
  char     comment[FLEN_COMMENT]="\0";
  //char     magstr[FLEN_VALUE]="\0";
  fitsfile *infits=NULL;

  /* Open input file as FITS file */
  if (fits_open_file(&infits,spec->file,READONLY,&status))
    errormsg("SDSS_rspec(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  fits_get_hdu_type(infits,&hdutype,&status);
  if (hdutype!=IMAGE_HDU)
    errormsg("SDSS_rspec(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("SDSS_rspec(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum != 12 && hdunum != 14)
    nferrormsg("SDSS_rspec(): Number of HDUs is %d instead of %d\n\
\tin file %s",hdunum,14,spec->file);

  /* Resolving power for spectrum */
  spec->R=RESPOW;
  /* Find plate ID */
  if(fits_read_key(infits,TINT,"PLATEID",&(spec->plate),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PLATEID",spec->file);
  /* Find Fiber ID */
  if(fits_read_key(infits,TINT,"FIBERID",&(spec->fiberID),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","FIBERID",spec->file);

  // TODO: Read in MJD from HDU 1 instead of HDU 3.
 
  /* Formulate object name from object RA and DEC */
  // RAOBJ and DECOBJ seems to have been renamed to PLUG_RA and PLUG_DEC in DR10.
  // TODO: Check if that is really correct!!

  if(fits_read_key(infits,TDOUBLE,"PLUG_RA",&(spec->ra),comment,&status))
  {
	  fits_report_error(stderr, status);
	  errormsg("SDSS_rspec(): Cannot read value of header card %s\n\\tfrom FITS file %s.","PLUG_RA",spec->file);
  }

  if(fits_read_key(infits,TDOUBLE,"PLUG_DEC",&(spec->dec),comment,&status))
  {
  	  fits_report_error(stderr, status);
	  errormsg("SDSS_rspec(): Cannot read value of header card %s\n\\tfrom FITS file %s.","PLUG_DEC",spec->file);
  }

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
    errormsg("SDSS_rspec(): Error returned from ast_eq2gal()");

  /* Find the fiber magnitudes */

  /* Overwrite these with PSF magnitudes */
  // Move to HDU 3, since it contains PSFMAG and PRIMTARG
  // TODO: Check hdutype

  int test=0;
  fits_get_hdu_num(infits, &test);
  if(set->debug) printf("hdunum?: %d\n", test);
  fits_movabs_hdu(infits, 3, &hdutype, &status);
  fits_get_hdu_num(infits, &test);
  if(set->debug) printf("hdunum?: %d\n", test);

  int colnum_psfmag = 0;
  int colnum_primtarg = 0;
  fits_get_colnum(infits, CASEINSEN, "PSFMAG", &colnum_psfmag, &status);
  double psfmag[5];
  fits_read_col(infits, TDOUBLE, colnum_psfmag, 1, 1, 5, &nulval, &psfmag, &anynul, &status);
  for(i = 0; i < 5; i++){
	  spec->mag[i] = psfmag[i];
	  if(set->debug) printf("psfmag: %f\n", psfmag[i]);
  }



  // Get Primtarg
  fits_get_colnum(infits, CASEINSEN, "PRIMTARGET", &colnum_primtarg, &status);
  fits_read_col(infits, TINT, colnum_primtarg, 1, 1, 1, &nulval, &(spec->tsf.primtarg), &anynul, &status);
  if (!SDSS_primtarg(&(spec->tsf)))
	  nferrormsg("SDSS_rspec(): Error returned from SDSS_primtarg()");


  /*if (fits_read_key(infits,TDOUBLE,"PSFUMAG",&(spec->u),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFUMAG",spec->file);
  if (fits_read_key(infits,TDOUBLE,"PSFGMAG",&(spec->g),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFGMAG",spec->file);
  if (fits_read_key(infits,TDOUBLE,"PSFRMAG",&(spec->r),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFRMAG",spec->file);
  if (fits_read_key(infits,TDOUBLE,"PSFIMAG",&(spec->i),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFIMAG",spec->file);
  if (fits_read_key(infits,TDOUBLE,"PSFZMAG",&(spec->z),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PSFZMAG",spec->file);*/



  /* Find SDSS collaborations estimate of redshift and redshift error  */
  int colnum_Z = 0;
  fits_get_colnum(infits, CASEINSEN, "Z", &colnum_Z, &status);
  fits_read_col(infits, TDOUBLE, colnum_Z, 1, 1, 1, &nulval, &(spec->zem), &anynul, &status);

  if(set->debug) printf("colnum_Z and Z: %d\t%f\n", colnum_Z, spec->zem);

  /*if (fits_read_key(infits,TDOUBLE,"Z",&(spec->zem),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","Z",spec->file);*/

  /* Read in MJD of observation */
  int colnum_mjd = 0;
  fits_get_colnum(infits, CASEINSEN, "MJD", &colnum_mjd, &status);
  fits_read_col(infits, TINT, colnum_mjd, 1, 1, 1, &nulval, &(spec->mjd), &anynul, &status);
    /*if (fits_read_key(infits,TINT,"MJD",&(spec->mjd),comment,&status))
      errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
        \tfrom FITS file %s.","MJD",spec->file);*/


  // Move back to HDU 2.
  fits_report_error(stderr, status);
  fits_movabs_hdu(infits, 2, &hdutype, &status);


  /* Get image dimensions */
  /*if (fits_get_img_param(infits,9,&bitpix,&naxis,naxes,&status))
    nferrormsg("SDSS_rspec(): Couldn't get image dimensions for\n\
\tfile %s",spec->file);
  if (!naxis) nferrormsg("SDSS_rspec(): Couldn't find input image extension\n\
\tfor file %s",spec->file);
  else if (naxis!=2) nferrormsg("SDSS_rspec(): Abnormal number of dimensions\n\
\tin image file %s",spec->file);*/

  /* Find number of spectral pixels */
  int statusin = 0;
  if (fits_read_key(infits,TINT,"NAXIS2",&(spec->np),comment,&statusin))
    nferrormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","NAXIS2",spec->file);

  if(set->debug) printf("spec->np: %d\n", spec->np);

  /* Allocate memory for wavelength, flux, error etc. arrays */
  if ((spec->wl=darray(spec->np))==NULL)
    errormsg("SDSS_rspec(): Cannot allocate memory for\n\
\tspec->wl array of size %d",spec->np);
  if ((spec->fl=darray(spec->np))==NULL)
    errormsg("SDSS_rspec(): Cannot allocate memory for\n\
\tspec->fl array of size %d",spec->np);
  if ((spec->er=darray(spec->np))==NULL)
    errormsg("SDSS_rspec(): Cannot allocate memory for\n\
\tspec->er array of size %d",spec->np);
  if ((spec->co=darray(spec->np))==NULL)
    errormsg("SDSS_rspec(): Cannot allocate memory for\n\
\tspec->co array of size %d",spec->np);
  if ((spec->pc=darray(spec->np))==NULL)
    errormsg("SDSS_rspec(): Cannot allocate memory for\n\
\tspec->pc array of size %d",spec->np);
  if ((spec->sn=darray(spec->np))==NULL)
    errormsg("SDSS_rspec(): Cannot allocate memory for\n\
\tspec->sn array of size %d",spec->np);
  if ((spec->st=iarray(spec->np))==NULL)
    errormsg("SDSS_rspec(): Cannot allocate memory for\n\
\tspec->st array of size %d",spec->np);

  /* Set default continuum level to Joe Liske's favourite value */
  for (i=0; i<spec->np; i++) { spec->co[i]=-10.0; spec->pc[i]=0.0; }

  /* Read in wavelength information */
  // Have to move to HDU 1 first
  fits_movabs_hdu(infits, 1, &hdutype, &status);
  if (fits_read_key(infits,TDOUBLE, deltwl_key, &(spec->deltwl),comment,&status)){
	  fits_report_error(stdout, status);
	  if(set->debug) printf("Deltawl: %f\n", spec->deltwl);
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s!",deltwl_key,spec->file);
  }
  /*if (fits_read_key(infits,TDOUBLE,cpix_key,&(spec->cpix),comment,&status)){
   nferrormsg("SDSS_rspec(): Cannot read value of header card %s\n\
\tfrom FITS file %s!",cpix_key,spec->file);
   // TODO: Check, if this makes sense.
   // If starting pixel cannot be found, assume 1.
   spec->cpix = 1;
  }*/
  spec->cpix = 1;
  if (fits_read_key(infits,TDOUBLE,beginwl_key,&(spec->beginwl),comment,&status))
    errormsg("SDSS_rspec(): Cannot read value of header card %s\n\
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
  // Flux information is stored in HDU 2 in case of spec files.

  /*if (fits_read_img(infits,TDOUBLE,first,spec->np,&nulval,spec->fl,
		    &anynul,&status))
    errormsg("SDSS_rspec(): Cannot read flux array in FITS file\n\t%s",spec->file);*/
  // Read in number of pixels. Move to HDU 2 first
  int colnum_flux = 0;
  fits_movabs_hdu(infits, 2, &hdutype, &status);

  // Read in Flux information
  fits_get_colnum(infits, CASEINSEN, "flux", &colnum_flux, &status);

  if(fits_read_col(infits, TDOUBLE, colnum_flux, 1, 1, spec->np, &nulval, spec->fl, &anynul, &status))
	  errormsg("SDSS_rspec(): Cannot read line from flux column from FITS file %s", spec->file);

  //printf("number of pixels: %d\n", spec->np);
  /* Read in error information */
  // Error given as inverse variance.
  int colnum_err = 0;
  fits_get_colnum(infits, CASEINSEN, "ivar", &colnum_err, &status);
  fits_report_error(stdout, status);

  if( fits_read_col(infits, TDOUBLE, colnum_err, 1, 1, spec->np, &nulval, spec->er, &anynul, &status)){
	  fits_report_error(stdout, status);
	  errormsg("SDSS_rspec(): Cannot read line from error column from FITS file %s", spec->file);
  }

  // Since the error is given in inverse variance, we still have to convert to one STD:

  for(i = 0; i < spec->np; i++){
	  spec->er[i] = sqrt(1 / spec->er[i]);
  }

/*  if (fits_read_img(infits,TDOUBLE,first,spec->np,&nulval,spec->er,
		    &anynul,&status))
    errormsg("SDSS_rspec(): Cannot read error array in FITS file\n\t%s",
	     spec->file);*/


  /* Define signal-to-noise and status arrays */
  for (i=0; i<spec->np; i++) {
    spec->st[i]=(spec->er[i]>0.0) ? 1 : 0;
    spec->er[i]=spec->st[i] ? pow(spec->er[i],-.5) : 0;
    spec->sn[i]=(spec->st[i]) ? spec->fl[i]/spec->er[i] : 0.0;
  }

  /* Close input FITS file */
  fits_close_file(infits,&status);

  return 1;

}
