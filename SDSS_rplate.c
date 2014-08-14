/****************************************************************************
* Read in relevant SDSS FITS file header information and spectrum
* This program written if FITS file is an spPlate file.
****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "SDSS-spec.h"
#include "astron.h"
#include "memory.h"
#include "error.h"


int SDSS_rplate(spectrum *spec, setblock *set) {


  double   ra_sf=0.0,dec_sf=0.0;
  int      ra_h=0,ra_m=0,ra_s=0;
  int      dec_d=0,dec_m=0,dec_s=0,dec_sign=1;
  int      hdutype=0,hdunum=0,status=0,loglin=0,anynul=0;
  int      colnum=0;
  int      i=0;
  long     fpixel[2]={1,spec->fiberID};
  char     beginwl_key[FLEN_KEYWORD]="CRVAL1";
  char     cpix_key[FLEN_KEYWORD]="CRPIX1";
  char     deltwl_key[FLEN_KEYWORD]="CD1_1";
  char     comment[FLEN_COMMENT]="\0";
  char     buffer[HUGESTRLEN]="\0";
  char     *cptr=NULL;
  double   *s=NULL;
  int      *t=NULL;

  fitsfile *infits=NULL;

  if(access(spec->file,F_OK)!=0)
    errormsg("Cannot Access file %s",spec->file);
  /* Open input file as FITS file */
  if (fits_open_file(&infits,spec->file,READONLY,&status))
    errormsg("SDSS_rplate(): Cannot open FITS file\n\t%s",spec->file);

  /* Check HDU type */
  fits_get_hdu_type(infits,&hdutype,&status);
  if (hdutype!=IMAGE_HDU)
    errormsg("SDSS_rplate(): File not a FITS image: %s",spec->file);

  /* Check number of HDUs */
  if (fits_get_num_hdus(infits,&hdunum,&status))
    errormsg("SDSS_rplate(): Cannot find number of HDUs in file\n\
\t%s",spec->file);
  if (hdunum!=9)
    warnmsg("SDSS_rplate(): Number of HDUs is %d instead of %d\n\
\tin file %s",hdunum,9,spec->file);

  /* Resolving power for spectrum */
  spec->R=RESPOW;
  /* Find plate ID */
  if(fits_read_key(infits,TINT,"PLATEID",&(spec->plate),comment,&status))
    errormsg("SDSS_rplate(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","PLATEID",spec->file);
    /* Read in MJD of observation */
  if (fits_read_key(infits,TINT,"MJD",&(spec->mjd),comment,&status))
    errormsg("SDSS_rplate(): Cannot read value of header card %s\n\
      \tfrom FITS file %s.","MJD",spec->file);

  printf("FiberID: %d\n", spec->fiberID);

  /* Move to HDU 6 */
  if(fits_movabs_hdu(infits,6,&hdutype,&status))
    errormsg("Could Not Move to HDU 6 in file %s",spec->file);
  if(hdutype != 2) 
    errormsg("HDU 6 not a Binary Table in file %s",spec->file);
  /* Start Reading values out of the table */
  fits_get_colnum(infits,CASEINSEN,"RA",&colnum,&status);
  fits_read_col(infits,TDOUBLE,colnum,spec->fiberID,1,1,&s,&(spec->ra),&anynul,&status);
  fits_get_colnum(infits,CASEINSEN,"DEC",&colnum,&status);
  fits_read_col(infits,TDOUBLE,colnum,spec->fiberID,1,1,&s,&(spec->dec),&anynul,&status);
  fits_get_colnum(infits,CASEINSEN,"PRIMTARGET",&colnum,&status);
  fits_read_col(infits,TINT,colnum,spec->fiberID,1,1,&t,&(spec->tsf.primtarg),&anynul,&status);
  fits_report_error(stdout, status);
  if(status != 0)
    errormsg("Error occured reading Binary Table of file %s",spec->file);


  /* Convert Primtarg to targeting flags */  
  if (!SDSS_primtarg(&(spec->tsf)))
    errormsg("SDSS_rplate(): Error returned from SDSS_primtarg()");

  /* Formulate object name from object RA and DEC */
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
    errormsg("SDSS_rplate(): Error returned from ast_eq2gal()");
  spec->ra*=15.0;

  /* Move to Primary HDU */
  if(fits_movabs_hdu(infits,1,&hdutype,&status))
    errormsg("Could Not Move to HDU 1 in file %s",spec->file);
  if(hdutype != 0) 
    errormsg("HDU 1 not an Image in file %s",spec->file);

  /* Find number of spectral pixels */
  if (fits_read_key(infits,TINT,"NAXIS1",&(spec->np),comment,&status))
    errormsg("SDSS_rplate(): Cannot read value of header card %s\n\
\tfrom FITS file %s.","NAXIS1",spec->file);

  /* Allocate memory for wavelength, flux, error etc. arrays */
  if ((spec->wl=darray(spec->np))==NULL)   // Wavelength
    errormsg("SDSS_rplate(): Cannot allocate memory for\n\
\tspec->wl array of size %d",spec->np);
  if ((spec->fl=darray(spec->np))==NULL)   // Flux 
    errormsg("SDSS_rplate(): Cannot allocate memory for\n\
\tspec->fl array of size %d",spec->np);
  if ((spec->er=darray(spec->np))==NULL)   // Error
    errormsg("SDSS_rplate(): Cannot allocate memory for\n\
\tspec->er array of size %d",spec->np);

  /* Make the Median array, only if we are not using synflux */
  if(!set->synflux){
    if ((spec->co=darray(spec->np))==NULL)   // Median
      errormsg("SDSS_rplate(): Cannot allocate memory for\n\
\tspec->co array of size %d",spec->np);
  }

  if ((spec->pc=darray(spec->np))==NULL)  // Continuum
    errormsg("SDSS_rplate(): Cannot allocate memory for\n\
\tspec->pc array of size %d",spec->np);
  if ((spec->sn=darray(spec->np))==NULL)  // S/N 
    errormsg("SDSS_rplate(): Cannot allocate memory for\n\
\tspec->sn array of size %d",spec->np);
  if ((spec->st=iarray(spec->np))==NULL)  // Status
    errormsg("SDSS_rplate(): Cannot allocate memory for\n\
\tspec->st array of size %d",spec->np);

  /* Set default continuum level to Joe Liske's favourite value */
  for (i=0; i<spec->np; i++) { spec->co[i]=-10.0; spec->pc[i]=0.0; }


  /* Read in wavelength information */
  if (fits_read_key(infits,TDOUBLE,deltwl_key,&(spec->deltwl),comment,&status))
    errormsg("SDSS_rplate(): Cannot read value of header card %s\n\
\tfrom FITS file %s!",deltwl_key,spec->file);
  if (fits_read_key(infits,TDOUBLE,cpix_key,&(spec->cpix),comment,&status))
    errormsg("SDSS_rplate(): Cannot read value of header card %s\n\
\tfrom FITS file %s!",cpix_key,spec->file);
  if (fits_read_key(infits,TDOUBLE,beginwl_key,&(spec->beginwl),comment,&status))
    errormsg("SDSS_rplate(): Cannot read value of header card %s\n\
\tfrom FITS file %s!",beginwl_key,spec->file);
  if(fits_read_key(infits,TINT,"DC-FLAG",&loglin,comment,&status))
    errormsg("SDSS_rplate(): Cannot read value of header card %s\n\
\tfrom FITS file %s!","DC-FLAG",spec->file);

  /* Create wavelength scale from header cards */
  if (loglin == 1)
    /* Logarithmic wavelength scale */
    for (i=0; i<spec->np; i++)
      spec->wl[i]=pow(10.0,spec->beginwl+((double)(i+1)-spec->cpix)*spec->deltwl);
  else
    /* Linear wavelength scale */
    for (i=0; i<spec->np; i++)
      spec->wl[i]=spec->beginwl+((double)(i+1)-spec->cpix)*spec->deltwl;

  /* Record dispersion in km/s */
  spec->disp=(spec->wl[1]-spec->wl[0])/spec->wl[0]*C_C_K;

  /* Read in flux information */
  if(fits_read_pix(infits,TDOUBLE,fpixel,spec->np,NULL,spec->fl,NULL,&status)){
    fits_get_errstatus(status,buffer);
    errormsg("%s on flux of file %s",buffer,spec->file);
  }

  /* Move to HDU 2*/
  if(fits_movabs_hdu(infits,2,&hdutype,&status))
    errormsg("Could Not Move to HDU 2 in file %s",spec->file);
  if(hdutype != 0) 
    errormsg("HDU 2 not an Image in file %s",spec->file);

  /* Read in error information */
  if(fits_read_pix(infits,TDOUBLE,fpixel,spec->np,NULL,spec->er,NULL,&status)){
    fits_get_errstatus(status,buffer);
    errormsg("%s on error of file %s",buffer,spec->file);
  }

  /* Define signal-to-noise and status arrays */
  /* Since HDU 2 actually contains invvar, we convert to er */
  for (i=0; i<spec->np; i++) {
    spec->st[i]=(spec->er[i]>0.0) ? 1 : 0; 
    spec->er[i]=spec->st[i] ? pow(spec->er[i],-.5) : 0;   
    spec->sn[i]=(spec->st[i]) ? spec->fl[i]/spec->er[i] : 0.0;
  }

  /* Close input FITS file */
  fits_close_file(infits,&status);

  /* Because spPlate files dont contain Z information, we construct the name of where to look
     Which in this case is an spZbest file                                                   */

  cptr=getenv("SPECTRO_DATA");
  //sprintf(buffer,"%s/%04i/spZbest-%04i-%i.fits",cptr,spec->plate,spec->plate,spec->mjd);
  // TODO: Change directory to correct one for real data
  sprintf(buffer,"../Fits/spZbest-%04i-%i.fits",spec->plate,spec->mjd);

  /* Open FITS file */

  printf("spZBest Filename: %s\n", buffer);
  if(fits_open_file(&infits,buffer,READONLY,&status)){;
    fits_get_errstatus(status,buffer);
    errormsg("%s on spZbest file of %s",buffer,spec->file);
  }
  /* Move to HDU 2 */
  if(fits_movabs_hdu(infits,2,&hdutype,&status))
    errormsg("Could Not Move to HDU 2 in file %s",spec->file);
  if(hdutype != 2) 
    errormsg("HDU 2 not a Binary Table in file %s",strrchr(buffer,'/'));

  /* READ Z and Zerr Data */
  fits_get_colnum(infits,CASEINSEN,"Z",&colnum,&status);
  fits_read_col(infits,TDOUBLE,colnum,spec->fiberID,1,1,&s,&(spec->zem),&anynul,&status);
//  fits_get_colnum(infits,CASEINSEN,"Z_ERR",&colnum,&status);
//  fits_read_col(infits,TDOUBLE,colnum,spec->fiberID,1,1,&s,&(spec->zem),&anynul,&status);

  /* If Requested read in synflux */
  if(set->synflux){
      /* Move to HDU 3 */
    if(fits_movabs_hdu(infits,3,&hdutype,&status))
      errormsg("Could Not Move to HDU 3 in file %s",spec->file);
    if(hdutype != 0) 
      errormsg("HDU 3 not a Binary Table in file %s",strrchr(buffer,'/'));
      /* Read in flux information */
    if(fits_read_pix(infits,TDOUBLE,fpixel,spec->np,NULL,spec->pc,NULL,&status)){
      fits_get_errstatus(status,buffer);
      errormsg("%s on spZbest of %s",buffer,spec->file);
    }
  }
  fits_close_file(infits,&status);

  return 1;
}
