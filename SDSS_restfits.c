#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SDSS-restspec.h"

void fitserrormsg(int fitstatus) {

  char comment[32]="\0";

  if(!fitstatus) return;
  else  {
    fits_get_errstatus(fitstatus,comment);
    errormsg("FITS ERROR: %s",comment);
  }
}

int SDSS_restfits(compspec *cspec,int nrun) {

  int i=0;
  time_t tm;
  double dtmp=0.0;
  char buffer[64]="\0",meantype[2][64]={"Arithmetic","Geometric"};
  extern char *progname;
  extern setblock set;

  /* FITS WRITE DECLARATIONS */
  fitsfile *outfits=NULL;
  char filename[129]="\0";
  int fitstatus=0;
  int naxis=2;
  int crpix=1;
  long naxes[2]={0,2};
  long fpixel[3]={1,3};

  tm=time(NULL);

  if(set.force) strcpy(filename,"!");
  strcat(filename,set.outfile);

  if(fits_create_file(&outfits,filename,&fitstatus))
    fitserrormsg(fitstatus);

  if(fits_create_img(outfits,-32,naxis,naxes,&fitstatus))
    fitserrormsg(fitstatus);
  sprintf(buffer,"%s",strtok(asctime(gmtime(&tm)),"\n"));
  if(fits_write_key(outfits,TSTRING,"AUTHOR",AUTHOR,"File Creator",&fitstatus) ||
     fits_write_key(outfits,TSTRING,"DATE",buffer,"File Created On (UTC)",&fitstatus) ||
     fits_write_key(outfits,TSTRING,"PROGRAM",progname,"Created by program",&fitstatus) ||
     fits_write_key(outfits,TINT,"NSPEC",&nrun,"Number of Composite Spectra",&fitstatus))
    fitserrormsg(fitstatus);

  for(i=0;i<nrun;i++) {

    naxes[0]=cspec->np;
    if(fits_create_img(outfits,-32,naxis,naxes,&fitstatus))
      fitserrormsg(fitstatus);

    fpixel[1]=1;
    sprintf(buffer,"[%4.3f,%4.3f]",cspec[i].minz,cspec[i].maxz);

    dtmp=cspec[i].maxz-cspec[i].minz;
    fits_write_key(outfits,TDOUBLE,"MEANZ",&(cspec[i].meanz),
		   "Mean Redshift used in this spectrum",&fitstatus);
    fitserrormsg(fitstatus);
    fits_write_key(outfits,TDOUBLE,"MINZ",&(cspec[i].minz),
		   "Min Redshift used in this spectrum",&fitstatus);
    fitserrormsg(fitstatus);
    fits_write_key(outfits,TDOUBLE,"MAXZ",&(cspec[i].maxz),
		   "Max Redshift used in this spectrum",&fitstatus);
    fitserrormsg(fitstatus);
    fits_write_key(outfits,TSTRING,"ZINTRVL",buffer,
		   "Min/Max Redshift",&fitstatus);
    fitserrormsg(fitstatus);
    fits_write_key(outfits,TDOUBLE,"DZSIZE",&dtmp,
		   "Size of redshift interval",&fitstatus);
    fitserrormsg(fitstatus);
    fits_write_key(outfits,TINT,"NQSOs",&(cspec[i].nqso),
		   "Number of QSOs contributing to this spectrum",&fitstatus);
    fitserrormsg(fitstatus);
    fits_write_key(outfits,TDOUBLE,"CRVAL1",&(cspec[i].beginwl),
		   "log10 of first pixel",&fitstatus);
    fitserrormsg(fitstatus);
    fits_write_key(outfits,TDOUBLE,"CD1_1",&(cspec[i].deltwl),
		   "log10 of dispersion",&fitstatus);
    fitserrormsg(fitstatus);
    fits_write_key(outfits,TINT,"CRPIX1",&crpix,"Indexing Value",&fitstatus);
    fitserrormsg(fitstatus);    
    fits_write_key(outfits,TSTRING,"MEANTYP",meantype[set.gear],
		   "Type of Mean Used",&fitstatus);
    
    fitserrormsg(fitstatus);


    if(fits_write_pix(outfits,TDOUBLE,fpixel,cspec[i].np,cspec[i].fl,&fitstatus))
      fitserrormsg(fitstatus);
    fpixel[1]=2;
    if(fits_write_pix(outfits,TDOUBLE,fpixel,cspec[i].np,cspec[i].er,&fitstatus))
      fitserrormsg(fitstatus);
  }
  fits_close_file(outfits,&fitstatus);

  return 0;
}
