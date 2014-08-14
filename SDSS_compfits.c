#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "SDSS-spec.h"
#include "astron.h"
#include "memory.h"
#include "error.h"

int SDSS_compfits(compspec *cspec, spectrum *spec, int nspec,int speccount,
		  setblock *set){

  fitsfile *outfits=NULL;
  int i=0,j=0,status=0;
  int sidx=0;
  int *p=NULL;
  double *s=NULL,*t=NULL,**u=NULL;
  double Restrange[][2]=RESTRANGE;
  extern char *progname;
  char buffer[256]="!";
  time_t tm;


  /* FITS Header Info */
  int naxis=2,DR=5;
  long npix=0, naxes[2]={0,3};
  long fpixel[2]={1,1};
  int  crpix1=1,dc_flag=1,tru=1,fals=0; 
  double coeff1=1.00000000000000E-04,coeff0=0;
  char *ttype[]={"MJD","PLATE","FIB",
		 "RA","DEC","ALPHA","EALPHA","l","b","ZEM","MAGSPEC","FLAG"},
    *tform[]={"J","J","J","D","D","D","D","D","D","D","5D","I"},
      *tunit[]={"","","","","","","","","","","",""},
	extname[]="DataSet";

  tm=time(NULL);

  strcat(buffer,set->outfile);
  fits_create_file(&outfits,buffer,&status);

  for(j=0;j<set->nforest;j++) {
    i=0;
    while(cspec[j].nhist[++i] == 0) ;
    sidx=i;
    if(sidx>=2999) continue;

    for(i=sidx;i<3000;i++) 
      if(cspec[j].nhist[i] > 0)
	npix=i+1;
    naxes[0]=npix;

    coeff0=log10(cspec[j].wl[sidx]);

    if((s=darray(npix))==NULL)
      errormsg("Could not allocate memory");

    fits_create_img(outfits,-32,naxis,naxes,&status);
    fits_write_key(outfits,TSTRING,"PROGRAM",progname,"Generation program",&status);
    fits_write_key(outfits,TSTRING,"AUTHOR",AUTHOR,"Created by",&status);
    sprintf(buffer,"%s",strtok(asctime(gmtime(&tm)),"\n"));
    fits_write_key(outfits,TSTRING,"DATE",buffer,"DATE CREATED (UTC)",&status);
    fits_write_key(outfits,TSTRING,"ARRAY0","DA","1-Fo/Fc",&status);
    fits_write_key(outfits,TSTRING,"ARRAY1","ERR","1 SIGMA",&status);
    fits_write_key(outfits,TSTRING,"ARRAY2","NPOINTS","NUMBER OF CONTRIBUTING PIXELS",&status);
    fits_write_key(outfits,TINT,"NSPECTRA",&speccount,"Total number of contributing Quasars",&status);
    fits_write_key(outfits,TINT,"NSPEC",&(set->nforest),"Number of spectra (extensions) in this file",&status);
    fits_write_key(outfits,TLOGICAL,"VACUUM",&tru,"Wavelengths are in vacuum",&status);
    fits_write_key(outfits,TINT,"DC_FLAG",&dc_flag,"Log-Linear Flag",&status);
    fits_write_key(outfits,TINT,"CRPIX1",&crpix1,"Starting pixel (1-indexed)",&status);
    fits_write_key(outfits,TDOUBLE,"COEFF0",&coeff0,"Center wavelength (log10) of first pi",&status);
    fits_write_key(outfits,TDOUBLE,"COEFF1",&coeff1,"Log10 dispersion per pixel",&status);
    fits_write_key(outfits,TDOUBLE,"CRVAL1",&coeff0,"Iraf zero point",&status);
    fits_write_key(outfits,TDOUBLE,"CD1_1",&coeff1,"Iraf dispersion",&status);
    fits_write_key(outfits,TINT,"DR",&DR,"Data release data used",&status);
    sprintf(buffer,"[%5.1f,%5.1f]",Restrange[j][0],Restrange[j][1]);
    fits_write_key(outfits,TSTRING,"FOREST",buffer,"RF Wavelength Measured",&status);
    fits_write_key(outfits,TLOGICAL,"METAL",&fals,"Not Corrected for Metals",&status);
    fits_write_key_str(outfits,"COMMENT","Webb JK, Ouellet JL, Murphy M (2008)","",&status);

    fpixel[1]=0;

    /* apply metal correction if given as command line argument. Based on
    # \tau_corr(z) = 0.72(1+z)^(0.17) \tau_uncorr(z)
    # Schaye et al. (2003)
    # we need to assign a redshift to each pixel, thus we create the
    # zem array. Then, we first calculate the real optical depth from
    # hdu0_row1 (which technically is 1 - average transmittance) by
    # tau = - log10(average transmittance)
    # apply the correction and revert the optical depth back to
    # 1 - average transmittance
     */
    if (set->metal_corr == 1){
    	// create an array of redshifts:
    	//double tau_corr[npix];
    	double zem_temp = 0;
    	int q = 0;
    	double temp = 0;
    	for(q = 0; q < npix; q++){
    		// 1215.67: Ly alpha emission in rest frame
    		zem_temp = pow(10, coeff0 + coeff1*q)/1215.67 - 1;
    		printf("%f\t%f\t%i\t%f\t%i\n\n", coeff0, coeff1, npix, zem_temp, sidx);
    		temp = cspec[j].fl[sidx+q];
    		temp = -log10(temp);
            temp = 0.72*pow(1+zem_temp,0.17)*temp;
            cspec[j].fl[sidx+q] = pow(10, -temp);
    	}
    }

    for(i=0;i<npix;i++) s[i]=1-cspec[j].fl[sidx+i];
    fpixel[1]++;
    fits_write_pix(outfits,TDOUBLE,fpixel,npix,s,&status);

    for(i=0;i<npix;i++) s[i]=cspec[j].er[sidx+i];
    fpixel[1]++;
    fits_write_pix(outfits,TDOUBLE,fpixel,npix,s,&status);

    for(i=0;i<npix;i++) s[i]=1.0*cspec[j].nhist[sidx+i];
    fpixel[1]++;
    fits_write_pix(outfits,TDOUBLE,fpixel,npix,s,&status);
    free(s);
  }
  /* Quasar Statistics */
  if((t=darray(nspec))==NULL)
    errormsg("Could not allocate memory");

  fits_create_tbl(outfits,BINARY_TBL,0,12,ttype,tform,tunit,extname,&status);
  fits_write_key(outfits,TINT,"NSPEC",&nspec,"Total number of Quasars",&status);
  fits_write_key(outfits,TDOUBLE,"MEAN_A",&(cspec[0].mean_a),"Mean Spectral Index",&status);
  fits_write_key(outfits,TDOUBLE,"SIGMA_A",&(cspec[0].sigma_a),"STDDEV on alpha",&status);
  fits_write_key(outfits,TDOUBLE,"MED_A",&(cspec[0].med_a),"Median Spectral Index",&status);
  fits_write_key(outfits,TDOUBLE,"SIQR_A",&(cspec[0].siqr_a),"68% semi-interquartile range on alpha",&status);

  
  for(i=0;i<nspec;i++) t[i]=spec[i].ra;
  fits_write_col(outfits,TDOUBLE,4,1,1,nspec,t,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++) t[i]=spec[i].dec;
  fits_write_col(outfits,TDOUBLE,5,1,1,nspec,t,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++)  t[i]=spec[i].alpha; 
  fits_write_col(outfits,TDOUBLE,6,1,1,nspec,t,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++) t[i]=spec[i].ealpha;
  fits_write_col(outfits,TDOUBLE,7,1,1,nspec,t,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++) t[i]=spec[i].l;
  fits_write_col(outfits,TDOUBLE,8,1,1,nspec,t,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++) t[i]=spec[i].b;
  fits_write_col(outfits,TDOUBLE,9,1,1,nspec,t,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++) t[i]=spec[i].zem;
  fits_write_col(outfits,TDOUBLE,10,1,1,nspec,t,&status);
  fitserrmsg(status);
  free(t);
  if((p=iarray(nspec))==NULL)
    errormsg("Could not allocate memory for int array size %d\n",nspec);

  for(i=0;i<nspec;i++) p[i]=spec[i].mjd;
  fits_write_col(outfits,TINT,1,1,1,nspec,p,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++) p[i]=spec[i].plate;
  fits_write_col(outfits,TINT,2,1,1,nspec,p,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++) p[i]=spec[i].fiberID;
  fits_write_col(outfits,TINT,3,1,1,nspec,p,&status);
  fitserrmsg(status);
  for(i=0;i<nspec;i++) { 
    p[i]=(spec[i].usespec == 1 ? 0 : (spec[i].alpha<-900 ? 1 : 2 ));
  }
  fits_write_col(outfits,TINT,12,1,1,nspec,p,&status);
  fitserrmsg(status);

  u=dmatrix(nspec,5);
  for(i=0;i<nspec;i++) { 
    for(j=0;j<5;j++) 
      u[i][j]=spec[i].smag[j];
      printf("%s %lf %lf %lf %lf\n",spec[i].abfile,spec[i].zem,spec[i].smag[0],spec[i].smag[1],spec[i].smag[2]);
  }

  fits_write_col(outfits,TDOUBLE,11,1,1,nspec,u,&status);
  fitserrmsg(status);

  fits_close_file(outfits,&status);

  free(*u);
  free(u);
  free(p);
  return 1;
}


