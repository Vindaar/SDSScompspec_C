#include <stdlib.h>
#include <stdio.h>
#include "SDSS-spec.h"
#include <string.h>
#include <math.h>
#include "memory.h"
#include "error.h"

int SDSS_restinit(restspec **rspec,setblock *set) {

  int i=0,nhdu=0,np=0;
  char *cptr=NULL;
  double beginwl=0,deltawl=0;

  /* Fits declarations */
  int      status=0;
  long     first[2]={1,1};
  fitsfile *infits=NULL;


  /* Parse the file name */
  cptr=strtok(set->restspec,"[,]");
  set->hdu = atoi(strtok(NULL,"[,]"));
  strcpy(set->restspec,cptr);

  fits_open_file(&infits,set->restspec,READONLY,&status);
  fitserrmsg(status);
  fits_read_key(infits,TINT,"NSPEC",&nhdu,NULL,&status);
  fitserrmsg(status);

  if(set->hdu > nhdu) 
    errormsg("HDU number exceeds Number of HDUs in file");

  fits_movabs_hdu(infits,set->hdu,NULL,&status);
  fitserrmsg(status);
  fits_read_key(infits,TINT,"NAXIS1",&np,NULL,&status);
  fits_read_key(infits,TDOUBLE,"CRVAL1",&beginwl,NULL,&status);
  fits_read_key(infits,TDOUBLE,"CD1_1",&deltawl,NULL,&status);
  fitserrmsg(status);

  (*rspec)=(restspec *)malloc((size_t)sizeof(restspec));

  (*rspec)->np=np;

  if(((*rspec)->wl=darray(np))==NULL)
    errormsg("Could not allocate for restspec array");
  if(((*rspec)->fl=darray(np))==NULL)
    errormsg("Could not allocate for restspec array");

  for(i=0;i<np;i++) 
    (*rspec)->wl[i]=pow(10,beginwl+deltawl*i);

  (*rspec)->minwl=(*rspec)->wl[0];(*rspec)->maxwl=(*rspec)->wl[np-1];

  fits_read_pix(infits,TDOUBLE,first,np,0,(*rspec)->fl,NULL,&status);
  fitserrmsg(status);


  fits_close_file(infits,&status);


  return 1;
}
