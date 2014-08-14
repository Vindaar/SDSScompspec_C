#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SDSS-restspec.h"
#include "stats.h"

int SDSS_reststat(compspec *cspec,int *nrun){

  int i=0,indent=0,sidx=0,sidx2=0,eidx=0,irun=0;
  double norm=0;
  statset stats;
  extern setblock set;


  for(irun=0;irun<(*nrun);irun++) {
    sidx=0;

    while(cspec[irun].n[sidx]<MINSIGHT) sidx++;

    for(i=0;i<cspec[irun].np-sidx;i++) {
      cspec[irun].wl[i]=cspec[irun].wl[i+sidx];
      cspec[irun].fl[i]=cspec[irun].fl[i+sidx];
      cspec[irun].er[i]=cspec[irun].er[i+sidx];
      cspec[irun].n[i]=cspec[irun].n[i+sidx];

      if(cspec[irun].n[i]>MINSIGHT) eidx=i;
    }

    cspec[irun].np=eidx+1;
    if((cspec[irun].wl=(double *)
	realloc(cspec[irun].wl,(size_t)cspec[irun].np*sizeof(double)))==NULL)
      errormsg("Could not reallocate memory for compspec wavelength array");
    if((cspec[irun].fl=(double *)
	realloc(cspec[irun].fl,(size_t)cspec[irun].np*sizeof(double)))==NULL)
      errormsg("Could not reallocate memory for compspec flux array");
    if((cspec[irun].er=(double *)
	realloc(cspec[irun].er,(size_t)cspec[irun].np*sizeof(double)))==NULL)
      errormsg("Could not reallocate memory for compspec error array");
    if((cspec[irun].n=(int *)
	realloc(cspec[irun].n,(size_t)cspec[irun].np*sizeof(int)))==NULL)
      errormsg("Could not reallocate memory for compspec histogram array");

    cspec[irun].meanz/=cspec[irun].nqso;
    cspec[irun].beginwl=.0001*rint(log10(cspec[irun].wl[0])/.0001);

    if(irun>0) {
      sidx=rint((MAX(cspec[irun].beginwl,log10(LYA))-cspec[0].beginwl)
		/cspec[0].deltwl);

      sidx2=rint((MAX(log10(LYA),cspec[0].beginwl)-cspec[irun].beginwl)
		 /cspec[0].deltwl);

      sidx=MAX(sidx,0);
      sidx2=MAX(sidx2,0);

      median(&(cspec[0].fl[sidx]),&(cspec[0].st[sidx]),
	     MIN((cspec[0].np-sidx),(cspec[irun].np-sidx2)),&stats,0);
      norm=stats.med;
      median(&(cspec[irun].fl[sidx2]),&(cspec[irun].st[sidx2]),
	     MIN((cspec[0].np-sidx),(cspec[irun].np-sidx2)),&stats,0);
      norm/=stats.med;
    }
    else norm=1;

    for(i=0;i<cspec[irun].np;i++) {
      /* Normalize and calculate error for Geometric Mean */
      if(set.gear) {
	cspec[irun].er[i]=cspec[irun].fl[i]
	  *exp(sqrt(cspec[irun].er[i]-pow(log(cspec[irun].fl[i]),2))
	       /sqrt(cspec[irun].n[i]-1));
	cspec[irun].fl[i]*=norm;
	cspec[irun].er[i]*=norm;
      }
      /* Or for the Arithmetic Mean */
      else {
	cspec[irun].er[i]=sqrt(cspec[irun].er[i]
			       -pow(cspec[irun].fl[i],2)
			       /(cspec[irun].n[i]-1));
	cspec[irun].fl[i]*=norm;
	cspec[irun].er[i]*=norm;
      }
    }
  }
  return 0;

}
