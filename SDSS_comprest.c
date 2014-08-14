/* 
   THis is the meat of the program, it shifts the wavelength array,
   rebins the data, and renormalizes.  The renormalization procedure,
   is to compare the median fluxes in the overlap of the new curve
   with the existing composite (excluding the forest regions).  Using
   the median requires sorting arrays, which is very time
   consuming. This function does not alter the data.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SDSS-restspec.h"
#include "stats.h"

int SDSS_comprest(spectrum *specarray,compspec *cspec,int nrun,int idx) {

  int i=0,np=0,msidx=0,meidx=0;
  int indent1=0,indent2=0;
  int irun=0;
  double norm=0,weight=0;
  double beginwl=0,c1=.0001;
  double *fl=NULL,*er=NULL,*wl=NULL;
  double dtemp1=0,dtemp2=0;
  spectrum *spec;
  statset tempstats;
  extern setblock set;

  spec=&(specarray[idx]);
  np=spec->np;
  wl=(double *)malloc((size_t)np*sizeof(double));
  fl=(double *)malloc((size_t)np*sizeof(double));
  er=(double *)malloc((size_t)np*sizeof(double));

  /* PART I: Shift to rest-frame */
  beginwl=c1*rint((spec->beginwl-log10(1+spec->zem))/c1);
  for(i=0;i<np;i++)
    wl[i]=pow(10,beginwl+c1*i);
  
  /* PART II: Normalize */
  /* Select the Overlap region in the spectrum */
  msidx=(int)(rint(log10(MAX(MAX(LYA,wl[0]),cspec->wl[0]))/c1)-beginwl/c1);
  meidx=(int)(rint(log10(MIN(wl[np-1],cspec->wl[cspec->np-1]))/c1)-beginwl/c1);

  /* Find the Median */
  median(&(spec->fl[msidx]),&(spec->st[msidx]),meidx-msidx+1,&tempstats,0);
  norm=tempstats.med;
  /* Select the overlap region in the BIG composite */
  msidx=(int)(rint(log10(MAX(MAX(LYA,wl[0]),cspec->wl[0]))/c1)
	      -cspec->beginwl/c1);
  meidx=(int)(rint(log10(MIN(wl[np-1],cspec->wl[cspec->np-1]))/c1)
	      -cspec->beginwl/c1);
  /* Find the Median */
  median(&(cspec[0].fl[msidx]),&(cspec[0].st[msidx]),
	 meidx-msidx+1,&tempstats,0);
  norm/=(tempstats.med==0 ? 1 : tempstats.med);

  /* PART III: Combine the spectra */
  /* Pick which composite we want to add the data to */
  irun=0;
  for(irun=1;irun<nrun;irun++) {
    if(idx>cspec[irun].sidx && idx<cspec[irun].eidx)
      break;
    else ;
  }
  if(irun == nrun) irun=0;
  /* Shift our focus on the BIG composite to the region of interest */
  indent1=rint((cspec[irun].beginwl-cspec[0].beginwl)/cspec[0].deltwl);
  cspec[0].fl+=indent1;cspec[0].s+=indent1;cspec[0].n+=indent1;
  cspec[0].beginwl+=spec->deltwl*indent1;
  cspec[0].st+=indent1;cspec[0].wl+=indent1;cspec[0].np-=indent1;

  /* Align the Composites with the new spectrum */
  indent2=rint((beginwl-cspec[0].beginwl)/c1);
  cspec[0].fl+=indent2;cspec[0].s+=indent2;cspec[0].n+=indent2;
  cspec[0].st+=indent2;cspec[0].wl+=indent2;cspec[0].np-=indent2;
  cspec[0].beginwl+=c1*indent2;
  if(irun>0) {
    cspec[irun].fl+=indent2;cspec[irun].s+=indent2;cspec[irun].n+=indent2;
    cspec[irun].st+=indent2;cspec[irun].wl+=indent2;cspec[irun].np-=indent2;
    cspec[irun].beginwl+=c1*indent2;
  }

  /* Normalize Spectrum and rebin*/
  for(i=0;i<np;i++) {
    fl[i]=spec->fl[i]/norm;
    er[i]=spec->er[i]/norm;
    spec->fl[i]/=norm;
    spec->er[i]/=norm;
    weight=1;

    /* If we are below the bottom of the array, or have a bad pixel then
       continue*/
    if((spec->st[i]==0) || (fl[i]<=0 && set.gear))
      continue;
    
    /* Find value to add into mean */
    if(set.gear) 
      dtemp1=log(fl[i]);
    else
      dtemp1=fl[i];

    /* Calculate Geometric Mean */
    if(set.gear) {
      /* For the Big Composite */
      dtemp2=(cspec[0].s[i]>0 ? 
	      log(cspec[0].fl[i])*cspec[0].s[i] : 0);
      cspec[0].fl[i]=
	exp((dtemp2 + weight*dtemp1)/(cspec[0].s[i]+weight));
      if(irun>0) {
	/* For the ith composite */
	dtemp2=(cspec[irun].s[i]>0 ? 
		log(cspec[irun].fl[i])*cspec[irun].s[i] : 0);
	cspec[irun].fl[i]=
	  exp((dtemp2 + weight*dtemp1)/(cspec[irun].s[i]+weight));
      }
    }
    /* Or the Arithmetic Mean, if you like */
    else { 
      /* For the Big Composite */
      dtemp1=(cspec[0].s[i] ? 
	      cspec[0].fl[i]*cspec[0].s[i] : 0);
      cspec[0].fl[i]=
	(dtemp2 + weight*dtemp1)/(cspec[0].s[i]+weight);
      if(irun > 0) {
	/* For the ith composite */
	dtemp1=(cspec[irun].s[i] ? 
		cspec[irun].fl[i]*cspec[irun].s[i] : 0);
	cspec[irun].fl[i]=
	  (dtemp2 + weight*dtemp1)/(cspec[irun].s[i]+weight);
      }
    }
    cspec[0].n[i]+=1;
    cspec[0].s[i]+=weight;
    cspec[0].st[i]=(cspec[0].n[i] > 0 ? 1 : 0);
    if(irun>0) {
      cspec[irun].n[i]+=1;
      cspec[irun].s[i]+=weight;
      cspec[irun].st[i]=(cspec[irun].n[i] > 0 ? 1 : 0);
    }
  }

  cspec[0].meanz+=spec->zem;
  cspec[0].nqso++;
  if(irun>0) {
    cspec[irun].meanz+=spec->zem;
    cspec[irun].nqso++;
  }

  cspec[0].fl-=indent2;cspec[0].s-=indent2;cspec[0].n-=indent2;
  cspec[0].st-=indent2;cspec[0].wl-=indent2;cspec[0].np+=indent2;
  cspec[0].beginwl-=c1*indent2;
  if(irun>0) {
    cspec[irun].fl-=indent2;cspec[irun].s-=indent2;cspec[irun].n-=indent2;
    cspec[irun].st-=indent2;cspec[irun].wl-=indent2;cspec[irun].np+=indent2;
    cspec[irun].beginwl-=c1*indent2;
  }
  cspec[0].fl-=indent1;cspec[0].s-=indent1;cspec[0].n-=indent1;
  cspec[0].st-=indent1;cspec[0].np+=indent1;cspec[0].wl-=indent1;
  cspec[0].beginwl-=spec->deltwl*indent1;
  free(wl);free(fl);free(er);
  return 0;
}
