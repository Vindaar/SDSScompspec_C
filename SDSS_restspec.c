#include <stdio.h>
#include <stdlib.h>
#include "SDSS-spec.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int SDSS_restspec(spectrum *spec,restspec *rspec) {


  int n=0,i=0,bin=0;
  double wl=0,frac=0;
  double weight=0;
  double *temp=NULL;
  restspec tempspec;

  if((tempspec.wl=darray(rspec->np))==NULL)
    errormsg("Could not allocate Memory for temporary wavelength array");
  if((tempspec.fl=darray(rspec->np))==NULL)
    errormsg("Could not allocate Memory for temporary flux array");

  if((temp=darray(spec->np))==NULL)
    errormsg("Could not allocate Memory for temp array");


  for(i=0;i<spec->np;i++) {
    wl=spec->wl[i]/(spec->zem+1);
    frac+=(wl< 1430 && wl>1400 ? spec->fl[i] : 0);
    n+=(wl< 1430 && wl>1400 ? 1 : 0);
  }
  frac/=n;
  for(i=0;i<rspec->np;i++) {
    tempspec.fl[i]=frac*rspec->fl[i];
    tempspec.wl[i]=(1+spec->zem)*rspec->wl[i];
  }

  if(spec->wl[0] > tempspec.wl[0])
    while(tempspec.wl[bin]<=spec->wl[0]) bin++;
  bin--;

  for(i=0;i<spec->np;i++) {
    if(spec->wl[i+1]<tempspec.wl[0]) {
      spec->st[i]=0;
      continue;
    }
    while(tempspec.wl[bin] < spec->wl[i+1]) {
      weight=MIN(spec->wl[i+1],rspec->wl[bin+1])-MAX(spec->wl[i],rspec->wl[bin])/(rspec->wl[bin+1]-rspec->wl[bin]);

      spec->pc[i]+=tempspec.fl[bin]*weight;
      temp[i]+=weight;
      bin++;
    }
    spec->pc[i]/=temp[i];
  }

  free(temp); free(tempspec.fl); free(tempspec.wl);
  return 1;

}
