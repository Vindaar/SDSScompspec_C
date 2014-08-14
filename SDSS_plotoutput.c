#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SDSS-spec.h"
#include "SDSS_plotoutput.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int SDSS_plotoutput(compspec *cspec){

	return 0;
}
// Code UNFINISHED TODO: Check if important!
//
//  int pgset=0, pg_id=0;
//  int i=0,j=0;
//  float *x=NULL,*y=NULL, *z=NULL;
//  char title_det[NAMELEN]="\0";
//  char option[64]="\0";
//  char dummy1[64]="\0",dummy2[64]="\0";
//  stats stat[3];
//  dipole dip={.THETA=0,.BETA=0,.PHI=0,.shift=0,.split=0};
//  extern plotenv plenv;
//  /* PLOTTING FLAGS & OPTIONS */
//  _Bool split=0,clean=1,newstats=1,endit=0;
//  double xmin=2,xmax=5,ymin=0,ymax=1.4;
//
//  if(!pgset) {
//    cpgqid(&pg_id); SDSS_pgenv_init(); pg_open(&plenv,"?\0",progname,1);pgset++;
//  }
//
//  /* Allocate memory for data plotting arrays */
//  if ((x=farray(set->nzpart-1))==NULL)
//     errormsg("SDSS_plot_spec(): Cannot allocate memory\n\
//\tfor x plotting array of size %d",set->nzpart-1);
//  if ((y=farray(set->nzpart-1))==NULL)
//     errormsg("SDSS_plot_spec(): Cannot allocate memory\n\
//\tfor y plotting array of size %d",set->nzpart-1);
//  if ((z=farray(set->nzpart-1))==NULL)
//     errormsg("SDSS_plot_spec(): Cannot allocate memory\n\
//\tfor z plotting array of size %d",set->nzpart-1);
//  /* Clean the plotting surface */
//  cpgsvp(plenv.vpl,plenv.vpr,plenv.vpd,plenv.vpu);
//  plenv.nxsub=plenv.nysub=1; cpgsubp(plenv.nxsub,plenv.nysub); cpgpage();
//  /* Plotting Bounds & Labels*/
//      sprintf(title_det,"Comp Spec");
//      plenv.xmin[0]=xmin;plenv.xmax[0]=xmax;plenv.ymin[0]=ymin;plenv.ymax[0]=ymax;
//      sprintf(plenv.xlab[0],"%s","Redshift");
//      sprintf(plenv.ylab[0],"%s","Opacity");
//      sprintf(plenv.title[0],"%s   %s","Opacity vs Redshift",title_det);
//      /* Create custom made plotting surface */
//      cpgsci(0); cpgsls(1);
//      cpgswin(plenv.xmin[0],plenv.xmax[0],plenv.ymin[0],plenv.ymax[0]);
//      cpgsci(15);
//      cpgrect(plenv.xmin[0],plenv.xmax[0],plenv.ymin[0],plenv.ymax[0]);
//      cpgsci(plenv.cya); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
//      cpgsci(plenv.yel); cpgbox("N",0.0,0,"N",0.0,0);
//      cpgsci(plenv.gre); cpglab(plenv.xlab[0],plenv.ylab[0],"\0");
//      cpgsch(0.6*plenv.ch); cpglab("\0","\0",plenv.title[0]);
//      cpgsch(plenv.ch);
//      file:///usr/share/doc/HTML/index.html
//      clean=!clean;
//
//    /* Plot the data */
//    cpgsci(4); cpgslw(2*plenv.lw);
//    for(i=0; cspec->nhist[i] > 0;i++){x[i]=cspec->wl[i];y[i]=cspec->fl[i];}
//    cpgline(i,x,y);
//    cpgslw(1.5*plenv.lw);
//    for(i=0; cspec->nhist[i] > 0;i++){y[i]=cspec->er[i];
//    }
//    cpgline(i,x,y);
//
//  cpgclos();
//  /* Free memory space */
///*  free(x);free(y);free(u);free(v);free(yerrbn);free(yerrbx);free(verrbn);free(verrbx);*/
////  free(stat.mean);free(stat.nspectra);free(stat.stddev);free(stat.zbar);
//  return 1;
//}













