/****************************************************************************
* Plot out an SDSS spectrum using PGPLOT
****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "SDSS-spec.h"
#include "memory.h"
#include "error.h"
extern plotenv plenv;

int SDSS_plot_spec(spectrum *spec) {

  float          *x=NULL,*y=NULL;
  //float          xmin,xmax;
  float 		 ymin,ymax;
  float          ftmp,temp=0;
  float          pgx=0,pgy=0,pgxo=0,pgyo=0;
  float          vpd0,vpu0,vpd1,vpu1;
  //float 		 vpd2,vpu2;
  int            sidx=0,eidx=0;
  //int            sp=0,ep=0;
  //int 			 cnp=0;
  int            nplot=0;
  int            i=0;
  //int 			 j=0,k=0;
  //int     		 flag=0;
  //char           c='y';
  char			 pgch[8]="p";
  // char 			 option[64]="\0";
  // char  		 dummy2[64]="\0";
  //float          zlya=spec->zem;

  //int            pspec = set.pspec;
  //int 	         medcont = set.medcont;
  

  // print some information about the spectrum
  printf("alpha: %f\t zem: %f\n", spec->alpha, spec->zem);

  /* Find starting and ending indices for continuum */
  ymax=1.01*djmax(spec->fl,spec->np);
  if(isinf(ymax)) ymax=100;
  ymin=-0.1*ymax;
  sidx=0; eidx=spec->np-1; nplot=eidx-sidx+1;

  vpd0=.08;vpu0=.18;vpd1=.22;vpu1=.95;

  plenv.xmin[0]=(plenv.xmin[1]=(plenv.xmin[2]=spec->wl[sidx]));
  plenv.xmax[0]=(plenv.xmax[1]=(plenv.xmax[2]=spec->wl[eidx]));
  plenv.ymin[0]=(plenv.ymin[1]=ymin);
  plenv.ymax[0]=(plenv.ymax[1]=ymax);
  strcpy(plenv.xlab[0],"Vacuum Wavelength \\A");
  strcpy(plenv.ylab[1],"10\\u-17\\d ergs/cm\\u2\\d/\\A");
  sprintf(plenv.title[0],"Quasar Spectra %s, z=%4.3f",spec->obj,spec->zem);

  while(1) {
    cpgbbuf();
    cpgpage();
    /* Detect any changes in size of plotting window */
    cpgqvsz(1,&ftmp,&(plenv.wwidth),&ftmp,&(plenv.wasp)); plenv.wasp/=plenv.wwidth;
    
    /* Allocate memory for data plotting arrays */
    if ((x=farray(nplot))==NULL)
       errormsg("SDSS_plot_spec(): Cannot allocate memory\n\
    \tfor x plotting array of size %d",nplot);
    if ((y=farray(nplot))==NULL)
       errormsg("SDSS_plot_spec(): Cannot allocate memory\n\
    \tfor y plotting array of size %d",nplot);


    for(i=0;i<nplot;i++){
      x[i]=spec->wl[i];y[i]=spec->fl[i];
    }

    /* Plot Spectrum Index */  
    cpgsvp(plenv.vpl,plenv.vpr,vpd0,vpu0);cpgsci(0);
    cpgswin(plenv.xmin[0],plenv.xmax[0],plenv.ymin[0],plenv.ymax[0]);
    cpgsci(15);
    cpgrect(plenv.xmin[0],plenv.xmax[0],plenv.ymin[0],plenv.ymax[0]);
    cpgsci(5);
    cpgsch(0.3*plenv.ch); cpgbox("BCTS",0.0,0,"BC",0.0,0); cpgsch(plenv.ch);
    cpgsci(7); cpgsch(0.8*plenv.ch); cpgbox("N",0.0,0,"",0.0,0);
    cpgsci(3); cpgmtxt("B",2.5,0.5,0.5,plenv.xlab[0]);
    cpgsci(1); cpgsch(plenv.ch);cpgbin(nplot,x,y,1);
    cpgsci(1);cpgbin(nplot,x,y,1);
    cpgsci(5); cpgsls(4);cpgmove(plenv.xmin[0],0.0); cpgdraw(plenv.xmax[0],0.0);
    cpgsfs(2);cpgsci(3);cpgsls(1);
    if(plenv.xmin[0] != plenv.xmin[1] || plenv.xmax[0] != plenv.xmax[1] ||  
      plenv.ymin[0] != plenv.ymin[1] || plenv.ymax[0] != plenv.ymax[1])
      cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[1],plenv.ymax[1]);
    cpgsci(0);
    cpgslw(plenv.lw); cpgsfs(1);

     /* Main Window */
    cpgsvp(plenv.vpl,plenv.vpr,vpd1,vpu1);cpgsci(0);
    cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[1],plenv.ymax[1]);
    cpgsci(15);
    cpgrect(plenv.xmin[1],plenv.xmax[1],plenv.ymin[1],plenv.ymax[1]);
    cpgsci(5);
    cpgsch(0.3*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0); cpgsch(plenv.ch);
    cpgsci(7); cpgsch(0.8*plenv.ch); cpgbox("N",0.0,0,"N",0.0,0);
    cpgsci(3); cpgmtxt("L",2.5,0.5,0.5,plenv.ylab[1]);
    cpgsci(1); cpgsch(plenv.ch); cpgbin(nplot,x,y,1);
    for(i=0;i<nplot;i++) y[i]=spec->er[i];
    cpgsci(2);cpgbin(nplot,x,y,1);
    for(i=0;i<nplot;i++) y[i]=spec->co[i];
    cpgsci(3);cpgbin(nplot,x,y,1);
    for(i=0;i<nplot;i++) y[i]=spec->pc[i];
    cpgsci(12);cpgslw(1.5*plenv.lw);cpgbin(nplot,x,y,1);
    cpgsci(5); cpgsls(4);cpgmove(plenv.xmin[1],0.0); cpgdraw(plenv.xmax[1],0.0);
    cpgsls(1); cpgsci(0);
    cpgslw(plenv.lw); cpgsfs(1);



    /* Labels */
    cpgsci(3);
    cpgmtxt("T",.5,.5,.5,plenv.title[0]);
    cpgsci(0);

    cpgebuf();
    /* Clean up */
    free(x); free(y);

    cpgsvp(plenv.vpl,plenv.vpr,0,1);
    cpgswin(plenv.vpl,plenv.vpr,0,1);
    cpgband(0,0,0.0,0.0,&pgx,&pgy,pgch);

    if(!strncmp(pgch,"n",1)) return 1;
    else if(!strncmp(pgch,"y",1)) return 2;
    else if(!strncmp(pgch,"b",1)) return -1;
    else if(!strncmp(pgch,"q",1)) return 10;
    else if(!strncmp(pgch,"A",1) && pgy<vpu0 && pgy > vpd0 && pgx>plenv.vpl && pgx<plenv.vpr){
      /* Use low-dispersion spectrum to navigate combined spectrum */
      cpgsvp(plenv.vpl,plenv.vpr,vpd0,vpu0);
      cpgswin(plenv.xmin[0],plenv.xmax[0],plenv.ymin[0],plenv.ymax[0]);
      cpgsci(7); cpgslw(2.0*plenv.lw);
      cpgsch(0.3*plenv.ch); cpgbox("BCTS",0.0,0,"BC",0.0,0); cpgsch(plenv.ch);
      cpgslw(plenv.lw); cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
      if (!strncmp(pgch,"?",1)) {
        fprintf(stderr,"\
instructions for navigation using 'thumbnail' spectrum:\n\
 Left mouse  : Re-centre detailed spectrum (flux panel) at selected\n\
                 wavelength\n\
 Middle mouse: Re-centre detailed spectrum (as above) and re-size\n\
                 window to suit new local extrema\n\
 Right mouse : Select custom window region.\n\
                 Right mouse again selects second corner of new window,\n\
                 any other entry to abort selection.\n\n");
        cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
      }
      else if (!strncmp(pgch,"X",1)) {
        strcpy(pgch,"\0"); cpgsci(2); pgxo=pgx; pgyo=pgy;
        cpgband(2,0,pgxo,pgyo,&pgx,&pgy,pgch);
        if (!strncmp(pgch,"X",1)) {
          if (pgxo>pgx) { temp=pgx;pgx=pgxo;pgxo=temp; }
          if (pgyo>pgy) { temp=pgy;pgy=pgyo;pgyo=temp; }
          pgxo=(MAX(pgxo,spec->wl[0])); pgx=(MIN(pgx,spec->wl[spec->np-1]));
          plenv.xmin[1]=pgxo;
          plenv.xmax[1]=pgx;
          plenv.ymax[1]=MIN(pgy,ymax);
          plenv.ymin[1]=MAX(pgyo,ymin);
        }
      }
      else if (!strncmp(pgch,"A",1)) {
        pgx=(MIN(spec->wl[spec->np-1],(MAX(pgx,spec->wl[0]))));
        temp=(plenv.xmax[1]-plenv.xmin[1])/2;
        plenv.xmin[1]=pgx-temp; plenv.xmax[1]=pgx+temp;
      }
    }
    else if (!strncmp(pgch,"A",1) && (pgy>vpd1 && pgy<vpu1)) {
      cpgsvp(plenv.vpl,plenv.vpr,vpd1,vpu1);
      cpgswin(plenv.xmin[1],plenv.xmax[1],plenv.ymin[1],plenv.ymax[1]);
      cpgsci(7); cpgslw(2.0*plenv.lw);
      cpgsch(0.9*plenv.ch); cpgbox("BCTS",0.0,0,"BCTS",0.0,0);
      cpgsch(plenv.ch); cpgslw(plenv.lw); cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
      if (!strncmp(pgch,"?",1)) {
        fprintf(stderr,"\
Options for navigation/editing using detailed spectrum (flux panel):\n\
 Left mouse  : Re-centre panel at selected wavelength\n\
 Middle mouse: Re-centre panel (as above) and re-size window to suit\n\
                 new local extrema\n\
 Right mouse : Select custom window region.\n\
                 Right mouse again selects second corner of new window,\n\
                 any other entry to abort selection.\n\
 Space bar   : Write out to terminal information about cursor position.\n\
 b           : Define new bottom limit for plot\n\
 c           : Clip pixels from either combined or contributing spectra (act).\n\
                 c again selects second corner of clip window,\n\
               any other entry to abort selection.\n\
 l           : Define new left limit for plot\n\
 m           : Calculate normalized statistics over region (*)\n\
 r           : Define new right limit for plot\n\
 t           : Define new top limit for plot\n\
 u           : Unclip pixels from either combined or contributing spectra (act).\n\
                 u again selects second corner of clip window,\n\
                 any other entry to abort selection.\n\
 w           : Find equivalent width between this mark and next (*)\n\n");
        cpgband(7,0,0.0,0.0,&pgx,&pgy,pgch);
      }
      if (!strncmp(pgch,"D",1)) {
        pgx=(MIN(spec->wl[spec->np-1],(MAX(pgx,spec->wl[0]))));
        temp=.5*(plenv.xmax[1]-plenv.xmin[1]);
        plenv.xmin[1]=MAX(spec->wl[0],pgx-temp); plenv.xmax[1]=MIN(pgx+temp,spec->wl[spec->np-1]);
        plenv.ymax[1]=0;
        for(i=0;i<spec->np;i++) { 
          if(spec->wl[i]>plenv.xmin[1] && spec->wl[i]<plenv.xmax[1]) {
            plenv.ymax[1]=MAX(plenv.ymax[1],spec->fl[i]);
          }
        }
        plenv.ymin[1]=-0.1*plenv.ymax[1];
      }
      else if (!strncmp(pgch,"X",1)) {
        strcpy(pgch,"\0"); cpgsci(2); pgxo=pgx; pgyo=pgy;
        cpgband(2,0,pgx,pgy,&pgx,&pgy,pgch);
        if (!strncmp(pgch,"X",1)) {
          if (pgxo>pgx) { temp=pgx;pgx=pgxo;pgxo=temp; }
          if (pgyo>pgy) { temp=pgy;pgy=pgyo;pgyo=temp; }
          pgxo=(MAX(pgxo,spec->wl[0])); pgx=(MIN(pgx,spec->wl[spec->np-1]));
          plenv.xmin[1]=pgxo; plenv.xmax[1]=pgx;
          plenv.ymax[1]=MIN(ymax,pgy);
          plenv.ymin[1]=MAX(-0.1*plenv.ymax[1],pgyo);
        }
      }
      else if (!strncmp(pgch,"A",1)) {
        pgx=(MIN(spec->wl[spec->np-1],(MAX(pgx,spec->wl[0]))));
        temp=.5*(plenv.xmax[1]-plenv.xmin[1]);
        plenv.xmin[1]=MAX(spec->wl[0],pgx-temp); plenv.xmax[1]=MIN(pgx+temp,spec->wl[spec->np-1]);
        temp=.5*(plenv.ymax[1]-plenv.ymin[1]);
        plenv.ymax[1]=MIN(ymax,pgy+temp);
        plenv.ymin[1]=MAX(ymin,pgy-temp);
      }
      else if (!strncmp(pgch,"l",1)) {
        plenv.xmin[1]=MAX(plenv.xmin[1],pgx);

      }
      else if (!strncmp(pgch,"r",1)) {
        plenv.xmax[1]=MIN(pgx,plenv.xmax[1]);
      }
      else if (!strncmp(pgch,"t",1)) {
        plenv.ymax[1]=(pgy>plenv.ymin[1]) ? pgy : plenv.ymax[1];
      }
      else if (!strncmp(pgch,"b",1)) {
        plenv.ymin[1]=(pgy<plenv.ymax[1]) ? pgy : plenv.ymin[1];
      }
      else if (!strncmp(pgch," ",1)) {
	fprintf(stderr,"(%5.1f,%5.3f)\n",pgx,pgy);
      }
    }
    else if (!strncmp(pgch,".",1)) {
      temp=0.5*(plenv.ymax[1]-plenv.ymin[1]);
      plenv.ymax[1]+=temp; plenv.ymin[1]-=temp;
    }
    else if (!strncmp(pgch,",",1)) {
      temp=0.25*(plenv.ymax[1]-plenv.ymin[1]);
      plenv.ymax[1]-=temp; plenv.ymin[1]+=temp;
    }
    else if (!strncmp(pgch,">",1)) {
      temp=.5*(plenv.xmax[1]-plenv.xmin[1]);
      plenv.xmax[1]+=temp;plenv.xmin[1]-=temp;
      plenv.xmax[1]=MIN(plenv.xmax[1],spec->wl[spec->np-1]);
      plenv.xmin[1]=MAX(plenv.xmin[1],spec->wl[0]);
    }
    else if (!strncmp(pgch,"<",1)) {
      temp=.25*(plenv.xmax[1]-plenv.xmin[1]);
      plenv.xmax[1]-=temp;plenv.xmin[1]+=temp;
    }
    else if (!strncmp(pgch,"i",1)) {
      fprintf(stderr,"%05i\t%04i\t%03i\n",spec->mjd,spec->plate,spec->fiberID);
    }
  } /* End While Loop */

  return 1;
}
