#include <stdlib.h>
#include "SDSS-spec.h"
#include "const.h"
#include "memory.h"
#include "error.h"

int SDSS_robustspec(compspec *robustspec, compspec *cspec, spectrum *spec, int forest) {

  int i=0,j=0;
  int goodpix=0;
  double resfactor=0;

   spec->lambmin =
     MAX((spec->zem+1)*(Restrange[forest][0]),spec->wl[0]);
   spec->lambmax= (spec->zem+1)*(Restrange[forest][1]);

   if(spec->lambmax < spec->wl[0]){
     return 9;
   }
   //   spec->lambmin = 1130*(spec->zem+1);

  /* The first time through this loop */
  if(robustspec->sidx == 0){
    robustspec->sidx=1000; // Sets the start index to 300, by default
    /* Then adds the first set of values into the block arrays */
    for(i=0;i<spec->np;i++){
      goodpix = 
       (((spec->wl[i] !=0) &&
        (spec->wl[i] + spec->wl[i-1] > 2*spec->lambmin) && 
        (spec->wl[i] + spec->wl[i+1] < 2*spec->lambmax) &&
        (spec->er[i]!=0) && 
        (abs((spec->fl[i]/spec->pc[i] - cspec->fl[robustspec->sidx + i])/
           (cspec->er1[robustspec->sidx + i]*sqrt(cspec->nhist[robustspec->sidx+i]))) < SIGCLIP))
        ? 1 : 0);
        robustspec->wl[robustspec->sidx+i] = spec->wl[i]; 
//        robustspec->dist[robustspec->sidx+i] = lumdist(spec->wl[i],forest);

      if(goodpix){

        /* UNweighted */
	//robustspec->fl[robustspec->sidx+i] = goodpix*spec->fl[i]/spec->pc[i]; 
	//robustspec->er1[robustspec->sidx+i] = pow(spec->fl[i]/spec->pc[i],2)*goodpix;
	//robustspec->er2[robustspec->sidx+i] = goodpix*spec->er[i]; 

        /* Weighted */

        robustspec->fl[robustspec->sidx+i] = spec->fl[i]/(spec->pc[i]*pow(spec->er[i],2));   //
        robustspec->er1[robustspec->sidx+i] =                                                //
          pow(spec->fl[i]/(spec->pc[i]),2)/pow(spec->er[i],2);                     //
        robustspec->er2[robustspec->sidx+i] = 1/pow(spec->er[i],2);                          //

        robustspec->nhist[robustspec->sidx+i]++;        /* Count that point */
      }

      if(spec->wl[i]+spec->wl[i+1] > 2*spec->lambmax) break;
    }

  }

  /* Every subsequent time through the loop: */
  else {
    /* If the new spectrum, goes lower than the composite spectrum, add onto the
       bottom */

    if( robustspec->wl[robustspec->sidx] - spec->wl[0] > .05 ) { 
      while(robustspec->wl[robustspec->sidx] - spec->wl[j] > .05) j++; //Count How many bins its lower by
      robustspec->sidx=robustspec->sidx-j;  // Set the new start index
      for(i=0;robustspec->wl[robustspec->sidx+i]==0;i++) { // Fill in the missing wavelengths
        robustspec->wl[robustspec->sidx+i]=spec->wl[i];
//        robustspec->dist[robustspec->sidx+i] = lumdist(spec->wl[i],forest);
      }
      j=0;  //Reset j
    }

    /* If the new spectrum is higher than bottom of the composite spectrum, add on j bins to as buffer zone */    
    else if(spec->wl[0] - robustspec->wl[robustspec->sidx] > .05){
      for(j=0; spec->wl[0] - robustspec->wl[robustspec->sidx+j] > .05 ;j++){
        ;
      }
    }
    else j=0;

    /* If something dumb happened, then everything comes crashing down and stops */
    if(abs(robustspec->wl[robustspec->sidx+j] - spec->wl[0]) > .05){
      errormsg("Something stupid:\nshould be same: %f %f\n starts: %f %f\nAlignment Error",robustspec->wl[robustspec->sidx+j],spec->wl[0],robustspec->wl[robustspec->sidx],spec->wl[0]);
    }

    /* If its all good, then we add in the new spectrum */
    else {
      for(i=0;i<spec->np;i++,j++){
        if(spec->wl[i] + spec->wl[i+1] > 2*spec->lambmax) break;
        goodpix = 
         ((spec->wl[i] + spec->wl[i-1] > 2*spec->lambmin) && 
          (spec->er[i]!=0) && 
          (abs((spec->fl[i]/spec->pc[i] - cspec->fl[robustspec->sidx + i])/
            (cspec->er1[robustspec->sidx + i]*sqrt(cspec->nhist[robustspec->sidx+i]))) < SIGCLIP))
          ? 1 : 0;
        robustspec->wl[robustspec->sidx+j] = spec->wl[i];
//        if(robustspec->dist[robustspec->sidx+i]==0)
//          robustspec->dist[robustspec->sidx+i] = lumdist(spec->wl[i],forest);

        if(goodpix){

          /* UNweighted */
	  //  robustspec->fl[robustspec->sidx+j] += spec->fl[i]/spec->pc[i];
	  //  robustspec->er1[robustspec->sidx+j] += pow(spec->fl[i]/spec->pc[i],2);
	  //  robustspec->er2[robustspec->sidx+j] += spec->er[i];

          /* Weighted */

	  robustspec->fl[robustspec->sidx+j] +=                               //
	    spec->fl[i]/(spec->pc[i]*pow(spec->er[i],2));           //
	  robustspec->er1[robustspec->sidx+j] +=                              //
            pow(spec->fl[i]/(spec->pc[i]),2)/pow(spec->er[i],2);    //
          robustspec->er2[robustspec->sidx+j] +=                              //
            1/pow(spec->er[i],2);                                   //
 

          robustspec->nhist[robustspec->sidx+j]++;         /* Count that point */
        }
      }
    }
  }
  return 1;
}

