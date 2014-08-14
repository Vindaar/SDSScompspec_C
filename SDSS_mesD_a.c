/****************************************************
Dsum, deltawl, and meanwl are (refered to as the sums) running sum which
adds the Da values for points, the width of the points, and the wavelength of
the pixels. Each pixel is weighted by the width of the pixel, so the mean 
opacity would be something like Sum(Da_i * deltawl_i)/Sum(deltawl_i)
where _i indicates of the ith pixel, in the program this would be written
as Dsum/deltawl.

There is also a flag called badflag, this indicates that a pixel is either a 
no data pixel, or a cosmic ray, the check for the former is if the error array drops to 0, and the check for the latter is if the error at that pixel is > zstar*errbar
where errbar is the mean value of the error array for the surrounding 60 pixels, and
zstar is an arbitrary number decided by the user, 3 works well.

****************************************************/


#include "SDSS-spec.h"
#include "const.h"
#include "error.h"


int FindBin(double wlcut[], double wlx, int nzpart) {
  /* This program sorts through the bins looking for which bin the top
     pixel lies in    */
  int i=0;
  for(i=0;wlcut[i+1]<=wlx && i<=nzpart-2;i++)
    ;
  if (wlx >= wlcut[i])
    return i;
  else
    return -1;
}

int SDSS_mesD_a(spectrum *spec, double wlcut[], int nzpart, int forest){

   int i=0, j=0, k=0, m=0, topbin=0;
   long double Dsum=0;
   int badflag=0,bin=0,errbin=30;
   double deltawl=0,errbar=0,zstar=3,meanwl=0, p=0;

   
   /* Set the bottom of the calculated forest */
   spec->lambmin =
     MAX((spec->zem+1)*(LYlines[forest+1]+emislinebuffer[forest+1]),spec->wl[0]);
 
  /* Initialise arrrays */
  if((spec->opacity=darray(nzpart))==NULL)
    errormsg("Could not build D_a array for file \n\t%s",spec->file);

  if((spec->deltawl=darray(nzpart))==NULL)
    errormsg("Could not build delta lambda array for file \n\t%s",spec->file);
  if((spec->meanz=darray(nzpart))==NULL)
    errormsg("Could not build mean z array for file \n\t%s",spec->file);

  /* Exit the program if no bounds were measured */
  if(spec->errflg){
    nferrormsg("D_a not measured due to previous error in file \n\t%s", 
                spec->file);
    return 0;
    }
  /* Exit if a Null range is specified */
  if(spec->lambmin >= spec->lambmax){
    nferrormsg("Bounds specify a NULL range in file\n\t%s", spec->file);
    spec->errflg++;
    return 0;
  }     

  /* Exit the routine if forest lies above all the bins */
  if(spec->lambmin > wlcut[nzpart-1])
    return 1;

  topbin = FindBin(wlcut,spec->lambmax,nzpart);

 /* Exit the routine if the top pixel lies below all the bins */
  if(topbin == -1)
    return 1;


  bin=topbin;

  /* Start from the LYA line (excluding the PE zone) and count down pixel by pixel
     checking to see that each pixel has good data, (is not a Cosmic Ray) and 
     that the pixel has not switched bins */

  for(i=spec->lambmax_ind-errbin;i<=spec->lambmax_ind+errbin;i++) // See note below
    errbar+=spec->er[i]/(2*errbin+1);                             //*************** 



  for(i=spec->lambmax_ind;spec->wl[i] >= spec->lambmin; i--){

/*********************************************************************************
   This is a check on the variance to determine a cosmic ray or section of no data
   it flags an array called spec->bp, that aligns with the wl,fl etc arrays. This 
   routine is only run in the range of the D_a calculations. But with some editing
   can be applied to the entire spectrum. It is set up to handle running out of 
   data from the bottom, but not the top. It will also encounter serious difficulty
   in sky absorption regions.
**********************************************************************************/
    spec->bp[i]=0;
    /* Check for no data pixel */
    if(spec->er[i]==0) 
      spec->bp[i]++;
    /* Decent check for cosmic ray */
    if(spec->er[i]/errbar >= zstar)
      spec->bp[i]++;
    
/*********************************************************************************/


    /* If the pixel is a good data point, add the Da value, the width of the pixel,
       and the wl of the pixel to the sums, if a bad pixel, skip it */
    Dsum+=(badflag==0 ? (1-spec->fl[i]/spec->pc[i])*(spec->wl[i+1]-spec->wl[i]) : 0);
    deltawl += (badflag==0 ? spec->wl[i+1]-spec->wl[i] : 0);
    meanwl += (badflag==0 ? spec->wl[i]*(spec->wl[i+1]-spec->wl[i]) : 0);

    /* If the point i-1 is in the lower bin do this: */
    if(spec->wl[i-1] < wlcut[bin]){

      /* If there is a pixel below, */
      if(i>0){ 
        i--; //Go to that pixel
        /* Add the proportion of that pixel that lies in the current bin in */
        p=(spec->wl[i+1]-wlcut[bin])/(spec->wl[i+1]-spec->wl[i]); 
        Dsum+=(badflag==0 ?(1-spec->fl[i]/spec->pc[i])*(spec->wl[i+1]-spec->wl[i])*p:0);
        deltawl += (badflag==0 ? (spec->wl[i+1]-spec->wl[i])*p : 0);
        meanwl += (badflag==0 ? spec->wl[i]*(spec->wl[i+1]-spec->wl[i])*p : 0);
      }  
      else p=1;

      /* If there is some good data in this bin, find the means, and save it to spec*/
      if(deltawl>0){
        spec->opacity[bin] = Dsum/(deltawl);
        spec->deltawl[bin]=deltawl;
        spec->meanz[bin]=(meanwl/deltawl)/LYlines[forest]-1;
      }

      /* Reset the sums to 0, but add the proportion of the current pixel that didnt
         get added to the previous bin, if this there is no pixel below the previous
         one, then p is set to 1, meaning 0% of the lower pixel, which doesnt exist
         will be added to the new bin. But it doesnt matter anyway, since this was the
         last bin, its just a housekeeping thing, since this loop should really be
         a do loop */

      deltawl=(badflag==0 ? (spec->wl[i+1]-spec->wl[i])*(1-p) : 0);
      Dsum=(badflag==0 ? (1-spec->fl[i]/spec->pc[i])*(spec->wl[i+1]-spec->wl[i])*(1-p)
        : 0);
      meanwl=(badflag==0 ? spec->wl[i]*(spec->wl[i+1]-spec->wl[i])*(1-p) : 0);
      bin--;

      /* If moving down one bins means we are out of bins, then end the routine */
      if(bin==-1)
        return 1;
    }
    /* The errbar is a measure of the mean error 30 pixels on either side,
       (it is unweighted), but when the 30 pixel bin runs into the bottom
       of the spectrum, it gets shifted up to 29 below 31 above, 28 below
       32 above. etc.. */ 
    if(i-errbin+k<0)
      k++;
    errbar+=spec->er[i-errbin+k]/(2*errbin+1)-spec->er[i+errbin-k]/(2*errbin+1);
  }

    /* If we run out of pixels before we run out of bins, this for loop will end
       and the remaining values in the sums will get swept into the correct bins
       before exiting the routine */
    spec->opacity[bin] = deltawl>0 ? Dsum/(deltawl) : 0;
    spec->deltawl[bin]= deltawl >0 ? deltawl : 0;// Does this line seem stupid?
    spec->meanz[bin]=deltawl>0 ? (meanwl/deltawl)/LYlines[forest]-1 : 0;

    /* This is to save the last set of points calculated, really this should 
       eventually be replaced with a do loop. */




   return 1;
}
