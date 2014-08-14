#include <stdlib.h>
#include <math.h>
#include "SDSS-spec.h"
#include "stats.h"
#include "error.h"
#include "memory.h"

int SDSS_statistics(compspec *cspec,spectrum *spec,int nspec, setblock *set){

  int i=0,j=0;
  int *s=NULL;
  double *u=NULL;
  statset bstat;
  FILE *params = NULL;
  params = fopen("parameters.txt", "a");

  for(i = 0; i<5763; i++){
	 fprintf(params, "%f\t%i\t%t\n", cspec[0].sum[i], cspec[0].nhist[i], cspec[0].sum2[i]);
  }

  fclose(params);

  for(j=0;j<set->nforest;j++){
    i=0;
    while(cspec[j].nhist[i] == 0) i++;
    for(;i<3000;i++){
      if(cspec[j].nhist[i]>0){
        /* UNweighted */
	cspec[j].fl[i]=cspec[j].sum[i]/cspec[j].nhist[i];
	cspec[j].er[i]=cspec[j].sum2[i]/cspec[j].nhist[i]-
	  pow(cspec[j].fl[i],2);
	cspec[j].er[i]=(cspec[j].nhist[i] > 1 ?
			sqrt(cspec[j].er[i]/(cspec[j].nhist[i]-1)) : 0);

        /* Wieghted */
	//        cspec[j].fl[i]=cspec[j].sum[i]/cspec[j].wsum[i];
	//	cspec[j].er[i]=cspec[j].sum2[i]/cspec[j].nhist[i]-
	//	  pow(cspec[j].fl[i],2);
	//	cspec[j].er[i]=(cspec[j].nhist[i] > 1 ?
	//	  sqrt(cspec[j].er[i]/(cspec[j].nhist[i]-1)) : 0);
      }
    }
  }

  printf("nspec: %i\n", nspec);
  u=darray(nspec);
  s=iarray(nspec);
  for(i=0;i<nspec;i++) {
    u[i]=spec[i].alpha;
    s[i]=(spec[i].alpha<-900 ? 0 : 1);
  }
  stats(u,NULL,NULL,NULL,s,nspec,0,&bstat);
  cspec[0].mean_a=bstat.mean;
  cspec[0].sigma_a=bstat.rms;
  median(u,s,nspec,&bstat,0);
  cspec[0].med_a=bstat.med;
  cspec[0].siqr_a=bstat.siqr;


  free(u);
  free(s);
  return 1;
}

