/* The command line input must specify a settings file, this program reads in the settings file which is CR separated lines of the formate [SETTING]=[VALUE],
a default settings file is created in SDSS_defaultgen.c. The settings get saved to a structure typedefed as setblock, and refered to as *set. 

If nothing is specified in [VALUE] then a default is used (set in SDSS-spec.h), the exception is the z partition and input files must be specified.
If nothing is specified in outfile then the output is to stdout.

All things work fine, though there is a bug: the most important thing read in is the zcut, this is the values of the endpoints of the zbins,
this will generally be a very long list. strtok(buffer,tokens) splits up the line buffer by the tokens, it returns the everything up to the first
 token, and then saves the rest of the line to an internal buffer, later calling strtok(NULL,tokens) calls this internal buffer. However the
 internal buffer is limited, and a portion of the zcut line might be lost. There is certainly a better way to read in, and parse this long line,
 but I can not think of it.
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "SDSS-spec.h"
#include "file.h"
#include "const.h"
#include "memory.h"
#include "stats.h"
#include "error.h"



int get_wlcut(setblock *set){
  char dummy[64]="\0";
  FILE *wlcuts;
  int i=0, max=0;
  double y[4000];

  if((wlcuts=fopen("/home/jlo/work/anisotropy/data/wlcuts.txt","r"))==NULL)
    errormsg("Could not open wlcuts file");

  while(fgets(dummy,64,wlcuts)!=NULL){
    y[i]=atof(dummy);
    i++;
  }
  (*set).nzpart=i-1;
  if(((*set).zcut=darray((*set).nzpart))==NULL)
    errormsg("Could not build Z partition in SDSS_settings");
  if(((*set).wlcut=darray((*set).nzpart))==NULL)
    errormsg("Could not build WL partition in SDSS_settings");

  for(i=0;i<(*set).nzpart;i++){
    set->wlcut[i]=y[i];
    set->zcut[i]=y[i]/LYlines[set->forest]-1;
  }
  return 1;
}



int SDSS_settings(char setfile[], setblock *set) {

  int      i=0,k=0;
  _Bool    j=0;
  char     buffer[VVHUGESTRLEN]="\0", *list=NULL,*end=NULL, *bufferLS, *bufferRS;
  char     dummy1[HUGESTRLEN]="\0",dummy2[HUGESTRLEN]="\0",*dummy3=NULL;
  double   x=0,y[1000];
  int      get_wlcut(setblock *set);
  FILE     *set_file=NULL;



  /** Open settings file for reading **/
  if ((set_file=faskropen("Valid Settings File?",setfile,5))==NULL)
    errormsg("SDSS_settings: Can not open file %s",setfile);

  while(fgets(buffer,HUGESTRLEN,set_file)!=NULL){ //Reads the line from the file


    bufferLS=strtok(buffer,"="); //Splits it up by the '='

  /* Start looking for certain phrases on the left side of the '=' */

    /* Determine which forest we are looking at */    
    if(!strcoll(bufferLS,"forest")){ 
      bufferRS=strtok(NULL,"=\n");
      if(!strcoll(bufferRS,"\0")){ set->forest=0; }
      else if(!strcoll(bufferRS,"B")){ set->forest=1; }
      else if(!strcoll(bufferRS,"G")){ set->forest=2; }
      else { set->forest=0; }
    }

    /* Set the Z partition */
    else if(!strcoll(bufferLS,"Zcut")){
      i=0;
      dummy3=strtok(NULL,"={,}\n");
      if(!strcoll(dummy3,"pixel")){
        get_wlcut(set);
        j++;
      }
      else{
        do {
          y[i] = (x=atof(dummy3));
          dummy3=strtok(NULL,"={,}\n");
          i++;
        } while(dummy3 != NULL);
        set->nzpart=i-1;
        if(set->nzpart>2) j++;
        if((set->zcut=darray(set->nzpart))==NULL)
          errormsg("Could not build Z partition in SDSS_settings");
        if((set->wlcut=darray(set->nzpart))==NULL)
          errormsg("Could not build WL partition in SDSS_settings");
        for(i=0;i<=set->nzpart-1;i++){
          set->zcut[i]=y[i];
          set->wlcut[i]=(y[i]+1)*LYlines[set->forest];
        }
      }
    }

    /* Set the Data file */
    else if(!strcoll(bufferLS,"infile")){
      bufferRS=strtok(NULL,"=\n");
      if(bufferRS==NULL)
        errormsg("No infile Specified");
      strcpy(set->infile,bufferRS);
    }

    /* Set the PE threshhold */
    else if(!strcoll(bufferLS,"Wlimit"))
      set->Wlimit = (x=atof(strtok(NULL,"=")))!=0 ? x : WLIMIT;

    /* Set the flag for plotting Spectra */
    else if(!strcoll(bufferLS,"pspec"))
      set->pspec = atoi(strtok(NULL,"="));

    /* Set the Medium Continuum */
    else if(!strcoll(bufferLS,"medcont"))
      set->medcont = (k=atoi(strtok(NULL,"="))) !=0 ? k : MEDCONT;

    /* Set Ncont Sigma */
    else if(!strcoll(bufferLS,"ncontsig"))
      set->ncontsig = (x=atof(strtok(NULL,"=")))!=0 ? x : NCONSIG;

    /* Specify an output file */
    else if(!strcoll(bufferLS,"outfile"))
      strcpy(set->outfile,(bufferRS=strtok(NULL,"=\n"))!=NULL ? bufferRS : "\0");

    /* Specifies output type */
    else if (!strcoll(bufferLS,"output"))
      set->output = atoi(strtok(NULL,"="));
  }

  if(!j)
    errormsg("No partition function specified in %s",setfile);
  if(set->output==-1){
    set->plotoutput++;
  }
  fclose(set_file);
  return 1;
}
    
