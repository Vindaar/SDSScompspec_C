/****************************************************************************
* Initialize plotting environment
****************************************************************************/

#include <string.h>
#include "SDSS-spec.h"

int SDSS_pgenv_init() {

  extern plotenv plenv;
  
  plenv.wwidth=W_WIDTH;
  plenv.wasp=W_ASP;
  plenv.ch=C_H;
  plenv.lw=L_W;
  plenv.nxsub=N_X_SUB;
  plenv.nysub=N_Y_SUB;
  plenv.vpu=VPU;
  plenv.vpd=VPD;
  plenv.vpl=VPL;
  plenv.vpr=VPR;

  strcpy(plenv.xlab[0],"");
  strcpy(plenv.ylab[0],"");
  strcpy(plenv.title[0],"");

  /* Arrow head style */
  cpgsah(1,45.0,0.3);

  return 1;
}
