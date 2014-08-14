/****************************************************************************

SDSS_primtarg: Assigns flags to target selection flags based on
already-set value for primtarg value.

****************************************************************************/

#include "SDSS-spec.h"
#include "error.h"

int SDSS_primtarg(targselflag *tsf) {

  //static int count=0;

  /* Check primtarg flag */
  if (tsf->primtarg==-1) {
    nferrormsg("SDSS_primtarg(): Primary target integer is not set");
    return 0;
  }

  /* Set integer selection flags from hex bit flags */
  tsf->qso_hiz=(tsf->primtarg & TARG_QSO_HIZ) ? 1 : 0;
  tsf->qso_cap=(tsf->primtarg & TARG_QSO_CAP) ? 1 : 0;		
  tsf->qso_skirt=(tsf->primtarg & TARG_QSO_SKIRT) ? 1 : 0;
  tsf->qso_first_cap=(tsf->primtarg & TARG_QSO_FIRST_CAP) ? 1 : 0;
  tsf->qso_first_skirt=(tsf->primtarg & TARG_QSO_FIRST_SKIRT) ? 1 : 0;
  tsf->qso_faint=(tsf->primtarg & TARG_QSO_FAINT) ? 1 : 0;
  tsf->qso_reject=(tsf->primtarg & TARG_QSO_REJECT) ? 1 : 0;
  tsf->galaxy_red=(tsf->primtarg & TARG_GALAXY_RED) ? 1 : 0;
  tsf->galaxy_red_ii=(tsf->primtarg & TARG_GALAXY_RED_II) ? 1 : 0;
  tsf->galaxy=(tsf->primtarg & TARG_GALAXY) ? 1 : 0;
  tsf->galaxy_big=(tsf->primtarg & TARG_GALAXY_BIG) ? 1 : 0;	
  tsf->galaxy_bright_core=(tsf->primtarg & TARG_GALAXY_BRIGHT_CORE) ? 1 : 0;
  tsf->rosat_a=(tsf->primtarg & TARG_ROSAT_A) ? 1 : 0;
  tsf->rosat_b=(tsf->primtarg & TARG_ROSAT_B) ? 1 : 0;
  tsf->rosat_c=(tsf->primtarg & TARG_ROSAT_C) ? 1 : 0;
  tsf->rosat_d=(tsf->primtarg & TARG_ROSAT_D) ? 1 : 0;
  tsf->rosat_e=(tsf->primtarg & TARG_ROSAT_E) ? 1 : 0;
  tsf->star_bhb=(tsf->primtarg & TARG_STAR_BHB) ? 1 : 0;
  tsf->star_carbon=(tsf->primtarg & TARG_STAR_CARBON) ? 1 : 0;
  tsf->star_brown_dwarf=(tsf->primtarg & TARG_STAR_BROWN_DWARF) ? 1 : 0;
  tsf->star_sub_dwark=(tsf->primtarg & TARG_STAR_SUB_DWARK) ? 1 : 0;
  tsf->star_caty_var=(tsf->primtarg & TARG_STAR_CATY_VAR) ? 1 : 0;
  tsf->star_red_dwarf=(tsf->primtarg & TARG_STAR_RED_DWARF) ? 1 : 0;
  tsf->star_white_dwarf=(tsf->primtarg & TARG_STAR_WHITE_DWARF) ? 1 : 0;
  tsf->star_pn=(tsf->primtarg & TARG_STAR_PN) ? 1 : 0;
  tsf->serendip_blue=(tsf->primtarg & TARG_SERENDIP_BLUE) ? 1 : 0;
  tsf->serendip_first=(tsf->primtarg & TARG_SERENDIP_FIRST) ? 1 : 0;
  tsf->serendip_red=(tsf->primtarg & TARG_SERENDIP_RED) ? 1 : 0;
  tsf->serendip_distant=(tsf->primtarg & TARG_SERENDIP_DISTANT) ? 1 : 0;
  tsf->serendip_manual=(tsf->primtarg & TARG_SERENDIP_MANUAL) ? 1 : 0;

  /* Is the QSO colour-selected? */
  tsf->colselect=(tsf->qso_hiz || tsf->qso_cap || tsf->qso_skirt) ? 1 : 0;

  /* Is the QSO selected as a FIRST QSO? */
  tsf->first=(tsf->qso_first_cap || tsf->qso_first_skirt) ? 1 : 0;

  return 1;

}
