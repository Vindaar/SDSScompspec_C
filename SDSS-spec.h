/***************************************************************************/
/* Definitions, structures and function prototypes for SDSS-spec		   */
/***************************************************************************/

/* INCLUDE FILES */
#include <fitsio.h>
#include <longnam.h>
#include "pg_plot.h"

/* DEFINITIONS */
/* Functions */
#define MAX(a,b)   a>b ? a : b
#define MIN(a,b)   a<b ? a : b
#define RMDR(a,b)  (b==0) ? b : a-((int)(a/b))*b
#define ABSF(a)    a<0 ? -a : a

/* General parameters */
#define    VERSION   2.10     /* Program version                             */
#define    AUTHOR    "Jon Ouellet" /* File Author                            */
#define    DIR_PERM  00755    /* Permission code for creating new directories*/
#define    INFIN     1.e32    /* Effectively infinite value                  */
#define    NFWHM     4.e18    /* Proportionality between N(HI) & FWHM        */
                              /* (Jenkins 1971)                              */
#define    MINSNR    .1      /* Min usable SNR for pixels */
/* Spectrum parameters */
#define    RESPOW    2000.00  /* Resolving power of SDSS spectrograph       */

/* SDSS effective band wavelengths, Richards et al. (2001, AJ, 121, 2308)   */
#define    SDSSEWLU  3651.0   /* Effective u'-band wavelength               */
#define    SDSSEWLG  4679.0   /* Effective g'-band wavelength               */
#define    SDSSEWLR  6175.0   /* Effective r'-band wavelength               */
#define    SDSSEWLI  7494.0   /* Effective i'-band wavelength               */
#define    SDSSEWLZ  8873.0   /* Effective z'-band wavelength               */

/* Emission lines */
#define    NEMFREE      4   /* Number of regions free of emission lines   */
#define    TOPCUT       1.5
#define    LOWCUT       -2
//#define    TOPCUT       15
//#define    LOWCUT       -20


/* Target selection flags */
#define TARG_QSO_HIZ            0x1
#define TARG_QSO_CAP            0x2
#define TARG_QSO_SKIRT          0x4
#define TARG_QSO_FIRST_CAP      0x8
#define TARG_QSO_FIRST_SKIRT    0x10
#define TARG_QSO_FAINT          0x2000000
#define TARG_QSO_REJECT         0x20000000
#define TARG_GALAXY_RED         0x20
#define TARG_GALAXY_RED_II      0x4000000
#define TARG_GALAXY             0x40
#define TARG_GALAXY_BIG         0x80
#define TARG_GALAXY_BRIGHT_CORE 0x100
#define TARG_ROSAT_A            0x200
#define TARG_ROSAT_B            0x400
#define TARG_ROSAT_C            0x800
#define TARG_ROSAT_D            0x1000
#define TARG_ROSAT_E            0x8000000
#define TARG_STAR_BHB           0x2000
#define TARG_STAR_CARBON        0x4000
#define TARG_STAR_BROWN_DWARF   0x8000
#define TARG_STAR_SUB_DWARK     0x10000
#define TARG_STAR_CATY_VAR      0x20000
#define TARG_STAR_RED_DWARF     0x40000
#define TARG_STAR_WHITE_DWARF   0x80000
#define TARG_STAR_PN            0x10000000
#define TARG_SERENDIP_BLUE      0x100000
#define TARG_SERENDIP_FIRST     0x200000
#define TARG_SERENDIP_RED       0x400000
#define TARG_SERENDIP_DISTANT   0x800000
#define TARG_SERENDIP_MANUAL    0x1000000

/* Continuum fitting */
#define    MEDCONT     41     /* Number of pixels for median continuum      */
#define    NCONSIG      2.0   /* Rejection threshold for cont. fit. [sigma] */
#define    SIGCLIP      1.0    /* Sig clip value for robust test */

/* Balnicity flags */
#define    BAL_REJ     0      /* Rejected because spectrum doesn't cover CIV*/
#define    BAL_BAD     1      /* Bad data: >2000km/s of bad data near CIV   */
#define    BAL_ACC     2      /* Accepted: balnicity calculation is possible*/
#define    BAL_SNR     3      /* 75% threshold has been mod. due to poor SNR*/

/* Dust maps */
#define    DUSTMAPS    "/home/basti/SDSS_indie/dust_maps/"
                              /* Dust maps installation directory           */
#define    DUSTNMAP    "maps/SFD_dust_4096_ngp.fits"
                              /* Nothern dust map location                  */
#define    DUSTSMAP    "maps/SFD_dust_4096_sgp.fits"
                              /* Southern dust map loc/home/hosseinf/usr/Anisotropy+Shelf/Machinery/SDSSprogram/dr7_list.lisation                 */
/* D_a Analysis */
#define    WLIMIT      .03    /* Default PE Threshhold                      */
#define LYMANS {1215.67,1025.72,972.537,949.743 }
#define RESTRANGE {{1095,1150},{982,1010},{955,968},{912,947}}
//#define RESTRANGE{{1440.0,1475.0}}

/* Emission lines */
#define    LYA       1215.67
#define    LYB       1025.72
#define    LYG        972.537
#define    LYD        949.743
#define    LYE        937.803
#define    LYL        912.00
#define    CIV       1545.86
#define    MGII      2800.32


/* Plotting parameters */
#define    W_WIDTH     11.0
#define    W_ASP        0.75
#define    C_H          1.4
#define    L_W          3.0
#define    VPU          0.90
#define    VPD          0.15
#define    VPL          0.10
#define    VPR          0.95
#define    N_X_SUB      1
#define    N_Y_SUB      1

/* STRUCTURES */
typedef struct TargSelFlag {
  int primtarg;                  /* Prim, target (decimal version of hex bit)*/
  int colselect;                 /* Is the QSO colour selected?              */
  int first;                     /* Is the QSO selected a FIRST QSO?         */
  int qso_hiz;                   /* Target selection flag: 1 (yes) or 0 (no) */
  int qso_cap;                   /* Target selection flag: 1 (yes) or 0 (no) */
  int qso_skirt;                 /* Target selection flag: 1 (yes) or 0 (no) */
  int qso_first_cap;             /* Target selection flag: 1 (yes) or 0 (no) */
  int qso_first_skirt;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int qso_faint;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int qso_reject;                /* Target selection flag: 1 (yes) or 0 (no) */
  int galaxy_red;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int galaxy_red_ii;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int galaxy;                    /* Target selection flag: 1 (yes) or 0 (no) */
  int galaxy_big;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int galaxy_bright_core;	 /* Target selection flag: 1 (yes) or 0 (no) */
  int rosat_a;			 /* Target selection flag: 1 (yes) or 0 (no) */
  int rosat_b;			 /* Target selection flag: 1 (yes) or 0 (no) */
  int rosat_c;			 /* Target selection flag: 1 (yes) or 0 (no) */
  int rosat_d;			 /* Target selection flag: 1 (yes) or 0 (no) */
  int rosat_e;			 /* Target selection flag: 1 (yes) or 0 (no) */
  int star_bhb;			 /* Target selection flag: 1 (yes) or 0 (no) */
  int star_carbon;               /* Target selection flag: 1 (yes) or 0 (no) */
  int star_brown_dwarf;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int star_sub_dwark;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int star_caty_var;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int star_red_dwarf;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int star_white_dwarf;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int star_pn;			 /* Target selection flag: 1 (yes) or 0 (no) */
  int serendip_blue;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int serendip_first;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int serendip_red;              /* Target selection flag: 1 (yes) or 0 (no) */
  int serendip_distant;		 /* Target selection flag: 1 (yes) or 0 (no) */
  int serendip_manual;		 /* Target selection flag: 1 (yes) or 0 (no) */
} targselflag;

/* Various Settings to be used in the program */
typedef struct Setblock {
  char        infile[NAMELEN];    /* The input file                           */
  char        outfile[NAMELEN];   /* The output file                          */
  char        notes[VLNGSTRLEN];  /* Optional Notes field                     */
  int         medcont;            /* Median Fitting length                    */
  double      ncontsig;
  int         pspec;              /* Flag for plotting Spectra                */
  int         pfilter;            /* Flag for plotting filter curves          */
  int         runtime;
  int         silent;             /* Supress status output                    */
  int         supress;            /* Supress the final output                 */
  int         forest;             /* Selects which forest we are measuring    */
  int         nocomp;             /* Flag for not building a Composite Spec   */
  int         restframe;          /* Flag for doing calculation in quasar frame */
  int         synflux;            /* Flag for using the synflux instead of PC */
  int         hdu;
  int         nforest;            /* Number of Forests derived in the program */
  int         noext;              /* Turn off Galactic reddening Correction   */
  int         metal_corr;         /* Turn on metal corrections                */
  char        keep[NAMELEN];
  char        restspec[NAMELEN];
  double      LYlines[4];
  int 	      debug;              /* Option to activate debugging information */
} setblock;

typedef struct Restspec {
  double minwl;
  double maxwl;
  double *wl;
  double *fl;
  int    np;
} restspec;

typedef struct Compspec {
  double      *wl;                /* Pointer to composite Wavelength Array */
  double      *sum;               /* Sum wi*xi */
  double      *sum2;              /* Sum wi*(xi^2) */
  double      *wsum;              /* Sum wi */
  double      *wsum2;             /* Sum wi^2 */
  double      *fl;                /* Pointer to composite Flux Array */
  double      *er;               /* Pointer to Composite Flux Error Array */
  double      mean_a;              /* Mean spectral index */
  double      sigma_a;            /* Standard dev of spectral indicies           */
  double      med_a;
  double      siqr_a;
  int      *nhist;
  int      sidx;
  int      eidx;
} compspec;

typedef struct Spectrum {
  int         plate;	          /* ID of the Plate                             */
  int         fiberID;            /* Optical Fiber Number                        */
  int         usespec;            /* Flag for contribution to DA                 */
  double      R;                  /* Resolving power of spectrum                 */
  double      ra;                 /* Right ascension of object                   */
  double      dec;                /* Declination of object                       */
  double      l;                  /* Galactic latitude of object                 */
  double      b;                  /* Galactic longitude of object                */
  double      mag[5];             /* PSF magnitudes in ugriz bands (corrected for*/
                                  /*   Galactic exinction)                       */
  double      emag[5];            /* Errors in PSF magnitudes in ugriz bands     */
  double      smag[5];            /* Spectroscopic Magnitudes in ugriz bands     */
  double      ext[5];             /* Galactic extinction in ugriz bands          */
  double      Ebv;                /* Galactic reddening or colour excess E(B-V)  */
  double      zem;                /* Emission redshift of QSO                    */
  double      zs;                 /* Starting redshift (used for abs. search     */
  double      ze;                 /* Ending redshift                             */
  double      beginwl;            /* Wavelength of first pixel (logarithmic)     */
  double      deltwl;             /* Dispersion from header (logaithmic)         */
  double      cpix;               /* Indexing pixel from header                  */
  double      disp;               /* Dispersion in km/s                          */
  double      alpha;              /* f_nu(nu) spectral index redwards of Lya     */
  double      ealpha;             /* Error in spectral index redwards of Lya     */
  double      chisq;              /* Reduced Chi^2 of fit                        */
  double      beta;               /* f_lambda(lamda) spectral index              */
  double      delta;              /* Offset from fit which provides beta         */
  double      cbalnic;            /* CIV balnicity                               */
  double      mgbalnic;           /* MgII balnicity                              */
  double      globmax;
  double      globmin;
  double      emfree[3][NEMFREE]; /* Spectrum of emission-line free regions      */
  double      lambmax;            /* Top of the range        */
  double      lambmin;            /* Bottom of range         */
  double      *wl;                /* Pointer to wavelength array                 */
  double      *fl;                /* Pointer to flux array                       */
  double      *er;                /* Pointer to error array                      */
  double      *co;                /* Pointer to median-filter continuum array    */
  double      *pc;                /* Pointer to power-law continuum array        */
  double      *sn;                /* Pointer to SNR array                        */
  double      *ps;                /* Pointer to power spectrum array             */
  double      *eps;               /* Pointer to Power spectrum error array       */
  int         fidx;               /* File index                                  */
  int         mjd;                /* MJD of Observation                          */
  int         np;                 /* Number of pixels in arrays                  */
  int         nem;                /* Number of emission lines                    */
  int         bal;                /* Is this QSO a BAL by eye or not?            */
  int         cbal;               /* Can CIV balnicity be determined?            */
  int         mgbal;              /* Can MgII balnicity be determined?           */
  int         *st;                /* Pointer to pixel status array               */
  int         lambmax_ind;        /* Index of lambmax                            */
  int         lambmin_ind;        /* Index of lambmin                            */
  _Bool       errflg;             /* A flag that says if an error occured        */
  int         PlateSpec;          /* 0 if file is spSpec, 1 if file is spPlate,
                                     2 if it is spec                             */
  char        file[LNGSTRLEN];    /* Name of file                                */
  char        abfile[NAMELEN];    /* Name of file without full path              */
  char        ra_str[NAMELEN];    /* String form of RA                           */
  char        dec_str[NAMELEN];   /* String form of DEC                          */
  char        obj[FLEN_KEYWORD];  /* Object name                                 */
//  char        time[NAMELEN];      /* Time                                        */
//  char        dummy1[HUGESTRLEN]; /* Used for completeness of input line         */
  targselflag tsf;                /* Target selection flags                      */
} spectrum;


// ==================== GLOBAL VARIABLE ==============================================

extern setblock set;
extern int abc;


/* FUNCTION PROTOTYPES */
#ifdef LINUX
void    absorption_lines__();
void    spabs_init__();
#else
void    absorption_lines_();
void    spabs_init_();
#endif
double djmax(double array[], int n);
double djmin(double array[], int n);
double dust_extinct(double l, double r, int opt);
int fhist(float *data, int ndat, float *hist, float *x, int nhist, int opt1,
          int opt2);
float fjmax(float array[], int n);
float fjmin(float array[], int n);
double gauss(double x);
void get_input(char *query, char *fmt, ...);
int idxdmax(double *array, int n);
int idxdmin(double *array, int n);
int idxdval(double *array, int n, double val);
int idxfmax(float *array, int n);
int idxfval(float *array, int n, float val);
double ran(long *idum);
double res_func(double wl,int spectro);
int SDSS_CIVbalnicity(spectrum *spec);
int SDSS_colors(spectrum *spec, int pfilter);
int SDSS_compfits(compspec *cspec, spectrum *spec,int nspec,int speccount,
		  setblock *set);
int SDSS_compspec(compspec *cspec, spectrum *spec, setblock *set);
int SDSS_compinit(compspec **cspec,setblock *set);
int SDSS_Ebv(spectrum *spec);
int SDSS_Gal_extinct_correct(spectrum *spec);
int SDSS_MgIIbalnicity(spectrum *spec);
int SDSS_pgenv_init();
int SDSS_plot_spec(spectrum *spec);
int SDSS_plotoutput(compspec *cspec);
int SDSS_primtarg(targselflag *tsf);
int SDSS_red_powerlaw(spectrum *spec);
int SDSS_restspec(spectrum *spec,restspec *rspec);
int SDSS_restinit(restspec **rspec,setblock *set);
int SDSS_rinputfile(spectrum **spec, int *nspec, char *infile);
int SDSS_robustspec(compspec *robustspec,compspec *cspec, spectrum *spec, int forest);
int SDSS_rspSpec(spectrum *spec,setblock *set);
int SDSS_rplate(spectrum *spec, setblock *set);
int SDSS_rspec(spectrum *spec,setblock *set);
int SDSS_sdss_cont(spectrum *spec);
int SDSS_settings(char setfile[], setblock *set);
int SDSS_statistics(compspec *cspec,spectrum *spec,int nspec,setblock *set);



