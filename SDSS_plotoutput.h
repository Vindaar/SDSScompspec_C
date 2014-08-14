/* INCLUDE FILES */

/* DEFINITIONS */


/* STRUCTURES */
typedef struct Stats{
  double *mean;
  double *zbar;
  double *stddev;
  int *nspectra;
  int npoints;
} stats;

typedef struct Dipole{
  _Bool shift;
  _Bool split;
  double THETA;
  double PHI;
  double BETA;
} dipole;


/* Prototypes */
  int findstats(spectrum *spec[],int nspec,setblock *set,stats stat[], dipole dip);
