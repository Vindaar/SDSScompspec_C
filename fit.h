/***************************************************************************
FIT.H: Include file for various fitting functions and routines
***************************************************************************/

/* INCLUDE FILES */

/* DEFINITIONS */
/* Functions */
#define MAX(a,b)   a>b ? a : b
#define MIN(a,b)   a<b ? a : b
#define MOD(a,b) (b==0.0) ? b : a-((int)(a/b))*b
#define SQR(a) ((a) == 0.0 ? 0.0 : a*a)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* Parameters */
#define MRQNMAX     100       /* Maximum number of mrqmin iterations */ 
#define FPMIN         1.e-30  /* Number near smallest representable float */ 
#define SVDNMAX     500       /* Maximum number of svdfit iterations */
#define SVDFIT_TOL    1.e-20  /* Tolerance level for svdfit() */
#define FIT_INFIN     1.e32   /* Effective infinity */
#define NBTOL         3       /* Bad pixel threshold for EW.c doublet fits */

/* Flags */
#define GJSNGLR      -1       /* Flag denoting that a singular gaussj() matrix 
				 was formed */
#define MRQITMAX     -2       /* Flag indicating that over MRQNMAX
				 iterations were required for fit */

/* STRUCTURES */

/* FUNCTION PROTOTYPES */
double chebyshev_eval(double x, double *coef, int ord);
int covsrt(double **covar, int ma, int *ia, int mfit);
int dpolcof(double *xa, double *ya, int n, double *cof);
int dpolint(double *xa, double *ya, int n, double x, double *y, double *dy);
int gaussj(double **a, int n, double **b, int m);
double legendre_eval(double x, double *coef, int ord);
int linfit(double *x, double *y, double *sig, int ndata, int opt1, int opt2,
	   double *a, double *b, double *siga, double *sigb, double *chi2,
	   double *q);
int linreg(double *x, double *y, int n, double *a, double *b);
int mrqfit_erffn(double *x, double *a, double **lim, int *ia, int **inf,
		 double *y, double *dyda, int nx, int na, int pix);
int mrqfit_gauss(double *x, double *a, double **lim, int *ia, int **inf,
		 double *y, double *dyda, int nx, int na, int pix);
int mrqfit_multierffn(double *x, double *a, double **lim, int *ia, int **inf,
		      double *y, double *dyda, int nx, int na, int pix);
int mrqcof(double **x, double *y, double *sig, int ndata, int nx, double *a,
	   double **lim, int *ia, int **inf, int ma, double **alpha,
	   double *beta, double *chisq,
	   int (*funcs)(double *, double *, double **, int *, int **, double *,
			double *, int, int, int));
int mrqmin(double **x, double *y, double *sig, int ndata, int nx, double *a,
	   double **lim, int *ia, int **inf, int ma, double **covar, 
	   double **alpha, double *chisq, double *alamda,
	   int (*funcs)(double *, double *, double **, int *, int **, double *,
			double *, int, int, int));
double poly_eval(double x, double *coef, int ord);
double poly_nooffset_eval(double x, double *coef, int ord);
double pythag(double a, double b);
int spline(double *x, double *y, int n, double yp1, double ypn, double *y2,
	   int opt);
double splint(double *xa, double *ya, double *y2a, int n, double x);
int svbksb(double **u, double *w, double **v, int m, int n, double *b,
	   double *x);
int svdcmp(double **a, int m, int n, double *w, double **v);
int svdfit(double *x, double *y, double *sig, int ndata, double *a, int ma,
	double ***u, double ***v, double **w, double *chisq,
	int (*funcs)(double, double *, int));
int svdfit_chebyshev(double x, double *pt, int nt);
int svdfit_legendre(double x, double *pl, int nl);
int svdfit_poly(double x, double *p, int np);
int svdfit_poly_nooffset(double x, double *p, int np);
int svdvar(double **v, int ma, double *w, double ***cvm);
