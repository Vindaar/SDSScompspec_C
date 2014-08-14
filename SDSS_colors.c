#include "SDSS-spec.h"
#include <stdlib.h>
#include "error.h"
#include "colorcurves.h"
#include "memory.h"
#include <math.h>

#define sround(a) .0001*rint(10000.0*log10(a))
#define idx(a,b) (int) rint((log10(a)-b)/.0001)

// changed by S. Schmidt

int SDSS_colors(spectrum *spec, int pfilter) {
	/* Declarations for every time through */
	int i = 0;
	int j = 0;
	int np = 5763; //5762
	int i_temp = 0;
	double usum = 0, gsum = 0, rsum = 0, isum = 0, zsum = 0;

	/* Declarations for the First time through */
	double ucurve[NU][3] = UCURVE, gcurve[NG][3] = GCURVE, rcurve[NR][3] = RCURVE,
			icurve[NI][3] = ICURVE, zcurve[NZ][3] = ZCURVE;
			static double **a=NULL, beginwl	= 3.4742;
	double *wl = NULL;
	float *pgx = NULL, *pgy = NULL;

	if (a == NULL ) {
		a = dmatrix(np, 5);
		wl = darray(np);
		for (i = 0; i < np; i++)
			wl[i] = pow(10, beginwl + .0001 * i);
		/* Build U Curve */
		i = idx(ucurve[0][0],beginwl);
		for (j = 0; j < NU - 1; j++) {
			for (; wl[i] < ucurve[j + 1][0]; i++) {
				a[i][0] = (ucurve[j + 1][1] - ucurve[j][1])
						/ (ucurve[j + 1][0] - ucurve[j][0])
						* (wl[i] - ucurve[j][0]) + ucurve[j][1];
				//printf("%i\t%lf\n", i, a[i][0]);
			}
		}
		/* Build G Curve */
		i = idx(gcurve[0][0],beginwl);
		for (j = 0; j < NG - 1; j++) {
			for (; wl[i] < gcurve[j + 1][0]; i++) {
				a[i][1] = (gcurve[j + 1][1] - gcurve[j][1])
						/ (gcurve[j + 1][0] - gcurve[j][0])
						* (wl[i] - gcurve[j][0]) + gcurve[j][1];
			}
		}
		/* Build R curve */
		i = idx(rcurve[0][0],beginwl);
		for (j = 0; j < NR - 1; j++) {
			for (; wl[i] < rcurve[j + 1][0]; i++) {
				a[i][2] = (rcurve[j + 1][1] - rcurve[j][1])
						/ (rcurve[j + 1][0] - rcurve[j][0])
						* (wl[i] - rcurve[j][0]) + rcurve[j][1];
			}
		}
		/* Build I Curve */
		i = idx(icurve[0][0],beginwl);
		for (j = 0; j < NI - 1; j++) {
			for (; wl[i] < icurve[j + 1][0] && i < np; i++) {
				a[i][3] = (icurve[j + 1][1] - icurve[j][1])
						/ (icurve[j + 1][0] - icurve[j][0])
						* (wl[i] - icurve[j][0]) + icurve[j][1];
			}
		}
		/* Build Z Curve */
		i = idx(zcurve[0][0],beginwl);
		for (j = 0; j < NZ - 1; j++) {
			for (; wl[i] < zcurve[j + 1][0] && i < np; i++) {
				a[i][4] = (zcurve[j + 1][1] - zcurve[j][1])
						/ (zcurve[j + 1][0] - zcurve[j][0])
						* (wl[i] - zcurve[j][0]) + zcurve[j][1];
			}
		}
		//    free(wl);
		if (pfilter) {
			cpgopen("?");
			cpgenv(wl[0], wl[np - 1], 0, .6, 0, 0);
			cpglab("Wavelength (\\A)", "Throughput",
					"Throughput For SDSS Filters");
			pgx = farray(np);
			pgy = farray(np);
			for (j = 0; j < 5; j++) {
				for (i = 0; i < np; i++) {
					pgx[i] = wl[i];
					pgy[i] = a[i][j];
				}
				cpgsci(j + 3);
				cpgline(np, pgx, pgy);
			}
			cpgclos();
			}
	}
	// Changed from else if( spec == NULL ) {
	if ( spec == NULL ) {
		free(*a);
		free(a);
		free(wl);
		return 1;
	}
	j = idx(spec->wl[0],beginwl);
	for (i = 0; i < spec->np && j < np; i++, j++) {
		usum += a[j][0] * spec->fl[i];
		gsum += a[j][1] * spec->fl[i];
		rsum += a[j][2] * spec->fl[i];
		isum += a[j][3] * spec->fl[i];
		zsum += a[j][4] * spec->fl[i];
	}

	spec->smag[0] = 22.5 - 2.5 * log10(usum);
	spec->smag[1] = 22.5 - 2.5 * log10(gsum);
	spec->smag[2] = 22.5 - 2.5 * log10(rsum);
	spec->smag[3] = 22.5 - 2.5 * log10(isum);
	spec->smag[4] = 22.5 - 2.5 * log10(zsum);

	free(wl); //memory leak, fixed

	return 1;
}
