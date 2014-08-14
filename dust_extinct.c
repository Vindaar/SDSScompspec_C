/******************************************************************************
 * Mean extinction law A(l) / A(V) (l in Angstrom)
 * R = A(V) / E(B-V)
 *
 * Cardelli et al., 1989, ApJ, 345, 245 (MW: Far-UV, UV and NIR)
 * O'Donnell, 1994, ApJ, 422, 158 (MW: optical)
 * Pei, 1992, ApJ, 395, 130 (MW, LMC, SMC)
 *
 * Options:
 * 0 = MW from CCM + O'Donnell
 * 1 = MW from Pei  -|
 * 2 = LMC from Pei  |- R is not a free parameter and the input value is ignored.
 * 3 = SMC from Pei -|
 * 4 = empirical MW from Pei   R should be = 3.08
 * 5 = empirical LMC from Pei  R should be = 3.16
 * 6 = empirical SMC from Pei  R should be = 2.93

 For option 0, there are not empirical data below 1000A or above 33um
 so the value of A(V)/A(l) for any wavelengths entered beyond these
 ranges simply represent extrapolations of the CCM fitting
 functions. They should not be treated as reliable.
 ******************************************************************************/

#include <math.h>
#include "error.h"

double dust_extinct(double l, double r, int opt) {

	double x = 0.0, y = 0.0, a = 0.0, b = 0.0;
	double xx1 = 0.0, xx2 = 0.0, xx3 = 0.0, xx4 = 0.0, xx5 = 0.0, xx6 = 0.0;
	int i = 0, j = 0;
#include "dust_extinct.h"

	/* 1/lambda in um^-1 */
	x = 1.0e4 / l;

	if (opt == 0) {
		/* CCM and O'Donnell */
		if (r == 0.0) {
			nferrormsg("dust_extinct(): R must be non-zero: %lf", r);
			return -1.0;
		}
		if (x <= 1.1) {
			/* NIR */
			xx1 = pow(x, 1.61);
			a = 0.574 * xx1;
			b = -0.527 * xx1;
		} else if (x <= 3.3) {
			/* Optical */
			y = x - 1.82;
			xx1 = y * y * y;
			xx2 = xx1 * y;
			xx3 = xx2 * y;
			xx4 = xx3 * y;
			xx5 = xx4 * y;
			xx6 = xx5 * y;
			a = 1.0 + 0.104 * y - 0.609 * y * y + 0.701 * xx1 + 1.137 * xx2
					- 1.718 * xx3 - 0.827 * xx4 + 1.647 * xx5 - 0.505 * xx6;
			b = 1.952 * y + 2.908 * y * y - 3.989 * xx1 - 7.985 * xx2
					+ 11.102 * xx3 + 5.491 * xx4 - 10.805 * xx5 + 3.347 * xx6;
		} else if (x <= 8.0) {
			/* UV */
			xx1 = x - 4.67;
			xx2 = xx1 * xx1;
			xx3 = x - 4.62;
			xx4 = xx2 * xx2;
			a = 1.752 - 0.316 * x - 0.104 / (xx2 + 0.341);
			b = -3.09 + 1.825 * x + 1.206 / (xx4 + 0.263);
			if (x >= 5.9) {
				xx3 = x - 5.9;
				xx4 = xx3 * xx3;
				xx5 = xx4 * xx3;
				a += -0.04473 * xx4 - 0.009779 * xx5;
				b += 0.213 * xx4 + 0.1207 * xx5;
			}
		} else {
			/* Far UV */
			xx1 = x - 8.0;
			xx2 = xx1 * xx1;
			xx3 = xx2 * xx1;
			a = -1.073 - 0.628 * xx1 + 0.137 * xx2 - 0.07 * xx3;
			b = 13.67 + 4.257 * xx1 - 0.42 * xx2 + 0.374 * xx3;
		}

		return (a + b / r);

	}

	if (opt >= 1 && opt <= 3) {
		/* Pei */
		/* These parametrizations are a bit bollocks since they are not
		 normalised properly */

		j = opt - 1;
		a = 0.0;
		for (i = 0; i < 5; i++) {
			xx1 = pow(1.0 / (x * lp[j][i]), np[j][i]);
			xx2 = pow(x * lp[j][i], np[j][i]);
			a += ap[j][i] / (xx1 + xx2 + bp[j][i]);
		}

		return a;

	}

	if (opt >= 4 && opt <= 6) {
		/* Pei empirical laws */
		if (r == 0.0) {
			nferrormsg("dust_extinct(): R must be non-zero: %lf", r);
			return -1.0;
		}
		j = opt - 4;
		if (x < lam[j][0] || x > lam[j][29]) {
			nferrormsg("dust_extinct(): Wavelength out of range: %f", l);
			return -1.0;
		}
		i = 1;
		while (i < 30 && x > lam[j][i])
			i++;
		a = (eleb[j][i] - eleb[j][i - 1]) / (lam[j][i] - lam[j][i - 1]);
		y = eleb[j][i] + a * (x - lam[j][i]);

		return y / r + 1.0;

	}

	nferrormsg("dust_extinct(): Unknown option: %d", opt);
	return -1.0;

}
