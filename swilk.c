/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2000-2016   The R Core Team.
 *  Copyright (C) 2024        Barcelona Supercomputing Center.
 *
 *  Based on Applied Statistics algorithms AS181, R94
 *    (C) Royal Statistical Society 1982, 1995
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.

 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

/* swilk.f -- translated by f2c (version 19980913).
 * ------- and produced by f2c-clean,v 1.8 --- and hand polished: M.Maechler
 *
 * 2024-04-18: Rodrigo Arias Mallo -- Return -1 on error.
 */

#include "swilk.h"

#include <math.h>
#include <errno.h>
#include "common.h"
#include "qnorm.h"
#include "pnorm.h"

#define MIN(a, b) ((a) > (b) ? (b) : (a))

/* Algorithm AS 181.2	Appl. Statist.	(1982) Vol. 31, No. 2
 *
 * Calculates the algebraic polynomial of order nord-1 with array of
 * coefficients cc. Zero order coefficient is cc(1) = cc[0]
 */
static double
poly(const double *cc, int nord, double x)
{
	double p, ret_val;

	ret_val = cc[0];
	if (nord > 1) {
		p = x * cc[nord-1];
		for (int j = nord - 2; j > 0; j--)
			p = (p + cc[j]) * x;
		ret_val += p;
	}
	return ret_val;
}


/*
 * Calculates the Shapiro-Wilk W test and its significance level.
 * ALGORITHM AS R94 APPL. STATIST. (1995) vol.44, no.4, 547-551.
 *
 * Preconditions: The input x must be sorted.
 */
int
swilk(const double *x, int n, double *w, double *pw)
{
	int nn2 = n / 2;
	double a[nn2 + 1]; /* 1-based */


	double small = 1e-19;

	/* polynomial coefficients */
	double g[2] = { -2.273,.459 };
	double c1[6] = { 0.,.221157,-.147981,-2.07119, 4.434685, -2.706056 };
	double c2[6] = { 0.,.042981,-.293762,-1.752461,5.682633, -3.582633 };
	double c3[4] = { .544,-.39978,.025054,-6.714e-4 };
	double c4[4] = { 1.3822,-.77857,.062767,-.0020322 };
	double c5[4] = { -1.5861,-.31082,-.083751,.0038915 };
	double c6[3] = { -.4803,-.082676,.0030302 };

	/* Local variables */
	int i, j, i1;

	double ssassx, summ2, ssumm2, gamma, range;
	double a1, a2, an, m, s, sa, xi, sx, xx, y, w1;
	double fac, asa, an25, ssa, sax, rsn, ssx, xsx;

	*pw = 1.;
	if (n < 3) {
		errno = EDOM;
		err("too few samples");
		return -1;
	}

	/* Don't even attempt to provide a P-value with large samples */
	if (n > 5000) {
		errno = EDOM;
		err("too many samples");
		return -1;
	}

	an = (double) n;

	if (n == 3) {
		a[1] = sqrt(0.5);
	} else {
		an25 = an + .25;
		summ2 = 0.;
		for (i = 1; i <= nn2; i++) {
			a[i] = qnorm((i - 0.375) / an25, 0., 1., 1, 0);
			double r__1 = a[i];
			summ2 += r__1 * r__1;
		}
		summ2 *= 2.;
		ssumm2 = sqrt(summ2);
		rsn = 1. / sqrt(an);
		a1 = poly(c1, 6, rsn) - a[1] / ssumm2;

		/* Normalize a[] */
		if (n > 5) {
			i1 = 3;
			a2 = -a[2] / ssumm2 + poly(c2, 6, rsn);
			fac = sqrt((summ2 - 2. * (a[1] * a[1]) - 2. * (a[2] * a[2]))
					/ (1. - 2. * (a1 * a1) - 2. * (a2 * a2)));
			a[2] = a2;
		} else {
			i1 = 2;
			fac = sqrt((summ2 - 2. * (a[1] * a[1])) /
					( 1.  - 2. * (a1 * a1)));
		}
		a[1] = a1;
		for (i = i1; i <= nn2; i++)
			a[i] /= - fac;
	}

	/* Check for zero range */
	range = x[n - 1] - x[0];
	if (range < small) {
		errno = ERANGE;
		return -1;
	}

	/* Check for correct sort order on range - scaled X */
	xx = x[0] / range;
	sx = xx;
	sa = -a[1];
	for (i = 1, j = n - 1; i < n; j--) {
		xi = x[i] / range;
		if (xx - xi > small) {
			/* Fortran had:	 print *, "ANYTHING"
			 * but do NOT; it *does* happen with sorted x (on Intel GNU/linux 32bit):
			 *	shapiro.test(c(-1.7, -1,-1,-.73,-.61,-.5,-.24, .45,.62,.81,1))
			 */
			errno = EINVAL;
			err("not sorted");
			return -1;
		}
		sx += xi;
		i++;
		if (i != j)
			sa += copysign(1.0, i - j) * a[MIN(i, j)];
		xx = xi;
	}

	/* Calculate W statistic as squared correlation between data and
	 * coefficients */
	sa /= n;
	sx /= n;
	ssa = ssx = sax = 0.;
	for (i = 0, j = n - 1; i < n; i++, j--) {
		if (i != j)
			asa = copysign(1.0, i - j) * a[1 + MIN(i, j)] - sa;
		else
			asa = -sa;

		xsx = x[i] / range - sx;
		ssa += asa * asa;
		ssx += xsx * xsx;
		sax += asa * xsx;
	}

	/* W1 equals (1-W) calculated to avoid excessive rounding error for W
	 * very near 1 (a potential problem in very large samples) */
	ssassx = sqrt(ssa * ssx);
	w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
	*w = 1. - w1;

	/* Calculate significance level for W */
	if (n == 3) {
		/* exact P value : */
		double pi6 = 1.90985931710274, /* = 6/pi */
		       stqr = 1.04719755119660; /* = asin(sqrt(3/4)) */
		*pw = pi6 * (asin(sqrt(*w)) - stqr);
		if (*pw < 0.)
			*pw = 0.;
		return 0;
	}

	y = log(w1);
	xx = log(an);
	if (n <= 11) {
		gamma = poly(g, 2, an);
		if (y >= gamma) {
			*pw = 1e-99;/* an "obvious" value, was 'small' which was 1e-19f */
			return 0;
		}
		y = -log(gamma - y);
		m = poly(c3, 4, an);
		s = exp(poly(c4, 4, an));
	} else {
		/* n >= 12 */
		m = poly(c5, 4, xx);
		s = exp(poly(c6, 3, xx));
	}
	/*DBG printf("c(w1=%g, w=%g, y=%g, m=%g, s=%g)\n",w1,*w,y,m,s); */

	*pw = pnorm(y, m, s, 0/* upper tail */, 0);

	return 0;
}
