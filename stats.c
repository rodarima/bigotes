#include "stats.h"

#include <math.h>
#include <stdlib.h>
#include "common.h"

double
stats_mean(double *x, long n)
{
	double sum = 0.0;
	for (long i = 0; i < n; i++)
		sum += x[i];

	return sum / (double) n;
}

double
stats_stdev(double *x, long n, double mean)
{
	double sum = 0.0;
	for (long i = 0; i < n; i++) {
		double d = x[i] - mean;
		sum += d * d;
	}

	return sqrt(sum / (double) (n - 1));
}

double
stats_skewness(double *x, long n, double mean, double stdev)
{
	/* Taken from GSL gsl/statistics/skew_source.c */

	/* find the sum of the cubed deviations, normalized by the sd. */

	/* we use a recurrence relation to stably update a running value so
	   there aren't any large sums that can overflow */

	double skew = 0;
	for (long i = 0; i < n; i++)
	{
		double d = (x[i] - mean) / stdev;
		skew += (d * d * d - skew) / (i + 1);
	}

	return skew;
}

double
stats_kurtosis(double *x, long n, double mean, double stdev)
{
	/* Taken from GSL gsl/statistics/kurtosis_source.c */

	/* find the fourth moment the deviations, normalized by the sd */

	/* we use a recurrence relation to stably update a running value so
	   there aren't any large sums that can overflow */

	double k = 0.0;
	for (long i = 0; i < n; i++) {
		double d = (x[i] - mean) / stdev;
		k += (d * d * d * d - k) / (i + 1);
	}

	return k - 3.0;  /* makes kurtosis zero for a Gaussian */
}

double
stats_median(double *x, long n)
{
	if (n <= 0)
		return NAN;

	if (n % 2 == 0) {
		/* even, compute average of middle values */
		long right = n / 2;
		long left = right - 1;
		return 0.5 * (x[left] + x[right]);
	} else {
		/* odd, just take the middle value */
		return x[n/2];
	}
}

static int
cmp_double(const void *pa, const void *pb)
{
	double a = *(const double *) pa;
	double b = *(const double *) pb;

	if (a < b)
		return -1;
	else if (a > b)
		return 1;
	else
		return 0;
}

double
stats_mad(double *x, long n, double median)
{
	double *absdev = safe_calloc(n, sizeof(double));

	for (long i = 0; i < n; i++)
		absdev[i] = fabs(x[i] - median);

	qsort(absdev, n, sizeof(double), cmp_double);
	double mad0 = stats_median(absdev, n);

	free(absdev);

	return 1.482602218505602 * mad0;
}

double
stats_percentile(double *x, long n, double p)
{
	long i = p * n;
	if (i < 0)
		i = 0;
	else if (i >= n)
		i = n - 1;

	return x[i];
}

long
stats_outliers(double *x, long n, double q1, double q3, double k)
{
	double iqr = q3 - q1;
	double x0 = q1 - k * iqr;
	double x1 = q3 + k * iqr;
	long outliers = 0;
	for (long i = 0; i < n; i++) {
		if (x[i] < x0 || x[i] > x1)
			outliers++;
	}

	return outliers;
}

double
stats_sem(double stdev, long n)
{
	return stdev / sqrt((double) n);
}

double
stats_percent(double dev, double center)
{
	return 100.0 * dev / center;
}
