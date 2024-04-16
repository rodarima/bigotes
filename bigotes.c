/* Copyright (c) 2021-2024 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

static char *progname = "bigotes";

struct sampling {
	long nmax;
	long nmin;
	long n;
	double *samples;
	double rsem;
	double last;
	double wall;
	double min_rsem;
	const char *name;
	double t0;
	double min_time;
};

static void
vaerr(const char *prefix, const char *func, const char *errstr, va_list ap)
{
	if (progname != NULL)
		fprintf(stderr, "%s: ", progname);

	if (prefix != NULL)
		fprintf(stderr, "%s: ", prefix);

	if (func != NULL)
		fprintf(stderr, "%s: ", func);

	vfprintf(stderr, errstr, ap);

	int len = strlen(errstr);

	if (len > 0) {
		char last = errstr[len - 1];
		if (last == ':')
			fprintf(stderr, " %s\n", strerror(errno));
		else if (last != '\n' && last != '\r')
			fprintf(stderr, "\n");
	}
}

static void
verr(const char *prefix, const char *func, const char *errstr, ...)
{
	va_list ap;
	va_start(ap, errstr);
	vaerr(prefix, func, errstr, ap);
	va_end(ap);
}

//static void
//vdie(const char *prefix, const char *func, const char *errstr, ...)
//{
//	va_list ap;
//	va_start(ap, errstr);
//	vaerr(prefix, func, errstr, ap);
//	va_end(ap);
//	abort();
//}

#define rerr(...) fprintf(stderr, __VA_ARGS__)
#define err(...)  verr("ERROR", __func__, __VA_ARGS__)
//#define die(...)  vdie("FATAL", __func__, __VA_ARGS__)
#define info(...) verr("INFO", NULL, __VA_ARGS__)
#define finfo(...) verr("INFO", __func__, __VA_ARGS__)
#define warn(...) verr("WARN", NULL, __VA_ARGS__)

#define dbg(...) do { \
	if (unlikely(is_debug_enabled)) verr("DEBUG", __func__, __VA_ARGS__); \
} while (0);

#define likely(x) __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)
#define UNUSED(x) (void)(x)

/* Poison assert */
#pragma GCC poison assert

#define USE_RET __attribute__((warn_unused_result))

#define ARRAYLEN(x) (sizeof(x)/sizeof((x)[0]))

/* Returns the current time in seconds since some point in the past */
static double
get_time(void)
{
	struct timespec tv;
	if(clock_gettime(CLOCK_MONOTONIC, &tv) != 0)
	{
		perror("clock_gettime failed");
		exit(EXIT_FAILURE);
	}

	return (double)(tv.tv_sec) +
		(double)tv.tv_nsec * 1.0e-9;
}

static int
do_run(char *argv[], double *ptime)
{
	int ret = 0;
	FILE *p = popen(argv[0], "r");

	if (p == NULL) {
		err("popen failed:");
		return -1;
	}

	char line[4096];
	if (fgets(line, 4096, p) == NULL) {
		err("missing stdout line");
		ret = -1;
		goto bad_close;
	}

	char *nl = strchr(line, '\n');
	if (nl != NULL)
		*nl = '\0';

	/* Clean status line */
	//fprintf(stderr, "%s\n", line);

	double time;
	sscanf(line, "%le", &time);
	//printf("got %e\n", time);
	*ptime = time;

	/* Drain the rest of the stdout */
	while (fgets(line, 4096, p) != NULL) {
		//fprintf(stderr, "%s", line);
	}

bad_close:
	pclose(p);

	return ret;
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

//static void
//resample(double *values, long n, double *out)
//{
//	for (long i = 0; i < n; i++) {
//		/* FIXME: Not really uniform */
//		out[i] = values[rand() % n];
//		//printf("out[%ld] = %e\n", i, out[i]);
//	}
//}
//
//static double
//mad_bootstrap(double *values, long n)
//{
//	long m = 1000;
//
//	double *r = calloc(n, sizeof(double));
//	if (r == NULL) {
//		perror("calloc failed");
//		exit(1);
//	}
//
//	double *absdev = calloc(n, sizeof(double));
//	if (absdev == NULL) {
//		perror("calloc failed");
//		exit(1);
//	}
//
//	double *mad = calloc(m, sizeof(double));
//	if (mad == NULL) {
//		perror("calloc failed");
//		exit(1);
//	}
//
//	for (long sample = 0; sample < m; sample++) {
//		resample(values, n, r);
//
//		qsort(r, n, sizeof(double), cmp_double);
//		double median = r[n / 2];
//
//		for (long i = 0; i < n; i++) {
//			absdev[i] = fabs(r[i] - median);
//		}
//
//		qsort(absdev, n, sizeof(double), cmp_double);
//		mad[sample] = absdev[n / 2];
//		//printf("mad[%ld] = %e\n", sample, mad[sample]);
//	}
//
//	double sum = 0.0;
//	for (long i = 0; i < m; i++)
//		sum += mad[i];
//
//	double mean = sum / (double) m;
//	double sumsqr = 0.0;
//	for (long i = 0; i < m; i++) {
//		double dev = mad[i] - mean;
//		sumsqr += dev * dev;
//	}
//
//	double var = sumsqr / m;
//	double stdev = sqrt(var);
//	double sem = stdev / sqrt(m);
//	double rsem = 100.0 * sem * 1.96 / mean;
//
//	free(mad);
//	free(absdev);
//	free(r);
//
//	return rsem;
//}

static void
stats(struct sampling *s)
{
	if (s->n < 1)
		return;

	long outliers = 0;
	double last = s->samples[s->n - 1];
	double median = last;
	double mean = last;
	double var = NAN;
	double stdev = NAN;
	double rstdev = NAN;
	double sem = NAN;
	double rsem = NAN;
	double mad = NAN;
	//double mad_se = NAN;
	double q1 = NAN;
	double q3 = NAN;
	double iqr = NAN;
	double pol = NAN;
	double smin = s->samples[s->n - 1];
	double smax = s->samples[s->n - 1];

	/* Need at least two samples */
	if (s->n >= 2) {
		/* Sort samples to take the median */
		qsort(s->samples, s->n, sizeof(double), cmp_double);

		double *absdev = calloc(s->n, sizeof(double));
		if (absdev == NULL) {
			perror("calloc failed");
			exit(1);
		}

		smin   = s->samples[0];
		q1     = s->samples[s->n / 4];
		median = s->samples[s->n / 2];
		q3     = s->samples[(s->n * 3) / 4];
		smax   = s->samples[s->n - 1];

		//qcd = (q3 - q1) / (q3 + q1);
		iqr = q3 - q1;

		double sum = 0.0;
		for (long i = 0; i < s->n; i++)
			sum += s->samples[i];

		double n = s->n;
		mean = sum / n;
		double sumsqr = 0.0;
		for (long i = 0; i < s->n; i++) {
			double x = s->samples[i];
			double dev = x - mean;
			sumsqr += dev * dev;
			absdev[i] = fabs(s->samples[i] - median);
			//printf("absdev[%3ld] = %e\n", i, absdev[i]);
			if (x < q1 - 1.5 * iqr || x > q3 + iqr * 1.5)
				outliers++;
		}
		qsort(absdev, s->n, sizeof(double), cmp_double);
		mad = absdev[s->n / 2] * 1.4826;
		//mad_se = mad_bootstrap(s->samples, s->n);
		pol = (double) outliers * 100.0 / n;

		var = sumsqr / n;
		stdev = sqrt(var);
		rstdev = 100.0 * stdev / mean;
		sem = stdev / sqrt(n);
		rsem = 100.0 * sem * 1.96 / mean;
		s->rsem = rsem;
		free(absdev);
	}

	/* Print the header at the beginning only */
	if (s->n == 1) {
		//printf("# --- bench6 ---\n");
		//printf("# Min %ld runs, max %ld\n", s->nmin, s->nmax);
		//printf("# Cutoff %%RSEM value set to %f\n", s->min_rsem);
		//printf("# RUN    Number of run\n");
		//printf("# LAST   Value of last run\n");
		//printf("# MEDIAN Median of values until now\n");
		//printf("# AVG    Mean of values until now\n");
		//printf("# SD     Standard deviation\n");
		//printf("# %%RSD   Relative standard deviation to the mean\n");
		//printf("# %%RSEM  Relative standard error of the mean\n");
		printf("%4s  %5s"
				"  %8s  %8s  %8s  %8s  %8s"
				"  %8s  %8s  %8s"
				"  %5s  %5s"
				"  %5s\n",
				"RUN", "WALL",
				"MIN", "Q1", "MEDIAN", "Q3", "MAX",
				"MAD", "IQR", "SD",
				"%RSD", "%RSEM",
				"%OUTLIERS");
	}
//RUN   WALL       LAST     MEDIAN        AVG         SD   %RSD   %RSEM
// 89  125.5  5.085e-03  5.075e-03  5.303e-03  3.500e-03  66.00   7.611
//RUN    WALL       LAST     MEDIAN        AVG         SD   %RSD  %RSEM
//  34    3.0  5.110e-03  5.097e-03  5.121e-03  1.327e-04   2.59   0.87
	printf(
			"\r%4ld  %5.1f"
			"  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e"
			"  %8.2e  %8.2e  %8.2e"
			"  %5.2f  %5.2f"
			"  %5.1f ",
			s->n, s->wall,			/* progress */
			smin, q1, median, q3, smax,	/* centrality */
			mad, iqr, stdev,		/* dispersion */
			rstdev, rsem,			/* rel. dispersion */
			pol				/* outliers */
		);
	fflush(stdout);
}

static int
should_continue(struct sampling *s)
{
	stats(s);

	if (s->n < s->nmin)
		return 1;

	if (s->rsem > s->min_rsem)
		return 1;

	double dt = get_time() - s->t0;
	if (dt < s->min_time)
		return 1;

	return 0;
}

static void
add_sample(struct sampling *s, double metric, double walltime)
{
	if (s->n >= s->nmax) {
		err("too many samples");
		exit(1);
	}

	s->samples[s->n] = metric;
	s->n++;
	s->last = metric;
	s->wall += walltime;

	FILE *f = fopen("data.csv", "a");
	if (f == NULL) {
		err("fopen failed:");
		exit(1);
	}
	fprintf(f, "%e\n", metric);
	fclose(f);
}

static int
sample(char *argv[])
{
	struct sampling s = { 0 };
	s.nmax = 4000;
	s.nmin = 50;
	s.min_rsem = 0.5;
	s.min_time = 60.0; /* At least one minute */
	s.samples = calloc(s.nmax, sizeof(double));
	s.n = 0;
	s.name = argv[0];
	s.t0 = get_time();

	FILE *f = fopen("data.csv", "w");
	if (f == NULL) {
		err("fopen failed:");
		exit(1);
	}
	fclose(f);

	while (should_continue(&s)) {
		double t0 = get_time();
		double metric;
		if (do_run(argv, &metric) != 0) {
			err("failed to run benchmark");
			return 1;
		}
		double t1 = get_time();
		double walltime = t1 - t0;

		add_sample(&s, metric, walltime);
	}

	fprintf(stdout, "\n");

	free(s.samples);

	return 0;
}

int
main(int argc, char *argv[])
{
	(void) argc;

	//while (argc && argv && strcmp(argv, "--") != 0) {
	//	argv++;
	//	argc--;
	//}

	if (sample(argv+1) != 0) {
		err("failed to sample the benchmark");
		return 1;
	}

	return 0;
}