/* Copyright (c) 2021-2024 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: GPL-3.0-or-later */

#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include "swilk.h"
#include "common.h"

static int read_from_stdin = 0;
static int use_wall_clock = 0;

struct sampling {
	long nmax;
	long nmin;
	long n;
	double *samples;
	double rsem;
	double emad;
	double last;
	double wall;
	double min_rsem;
	double min_emad;
	const char *name;
	double t0;
	double min_time;
	double last_stats;
};

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
do_run(const char *cmd, double *ptime)
{
	int ret = 0;
	FILE *p = popen(cmd, "r");

	if (p == NULL) {
		err("popen failed:");
		return -1;
	}

	char line[4096];
	if (!use_wall_clock) {
		if (fgets(line, 4096, p) == NULL) {
			err("missing stdout line");
			ret = -1;
			goto close_pipe;
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
	}

	/* Drain the rest of the stdout */
	while (fgets(line, 4096, p) != NULL) {
		//fprintf(stderr, "%s", line);
	}

close_pipe:
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

/* Test for normality */
static void
shapiro_wilk_test(struct sampling *s)
{
	double W, p;
	if (swilk(s->samples, s->n, &W, &p) != 0) {
		err("swilk failed");
		return;
	}

	const char *msg = (p < 0.05) ? "NOT normal" : "may be normal";
	printf("Shapiro-Wilk: W=%.2e, p-value=%.2e (%s)\n", W, p, msg);

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
//	long m = 200;
//
//	double *r = safe_calloc(n, sizeof(double));
//	double *absdev = safe_calloc(n, sizeof(double));
//	double *mad = safe_calloc(m, sizeof(double));
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
//		//mad[sample] = median;
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
	//double rstdev = NAN;
	double sem = NAN;
	double rsem = NAN;
	double mad = NAN;
	double rsemad = NAN;
	double q1 = NAN;
	double q3 = NAN;
	double iqr = NAN;
	//double pol = NAN;
	double smin = s->samples[s->n - 1];
	double smax = s->samples[s->n - 1];

#define NMAD 10
	static double oldmad[NMAD];

	if (s->n == 1) {
		for (int i = 0; i < NMAD; i++) {
			oldmad[i] = NAN;
		}
	}

	/* Need at least two samples */
	if (s->n >= 2) {
		/* Sort samples to take the median */
		qsort(s->samples, s->n, sizeof(double), cmp_double);

		double *absdev = safe_calloc(s->n, sizeof(double));

		smin   = s->samples[0];
		q1     = s->samples[s->n / 4];
		median = s->samples[s->n / 2];
		q3     = s->samples[(s->n * 3) / 4];
		smax   = s->samples[s->n - 1];

		//qcd = (q3 - q1) / (q3 + q1);
		iqr = q3 - q1;

		double sum = 0.0;
		long ncorr = 0;
		for (long i = 0; i < s->n; i++) {
			double x = s->samples[i];
			if (x < q1 - 3.0 * iqr || x > q3 + iqr * 3.0)
				continue;
			sum += s->samples[i];
			ncorr++;
		}

		mean = sum / (double) ncorr;
		double sumsqr = 0.0;
		for (long i = 0; i < s->n; i++) {
			double x = s->samples[i];
			double dev = x - mean;
			absdev[i] = fabs(s->samples[i] - median);
			//printf("absdev[%3ld] = %e\n", i, absdev[i]);
			if (x < q1 - 3.0 * iqr || x > q3 + iqr * 3.0)
				outliers++;
			else
				sumsqr += dev * dev;
		}
		qsort(absdev, s->n, sizeof(double), cmp_double);
		//mad = absdev[s->n / 2] * 1.4826;
		mad = absdev[s->n / 2];
		//rsemad = mad_bootstrap(s->samples, s->n);
		//pol = (double) outliers * 100.0 / n;

		var = sumsqr / ncorr;
		stdev = sqrt(var);
		//rstdev = 100.0 * stdev / mean;
		sem = stdev / sqrt(ncorr);
		rsem = 100.0 * sem * 1.96 / mean;
		s->rsem = rsem;
		free(absdev);
	}

	double madsum = rsemad;
	for (int i = 0; i < NMAD; i++)
		madsum += oldmad[i];
	double emad = madsum / (NMAD + 1.0);
	s->emad = emad;
	double rmad = 100.0 * mad / median;

	/* Print the header at the beginning only */
	if (read_from_stdin || s->n == 1) {
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
		printf("%5s  %5s"
				"  %8s  %8s  %8s  %8s  %8s"
				"  %8s  %5s  %5s"
				"  %5s\n",
				"RUN", "WALL",
				"MIN", "Q1", "MEDIAN", "Q3", "MAX",
				"MAD", "%MAD", "%SEM",
				"OUTLIERS");
	}
//RUN   WALL       LAST     MEDIAN        AVG         SD   %RSD   %RSEM
// 89  125.5  5.085e-03  5.075e-03  5.303e-03  3.500e-03  66.00   7.611
//RUN    WALL       LAST     MEDIAN        AVG         SD   %RSD  %RSEM
//  34    3.0  5.110e-03  5.097e-03  5.121e-03  1.327e-04   2.59   0.87
	printf(
			"%s%5ld  %5.1f"
			"  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e"
			"  %8.2e  %5.2f  %5.2f"
			"  %8ld ",
			read_from_stdin ? "" : "\r",
			s->n, s->wall,			/* progress */
			smin, q1, median, q3, smax,	/* centrality */
			mad, rmad, s->rsem,		/* rel. dispersion */
			outliers			/* outliers */
		);
	fflush(stdout);

	for (int i = 1; i < NMAD; i++)
		oldmad[i - 1] = oldmad[i];
	oldmad[NMAD - 1] = rsemad;

	s->last_stats = get_time();
}

static int
should_continue(struct sampling *s)
{
	double dt = get_time() - s->last_stats;
	/* Update stats at 30 FPS */
	if (dt > 1.0 / 30.0)
		stats(s);

	if (s->n < s->nmin)
		return 1;

	if (isnan(s->rsem) || s->rsem > s->min_rsem)
		return 1;

//	if (isnan(s->emad) || s->emad > s->min_emad)
//		return 1;

	if (s->wall < s->min_time)
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

static void
compute_histogram(struct sampling *s, int nbins, long *count)
{
	qsort(s->samples, s->n, sizeof(double), cmp_double);

	double qmin = s->samples[0];
	double qmax = s->samples[s->n - 1];

	double binlen = (qmax - qmin) / (double) nbins;

	long j = 0;

	for (int i = 0; i < nbins; i++) {
		double binstart = qmin + (double) (i + 0) * binlen;
		double binend   = qmin + (double) (i + 1) * binlen;
		count[i] = 0;

		for (; j < s->n; j++) {
			double x = s->samples[j];
			if (x < binstart) {
				/* Left out */
				continue;
			}

			if (x > binend) {
				/* Switch bin */
				break;
			}

			/* Fits this bin */
			count[i]++;
		}
	}
}

static void
draw_cell(double q)
{

	char *utf8[] = { "▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"};

	int i = (q * 8.0) - 1;
	if (i > 7)
		i = 7;
	else if (i < 0)
		i = 0;

	printf("%s", utf8[i]);
}

static void
plot_histogram(struct sampling *s, int w, int h)
{
	long *count = safe_calloc(w, sizeof(long));

	compute_histogram(s, w, count);

	long maxcount = 0;
	for (int i = 0; i < w; i++) {
		if (count[i] > maxcount)
			maxcount = count[i];
	}

	double *barlen = safe_calloc(w, sizeof(double));
	for (int i = 0; i < w; i++) {
		double rel = (double) count[i] / (double) maxcount;
		barlen[i] = (double) h * rel;
	}

	for (int i = h-1; i >= 0; i--) {
		putchar(' ');
		for (int j = 0; j < w; j++) {
			double cell = barlen[j] - i;
			if (cell > 0.0)
				draw_cell(cell);
			else
				putchar(' ');
		}
		putchar('\n');
	}
}


/* Return -1 on error, 0 on success */
static int
do_read(FILE *f, double *metric, int *end)
{
	char line[4096];
	if (fgets(line, 4096, f) == NULL) {
		if (feof(f)) {
			*end = 1;
			return 0;
		}

		err("missing stdout line");
		return -1;
	}

	char *nl = strchr(line, '\n');
	if (nl != NULL)
		*nl = '\0';

	if (sscanf(line, "%le", metric) != 1) {
		err("cannot find number in: %s", line);
		return -1;
	}

	return 0;
}

static int
sample(const char *cmd)
{
	struct sampling s = { 0 };
	s.nmax = 100000;
	s.nmin = 30;
	s.min_rsem = 1.0;
	s.min_emad = 1.0;
	s.min_time = 30.0;
	s.samples = safe_calloc(s.nmax, sizeof(double));
	s.n = 0;
	s.name = cmd;
	s.t0 = get_time();

	FILE *f = fopen("data.csv", "w");
	if (f == NULL) {
		err("fopen failed:");
		exit(1);
	}
	fclose(f);

	while (read_from_stdin || should_continue(&s)) {
		double metric;
		double walltime = 0.0;

		if (read_from_stdin) {
			int end = 0;
			if (do_read(stdin, &metric, &end) != 0) {
				err("cannot read sample from stdin");
				return 1;
			}

			if (end)
				break;
		} else {
			double t0 = get_time();
			if (do_run(cmd, &metric) != 0) {
				err("failed to run benchmark");
				return 1;
			}
			double t1 = get_time();
			walltime = t1 - t0;
			if (use_wall_clock)
				metric = walltime;
		}

		add_sample(&s, metric, walltime);
	}

	/* Always recompute the stats with all samples */
	stats(&s);
	printf("\n"); /* Finish stat line */

	printf("\n"); /* Leave one empty before histogram */
	plot_histogram(&s, 70, 8);
	printf("\n"); /* Leave one empty after histogram */

	shapiro_wilk_test(&s);

	free(s.samples);

	return 0;
}

static void
usage(void)
{
	printf("%s command\n", progname_get());
	exit(1);
}

int
main(int argc, char *argv[])
{
	progname_set("bigotes");
	int opt;
	const char *cmd = "stdin";

	while ((opt = getopt(argc, argv, "hiw")) != -1) {
		switch (opt) {
			case 'i':
				read_from_stdin = 1;
				break;
			case 'w':
				use_wall_clock = 1;
				break;
			case 'h':
			default: /* '?' */
				usage();
		}
	}

	if (!read_from_stdin) {
		if (optind >= argc) {
			err("bad usage: missing command");
			usage();
		} else {
			cmd = argv[optind];
		}
	}

	if (sample(cmd) != 0) {
		err("failed to sample the benchmark");
		return 1;
	}

	return 0;
}
