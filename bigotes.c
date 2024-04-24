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
#include "stats.h"
#include "common.h"

static int read_from_stdin = 0;
static int use_wall_clock = 0;
static int use_exec = 0;
static int use_machine_output = 0;
static int be_quiet = 0;
static int trim_outliers = 0;
static long min_samples = 30;
static const char *output_fname = "bigotes.csv";

struct sampling {
	long nmax;
	long nmin;
	long n;
	double *samples;
	double last;
	double wall;
	char *const *cmd;
	char *const *argv;
	int argc;
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

static void
drain(FILE *f)
{
	char line[4096];
	while (fgets(line, 4096, f) != NULL)
		;
}

static int
process_stdout(FILE *f, double *metric)
{
	if (use_wall_clock) {
		drain(f);
		return 0;
	}

	char line[4096];
	if (fgets(line, 4096, f) == NULL) {
		err("missing stdout line");
		return -1;
	}

	char *nl = strchr(line, '\n');
	if (nl != NULL)
		*nl = '\0';

	if (sscanf(line, "%le", metric) != 1) {
		err("failed to read metric from stdout");
		return -1;
	}

	/* Drain the rest of the stdout */
	while (fgets(line, 4096, f) != NULL)
		;

	return 0;
}

static int
do_exec(char *const argv[], double *metric)
{
	int fd[2];
	const int R = 0, W = 1;

	if (pipe(fd) == -1) {
		err("pipe failed:");
		return -1;
	}

	pid_t pid = fork();
	if (pid < 0) {
		err("fork failed:");
		return -1;
	} else if (pid > 0) {
		/* Parent */
		close(fd[W]);
		FILE *f = fdopen(fd[R], "r");
		if (f == NULL) {
			err("fdopen failed:");
			return -1;
		}
		if (process_stdout(f, metric) != 0)
			return -1;
		fclose(f);
	} else {
		/* Child */
		/* Make the stdout go to the pipe */
		if (dup2(fd[W], STDOUT_FILENO) == -1) {
			err("dup2 failed:");
			return -1;
		}
		close(fd[W]);
		close(fd[R]);

		if (execvp(argv[0], argv) != 0) {
			err("execvp failed:");
			return -1;
		}

		/* Not reached */
		err("execvp unexpected return");
		return -1;
	}

	return 0;
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

	if (use_machine_output) {
		printf("%-10s %e\n", "sw_pvalue", p);
		printf("%-10s %d\n", "sw_normal", (p >= 0.05));
	} else {
		const char *msg = (p < 0.05) ? "NOT normal" : "may be normal";
		printf("    Shapiro-Wilk: W=%.2e, p-value=%.2e (%s)\n", W, p, msg);
	}

}

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
	double sem = NAN;
	double rsem = NAN;
	double mad = NAN;
	double q1 = NAN;
	double q3 = NAN;
	double iqr = NAN;

	/* Need at least two samples */
	if (s->n >= 2) {
		/* Sort samples to take the median */
		qsort(s->samples, s->n, sizeof(double), cmp_double);

		double *absdev = safe_calloc(s->n, sizeof(double));

		q1     = s->samples[s->n / 4];
		median = s->samples[s->n / 2];
		q3     = s->samples[(s->n * 3) / 4];

		iqr = q3 - q1;

		double sum = 0.0;
		long ncorr = 0;
		for (long i = 0; i < s->n; i++) {
			double x = s->samples[i];
			if (trim_outliers) {
				if (x < q1 - 3.0 * iqr || x > q3 + iqr * 3.0)
					continue;
			}
			sum += x;
			ncorr++;
		}

		mean = sum / (double) ncorr;
		double sumsqr = 0.0;
		for (long i = 0; i < s->n; i++) {
			double x = s->samples[i];
			double dev = x - mean;
			absdev[i] = fabs(s->samples[i] - median);
			//printf("absdev[%3ld] = %e\n", i, absdev[i]);
			if (x < q1 - 3.0 * iqr || x > q3 + iqr * 3.0) {
				outliers++;
				if (trim_outliers)
					continue;
			}
			sumsqr += dev * dev;
		}
		qsort(absdev, s->n, sizeof(double), cmp_double);
		mad = absdev[s->n / 2];

		var = sumsqr / ncorr;
		stdev = sqrt(var);
		sem = stdev / sqrt(ncorr);
		rsem = 100.0 * sem * 1.96 / mean;
		free(absdev);
	}

	double rmad = 100.0 * mad / median;

	if (!be_quiet) {
		fprintf(stderr, "\rn=%ld t=%.1fs median=%.2e rmad=%.2f%% rsem=%.2f%% far=%ld    ",
				s->n, s->wall, median, rmad, rsem, outliers);
		fflush(stderr);
	}

	s->last_stats = get_time();
}

static int
should_continue(struct sampling *s)
{
	double dt = get_time() - s->last_stats;
	/* Update stats at 10 FPS */
	if (dt > 1.0 / 10.0)
		stats(s);

	if (s->n >= s->nmax)
		return 0;

	if (s->n < s->nmin)
		return 1;

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

	if (!read_from_stdin) {
		FILE *f = fopen(output_fname, "a");
		if (f == NULL) {
			err("fopen failed:");
			exit(1);
		}
		fprintf(f, "%e\n", metric);
		fclose(f);
	}
}

static void
compute_histogram(struct sampling *s, int nbins, long *count, int cut_outliers)
{
	qsort(s->samples, s->n, sizeof(double), cmp_double);

	double qmin, qmax;
	if (cut_outliers) {
		double q1 = s->samples[s->n / 4];
		double q3 = s->samples[(s->n * 3) / 4];
		double iqr = q3 - q1;
		qmin = q1 - 3.0 * iqr;
		qmax = q3 + 3.0 * iqr;
	} else {
		qmin = s->samples[0];
		qmax = s->samples[s->n - 1];
	}

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
plot_histogram(struct sampling *s, int w, int h, int cut_outliers)
{
	long *count = safe_calloc(w, sizeof(long));

	compute_histogram(s, w, count, cut_outliers);

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
		if (i == 0 && cut_outliers) {
			printf("…");
		} else {
			putchar(' ');
		}
		for (int j = 0; j < w; j++) {
			double cell = barlen[j] - i;
			if (cell > 0.0)
				draw_cell(cell);
			else
				putchar(' ');
		}
		if (i == 0 && cut_outliers) {
			printf("…");
		}
		putchar('\n');
	}
}

static void
print_summary(struct sampling *s)
{
	qsort(s->samples, s->n, sizeof(double), cmp_double);

	double mean = stats_mean(s->samples, s->n);
	double stdev = stats_stdev(s->samples, s->n, mean);
	double skewness = stats_skewness(s->samples, s->n, mean, stdev);
	double kurtosis = stats_kurtosis(s->samples, s->n, mean, stdev);

	double xmin = stats_percentile(s->samples, s->n, 0.0);
	double q1 = stats_percentile(s->samples, s->n, 0.25);
	double median = stats_median(s->samples, s->n);
	double q3 = stats_percentile(s->samples, s->n, 0.75);
	double xmax = stats_percentile(s->samples, s->n, 1.0);

	double mad = stats_mad(s->samples, s->n, median);
	long n = s->n;
	long far = stats_outliers(s->samples, s->n, q1, q3, 3.0);

	double sem = stats_sem(stdev, s->n);
	double rel_stdev = stats_percent(stdev, mean);
	double rel_mad = stats_percent(mad, median);
	double rel_sem = stats_percent(sem * 1.96, mean);
	double rel_far = stats_percent(far, s->n);

	if (use_machine_output) {
		printf("%-10s %e\n", "min", xmin);
		printf("%-10s %e\n", "q1", q1);
		printf("%-10s %e\n", "median", median);
		printf("%-10s %e\n", "mean", mean);
		printf("%-10s %e\n", "q3", q3);
		printf("%-10s %e\n", "max", xmax);
		printf("%-10s %ld\n", "samples", n);
		printf("%-10s %e\n", "wall", s->wall);
		printf("%-10s %e\n", "mad", mad);
		printf("%-10s %e\n", "stdev", stdev);
		printf("%-10s %e\n", "skewness", skewness);
		printf("%-10s %e\n", "kurtosis", kurtosis);
		printf("%-10s %ld\n", "far", far);
		printf("%-10s %e\n", "rfar", rel_far);
		printf("%-10s %e\n", "rmad", rel_mad);
		printf("%-10s %e\n", "rstdev", rel_stdev);
		printf("%-10s %e\n", "sem", sem);
		printf("%-10s %e\n", "rsem", rel_sem);
	} else {
		printf("\n");
		printf("%10s %10s %10s %10s %10s %10s\n",
				"MIN", "Q1", "MEDIAN", "MEAN", "Q3", "MAX");
		printf("% 10.3e % 10.3e % 10.3e % 10.3e % 10.3e % 10.3e \n",
				xmin, q1, median, mean, q3, xmax);
		printf("\n");
		printf("%10s %10s %10s %10s %10s %10s\n",
				"N", "WALL", "MAD", "STDEV", "SKEW", "KURTOSIS");
		printf("%10ld %10.1f % 10.3e % 10.3e % 10.3e % 10.3e\n",
				n, s->wall, mad, stdev, skewness, kurtosis);
		printf("\n");
		printf("%10s %10s %10s %10s %10s %10s\n",
				"FAR", "%FAR", "%MAD", "%STDEV", "SEM", "%SEM");
		printf("%10ld %10.2f % 10.2f % 10.2f % 10.3e % 10.2f\n",
				far, rel_far, rel_mad, rel_stdev, sem, rel_sem);
		printf("\n");
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

static void
print_command(struct sampling *s)
{
	if (use_machine_output) {
		const char *runmode = read_from_stdin ? "stdin" : "exec";
		printf("%-10s %s\n", "runmode", runmode);

		if (!read_from_stdin) {
			printf("%-10s ", "command");
			for (char *const* p = s->cmd; *p; p++) {
				printf("%s ", *p);
			}
			printf("\n");
		}
		return;
	}

	if (read_from_stdin) {
		printf("    Read from stdin\n");
	} else {
		printf("    Cmd: ");
		for (char *const* p = s->cmd; *p; p++) {
			printf("%s ", *p);
		}
		printf("\n");
	}
}

static int
do_sample(char * const cmd[], char * const argv[], int argc)
{
	struct sampling s = { 0 };
	s.nmax = 5000;
	s.nmin = min_samples;
	s.min_time = 30.0;
	s.samples = safe_calloc(s.nmax, sizeof(double));
	s.n = 0;
	s.cmd = cmd;
	s.argv = argv;
	s.argc = argc;
	s.t0 = get_time();

	if (!read_from_stdin) {
		FILE *f = fopen(output_fname, "w");
		if (f == NULL) {
			err("fopen failed:");
			exit(1);
		}
		fclose(f);
	}

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
			if (do_exec(cmd, &metric) != 0)
				return 1;
			double t1 = get_time();
			walltime = t1 - t0;
			if (use_wall_clock)
				metric = walltime;
		}

		add_sample(&s, metric, walltime);
	}

	/* Always recompute the stats with all samples */
	if (!read_from_stdin && !be_quiet) {
		/* Clear stat line */
		fprintf(stderr, "\r                                                            \r");
		fflush(stderr);
	}

	print_summary(&s);
	print_command(&s);

	shapiro_wilk_test(&s);

	if (!use_machine_output) {
		printf("\n"); /* Leave one empty before histogram */
		plot_histogram(&s, 64, 4, trim_outliers);
		printf("\n"); /* Leave one empty after histogram */
	}

	free(s.samples);

	return 0;
}

static void
usage(void)
{
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "    %s [options] [--] COMMAND [ARGS...]\n", progname_get());
	fprintf(stderr, "    %s [options] -i\n", progname_get());
	fprintf(stderr, "See the %s(1) manual for more details.\n", progname_get());
	exit(1);
}

int
main(int argc, char *argv[])
{
	progname_set("bigotes");
	int opt;

	while ((opt = getopt(argc, argv, "imn:wo:qhX")) != -1) {
		switch (opt) {
			case 'i':
				read_from_stdin = 1;
				break;
			case 'w':
				use_wall_clock = 1;
				break;
			case 'o':
				output_fname = optarg;
				break;
			case 'q':
				be_quiet = 1;
				break;
			case 'm':
				use_machine_output = 1;
				break;
			case 'n':
				min_samples = atol(optarg);
				break;
			case 'X':
				trim_outliers = 1;
				break;
			default: /* '?' */
				err("unknown option '%c'", opt);
				/* fall-through */
			case 'h':
				usage();
		}
	}

	if (use_exec + read_from_stdin > 1) {
		err("bad usage: only one operation mode allowed");
		usage();
	}

	if (!use_exec && !read_from_stdin)
		use_exec = 1;

	char * const * cmd = NULL;
	if (!read_from_stdin) {
		if (optind >= argc) {
			err("missing command to run");
			usage();
		} else {
			cmd = &argv[optind];
		}
	}

	if (do_sample(cmd, argv, argc) != 0)
		return 1;

	return 0;
}
