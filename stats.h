#ifndef STATS_H
#define STATS_H

double stats_mean(double *x, long n);
double stats_stdev(double *x, long n, double mean);
double stats_skewness(double *x, long n, double mean, double stdev);
double stats_kurtosis(double *x, long n, double mean, double stdev);
double stats_median(double *x, long n);
double stats_mad(double *x, long n, double median);
double stats_percentile(double *x, long n, double p);
long stats_outliers(double *x, long n, double q1, double q3, double k);
double stats_sem(double stdev, long n);
double stats_percent(double dev, double center);

#endif /* STATS_H */
