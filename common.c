/* Copyright (c) 2021-2023 Barcelona Supercomputing Center (BSC)
 * SPDX-License-Identifier: MIT */

#include "common.h"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>

int is_debug_enabled = 0;

static char *progname = NULL;
static int use_prefix = 0;
static int use_func = 0;

void
progname_set(char *name)
{
	progname = name;
}

const char *
progname_get(void)
{
	return progname;
}

void
enable_debug(void)
{
	is_debug_enabled = 1;
}

static void
vaerr(const char *prefix, const char *func, const char *errstr, va_list ap)
{
	if (progname != NULL)
		fprintf(stderr, "%s: ", progname);

	if (use_prefix && prefix != NULL)
		fprintf(stderr, "%s: ", prefix);

	if (use_func && func != NULL)
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

void
verr(const char *prefix, const char *func, const char *errstr, ...)
{
	va_list ap;
	va_start(ap, errstr);
	vaerr(prefix, func, errstr, ap);
	va_end(ap);
}

void
vdie(const char *prefix, const char *func, const char *errstr, ...)
{
	va_list ap;
	va_start(ap, errstr);
	vaerr(prefix, func, errstr, ap);
	va_end(ap);
	abort();
}

void *
safe_calloc(size_t nmemb, size_t size)
{
	void *p = calloc(nmemb, size);
	if (p == NULL)
		die("calloc failed:");

	return p;
}
