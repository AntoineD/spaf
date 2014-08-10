/** \file
  Standard and common headers to include, some constants as well.
 */

#ifndef INCLUDES_H
#define INCLUDES_H

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>
#include <limits.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PI 3.14159265358979323846264338327950

#define INT_FMT       "%d"
#define FLT_FMT       "%g"
#define FLT_FMT_FULL  "%11.9le"
#define DBL_FMT       "%g"
#define DBL_FMT_FULL  "%20.18le"

#define STRING_LEN    100

#endif
