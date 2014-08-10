#ifndef STATISTICS_H
#define STATISTICS_H

#include "includes.h"

void
StatAdd(
const char*   name,
const double  sample );

void
StatFlush(
const char* FileName );

#endif
