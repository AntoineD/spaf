#ifndef STRINGS_H
#define STRINGS_H

#include "includes.h"

bool
StringIsNumber(
const char* string );

bool
StringCompare(
const char* string1,
const char* string2 );

bool
StringEndsWith(
const char* end,
const char* string );

char*
StringDuplicate(
const char *string );

void
StringAppend(
      char**  string1,
const char*   string2 );

char*
StringConcat(
const char* string1,
const char* string2 );

double
String2double(
const char*  string,
const char   low,
const double low_value,
const double high_value,
const char   high );

int
String2int(
const char* string,
const char  low,
const int   low_value,
const int   high_value,
const char  high );

#endif
