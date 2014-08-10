/*
 *  logging.h
 *  spaf
 *
 *  Created by Antoine Dechaume on 20/03/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "includes.h"

void
SetLog(
const char*           LogName,
const unsigned short  level );

void
info(
const char *fmt,
     ... );

void
warning(
const char *fmt,
        ... );

void
error(
const char *fmt,
      ... );

void
debug(
const char *fmt,
      ... );

void
CloseLog( void );

#define assert_error(condition,...)  \
((void) ((condition) ? 0 : \
__assert_error(#condition, __FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)))

#define __assert_error(condition, file, function, line, ...) \
((void)warning("%s:%s:%u: failed assertion `%s'\n", file, function, line, condition), \
(void)warning(__VA_ARGS__), \
(void)error("\n"), 0)

// same as assert_error but issue a warning only and keep running

#define assert_warn(condition,...)  \
((void) ((condition) ? 0 : \
__assert_warn(#condition, __FILE__, __FUNCTION__, __LINE__, __VA_ARGS__)))

#define __assert_warn(condition, file, function, line, ...) \
((void)warning("%s:%s:%u: failed assertion `%s'\n", file, function, line, condition), \
(void)warning(__VA_ARGS__), \
(void)warning("\n"), 0)
