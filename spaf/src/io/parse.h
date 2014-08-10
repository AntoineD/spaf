#ifndef PARSE_H
#define PARSE_H

#include "includes.h"

void
TokenSetComment(
const char string );

void
TokenSetSeparator(
const char string[] );

void
TokenSetNameValueSeparator(
const char string[] );

void
TokenSetCheck(
const bool check );

int
TokenCount( void );

int
TokenizeFileLine(
FILE* FileId );

int
TokenizeString(
char* string );

char*
Token(
const int TokenId );

bool
TokenIs(
const char*   string,
const int TokenId );

bool
TokenNameIs(
const char* string );

bool
TokenValueIs(
const char* string );

char*
TokenName( void );

char*
TokenValue(
const int TokenId );

bool
SectionIs(
const char* SectionName );

bool
SectionIsEnd( void );

void
TokenValueBoundsError(
const bool IssueError );

double
Token2double(
const int TokenId,
const char    low,
const double  low_value,
const double  high_value,
const char    high );

int
Token2int(
const int TokenId,
const char    low,
const int     low_value,
const int     high_value,
const char    high );

#endif
