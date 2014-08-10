/** \file
  \brief Functions that process string.
 
  This file contains the functions that are used to process strings such as:
  - check a string is a number
  - compare to strings
  - check the end of a string
  - create a string
  - append a string to another one
  - concatenate 2 string into a new one
  - convert a string to a numeric data type
*/

#include <ctype.h>
#include <string.h>
#include "logging.h"
#include "strings.h"

//=============================================================================
/// Check a string has only digits, return true if it is, false otherwise.
//=============================================================================
bool
StringIsNumber(
const char* string )
{
  if ( string == NULL || strlen(string) == 0 ) return false;
  
  for ( int i = 0 ; i < (int) strlen(string) ; i++ ) 
    if ( isdigit( string[i] ) == 0 ) return false;

  return true;
}
//==============================================================================
/** Compare 2 strings.
 
 Returns true if strings are equal, false otherwise.
 */
//==============================================================================
bool
StringCompare(
const char* string1,  ///< string
const char* string2 ) ///< string
{
  return strcmp( string1, string2 ) == 0;
}
//==============================================================================
/** Check a string ends with a defined substring.
 
  Returns true if it is, false otherwise.
*/
//==============================================================================
bool
StringEndsWith(
const char* end,      ///< sub string
const char* string )  ///< string
{
  // check end
  assert_warn( end != NULL, "no extension provided" );
  
  if ( string == NULL || strlen(string) < strlen(end) ) return false;
  
  return StringCompare( &string[ strlen(string) - strlen(end) ], end );
}
//==============================================================================
/** Create a string from another one.
 
 Allocate memory and copy passed string, a la strdup.
 Returns a pointer to created string.
 You'll have to manage memory by yourself afterwards.
 */
//==============================================================================
char*
StringDuplicate(
const char *string )  ///< string to be duplicated
{
  if ( string == NULL ) return NULL;

  char *newstring = (char*) calloc( strlen(string) + 1, sizeof(char) );

  assert_error( newstring != NULL, "Memory allocation error");
  
  strcpy(newstring, string);
  
  return newstring;
}
//==============================================================================
/** Append a string to another one.
 
 Reallocate memory of initial string and use strcat.
 */
//==============================================================================
void
StringAppend(
      char**  string1,  ///< initial string
const char*   string2 ) ///< string to be appended
{
  // do nothing and return if there is no string to append
  if ( string2 == NULL || strlen(string2) == 0 ) return;
  
  *string1 = (char*) realloc( *string1,
    (strlen(*string1) + strlen(string2) + 1)*sizeof(char) );
  
  assert_error( *string1 != NULL, "Memory allocation error");
  
  strcat(*string1, string2);
}
//==============================================================================
/** Concatenate 2 strings into a new one.
 
 Allocate memory and copy strings into the new one.
 Returns a pointer to reallocated string.
 You'll have to manage memory by yourself afterwards.
 */
//==============================================================================
char*
StringConcat(
const char* string1,  ///< string
const char* string2 ) ///< string
{
  char *string = StringDuplicate(string1);
  StringAppend(&string, string2);
  
  return string;
}

//==============================================================================
//==============================================================================
// Read and check double
//==============================================================================
//==============================================================================

/*==============================================================================
PURPOSE
  Convert a string to a double or integer value, also check if the value is within
  a defined range.
 
INPUT
  char* strings : parsed strings
  char  low : ] stands for >, [ for >=
  double low_value : low limit for the value
  double high_value : high limit for the value
  char high : ] stands for <=, [ for <
  
OUTPUT :
  return : value of converted string
  
Antoine, 10/11/06
==============================================================================*/
double
String2double(
const char*  string,
const char   low,
const double low_value,
const double high_value,
const char   high )
{  
  // check input
  assert( string != NULL );
  assert( low_value <= high_value );
  
  double value = (double) atof(string);
  
  if      ( low == ']' ) assert( value > low_value );
  else if ( low == '[' ) assert( value >= low_value );
  else
    error( "%s : lower bound switch should be [ or ] : %c", __FUNCTION__, low );
  
  if      ( high == ']' ) assert( value <= high_value );
  else if ( high == '[' ) assert( value < high_value );
  else
    error( "%s : higher bound switch is not [ or ] : %c", __FUNCTION__, high );
  
  return value;
}
//==============================================================================
// Read and check integers
//==============================================================================
int
String2int(
const char* string,
const char  low,
const int   low_value,
const int   high_value,
const char  high )
{  
  // check input
  assert( string != NULL );
  assert( low_value <= high_value );

  int value = atoi(string);
  
  if      ( low == ']' ) assert( value > low_value );
  else if ( low == '[' ) assert( value >= low_value );
  else
    error( "%s : lower bound switch should be [ or ] : %c", __FUNCTION__, low );
  
  if      ( high == ']' ) assert( value <= high_value );
  else if ( high == '[' ) assert( value < high_value );
  else
    error( "%s : higher bound switch is not [ or ] : %c", __FUNCTION__, high );
  
  return value;
}
