/** \file

  \brief Parse library, used to parse a line of a file.
  
  The main routine is TokenizeFileLine, it read the current line of a file and get tokens
  from this line. A token is a group of characters separated by a space, a tabulation
  or an end of line (this can be changed). Tokens are identified internaly by an
  index which is the position of the token in the parsed line. They can be accessed
  from this index by the different interfaces provided. These interface can provide
  the token value under any type of variable.
*/

#include <assert.h>
#include "parse.h"
#include "strings.h"
#include "logging.h"

/// length of a string line
#define LINE_LEN 500

static char
/// comment keyword
_comment = '#',
/// char to be removed : spaces, tabulations and newline
_separator[] = " \t\n",
/// read line
_NameValueSeparator[] = "=";

static bool _CheckNameValue = true; ///< TokenizeFileLine will check tokens syntax by default

/// number of tokens
#define TOKENS_NB LINE_LEN / 2
static char* _tokens[ TOKENS_NB ]; ///< array of tokens
static int _number = 0; ///< number of tokens in the current parsed line


//==============================================================================
/// set the comment keyword
//==============================================================================
void
TokenSetComment(
const char keyword ) ///< comment keyword
{
  _comment = keyword;
}
//==============================================================================
/// set the tokens separator keywords
//==============================================================================
void
TokenSetSeparator(
const char string[] ) ///< separator keywords
{
  assert_error(strlen(string) != 0, "Separator is empty" );
  strcpy( _separator, string );
}
//==============================================================================
/// set the name-value tokens separator string
//==============================================================================
void
TokenSet_NameValueSeparator(
const char string[] ) ///< name-value separator keywords
{
  strcpy( _NameValueSeparator, string );
}
//==============================================================================
/// set whether to check or not the tokens name-value syntax, on by default
//==============================================================================
void
TokenSetCheck(
const bool check )
{
  _CheckNameValue = check;
}
//==============================================================================
/// Print tokens to stdin
//==============================================================================
#ifdef TOKEN_PRINT
static void
TokenPrint( void )
{
  for ( int counter = 0 ; _tokens[counter] != NULL ; counter++ )
    info( "TokenPrint : token %u = %s\n", counter, _tokens[counter] );
}
#endif
//==============================================================================
/// Returns tokens number
//==============================================================================
int
TokenCount( void )
{
  return _number;
}
//==============================================================================
/** Check tokens name-value format.

  Second token must be the name-value separator, third token must not be empty.
  We look for lines of the form : name = value, 
*/
//==============================================================================
static bool
TokenCheck( void )
{
  if ( *_tokens[1] != *_NameValueSeparator )
  {
    info( "TokenCheck : bad name value separator, %s required\n", _NameValueSeparator );
    return false;
  }
  
  if ( _tokens[2] == NULL )
  {
    info( "TokenCheck : no value !" );
    return false;
  }
  
  return true;
}  
//==============================================================================
/// defined outside TokenizeFileLine, seems to be more stable with some compilers
static char *ptr = NULL,
            _line[LINE_LEN] = "";  ///< string which stores the line read from a file

/** \brief Tokenize the current line of a file.

  Line is read until at least one token is found. A token is any string which
  does not contain a separator char. A separation char is any char contained in
  the sep string. All char after a comment keyword are ignored, the comment
  keyword is defined by the comment char.
  The array of pointers to char, tokens, contains the parsed tokens.
  Returns the number of tokens in the parsed line.
*/
//==============================================================================
int
TokenizeFileLine(
FILE* FileId ) ///< file identifier whose current line is to be tokenized
{
  // add a check for FileId here
  
  _number = 0;  // counter,
  strcpy( _line, "" ); // line,
  ptr = NULL; // and the pointer used to parse line
  for ( int i = 0 ; i < TOKENS_NB ; i++ )
    _tokens[i] = NULL; // initialize tokens array
  
  if ( fgets( _line , LINE_LEN , FileId ) ) // read line
  {
    ptr = strtok(_line, _separator); // get first token
    
    // get tokens until last one or comment char
    while ( ( ptr != NULL ) && ( *ptr != _comment ) )
    {
      _number++; // update tokens count
      _tokens[_number-1] = ptr;
      ptr = strtok(NULL, _separator); // jump to next token
    }

    if ( _number == 0 ) // no token in current line : parse next line
      TokenizeFileLine( FileId );
    else
    {
#ifdef TOKEN_PRINT
      TokenPrint();
#endif
  
      // check the tokens
      if ( _CheckNameValue == true && ( _tokens[2] != NULL ) ) TokenCheck();
    }
  }

  return _number;
}
//==============================================================================
/// Tokenize a string
/**
 string is modified here !
 */
//==============================================================================
int
TokenizeString(
char* string ) ///< string to be tokenized
{
  _number = 0;  // counter,
  ptr = NULL; // and the pointer used to parse line
  for ( int i = 0 ; i < TOKENS_NB ; i++ )
    _tokens[i] = NULL; // initialize tokens array
  
  if ( string != NULL ) // read line
  {
    ptr = strtok(string, _separator); // get first token
    
    // get tokens until last one or comment char
    while ( ( ptr != NULL ) && ( *ptr != _comment ) )
    {
      assert_error( _number < TOKENS_NB, "Increase tokens number");
      _tokens[_number] = ptr;
      _number++; // update tokens count
      ptr = strtok(NULL, _separator); // jump to next token
    }
    
#ifdef TOKEN_PRINT
      TokenPrint();
#endif
      
    // check the tokens
    if ( _CheckNameValue == true && ( _tokens[2] != NULL ) ) TokenCheck();
  }
  
  return _number;
}
//==============================================================================
//==============================================================================
static void
CheckTokenId(
const int TokenId ) ///< token index
{
  assert_error( TokenId <= _number, "token index is greater than tokens number.");
}
//==============================================================================
/// Return a token as a string
//==============================================================================
char*
Token(
const int TokenId ) ///< token index
{
  CheckTokenId( TokenId );
  
  return _tokens[TokenId-1];
}
//==============================================================================
/// Compare a token to a string.
//==============================================================================
bool
TokenIs(
const char*   string,   ///< string to be compared to a token
const int TokenId ) ///< token index
{
  CheckTokenId( TokenId );
  
  if ( _tokens[TokenId-1] == NULL ) return false;
  
  return strncmp( _tokens[TokenId-1] , string ,(int) strlen(string) ) == 0;
}
//==============================================================================
/// Compare the name token of a parsed line to a string
//==============================================================================
bool
TokenNameIs(
const char* string ) ///< string to be compared to the name token
{
  return TokenIs( string , 1 );
}
//==============================================================================
/// Compare the value token of a parsed line to a string
//==============================================================================
bool
TokenValueIs(
const char* string ) ///< string to be compared to the value token
{
  return TokenIs( string , 3 );
}
//==============================================================================
/// Return the name token of a parsed line as a string
//==============================================================================
char*
TokenName( void )
{
  return _tokens[0];
}
//==============================================================================
/// Return the value token of a parsed line as a string
//==============================================================================
char*
TokenValue(
const int TokenId ) ///< token index
{
  CheckTokenId( TokenId );
  
  if ( _tokens[TokenId+1] == NULL )
    return false;

  return _tokens[TokenId+1];
}
//==============================================================================
/// Check the parsed line is a section
//==============================================================================
bool
SectionIs(
const char* SectionName )
{
  return TokenIs( "section"          , 1 ) &&
         TokenIs( _NameValueSeparator, 2 ) &&
         TokenIs( SectionName        , 3 );
}
//==============================================================================
/// Check the parsed line is a section ending
//==============================================================================
bool
SectionIsEnd( void )
{
  return SectionIs( "end" );
}
//******************************************************************************
//******************************************************************************
// Read and check double
//******************************************************************************
//******************************************************************************
/*******************************************************************************
PURPOSE
  Convert a string to a double or integer value, also check if the value is within
  a defined range.
 
INPUT
  char* tokens : parsed tokens
  char  low : ] stands for >, [ for >=
  double low_value : low limit for the value
  double high_value : high limit for the value
  char high : ] stands for <=, [ for <
  
OUTPUT :
  return : value of converted string
  
Antoine, 10/11/06
*******************************************************************************/
//const int TokenId,    ///< token index
//const char    low,        ///< low bound include or exclude switch : '[' or ']'
//const type    low_value,  ///< low bound value
//const type    high_value, ///< high bound value
//const char    high )      ///< high bound include or exclude switch : '[' or ']'

// macro to expand all different cases
#define Token2(type) \
type \
Token2##type( \
const int TokenId, \
const char   low, \
const type   low_value, \
const type   high_value, \
const char   high ) \
{ CheckTokenId( TokenId ); \
  return String2##type( _tokens[TokenId-1], low, low_value, high_value, high ); }

Token2(int)
Token2(double)
