/** \file
 Logging tools.

 This file contains logging tools, for recording events in order to provide
 informations about a simulation. Four types of logging events are available,
 error (aborts the simulation if something goes wrong), warning (warn about
 something potentially bad), info (just gives info at pre-processing) and debug
 (gives info at timestepping). Two level of logging are available, by default
 events up to info are recorded, but by using -g flag when calling spaf, debug
 events are recorded as well. Logging events can be send to different streams,
 such as stdout or a file, the selection is achieved at logging initializing.
 These routines are basically wrappers around the dclog ones, which have been
 slightly modified to enable compiling on some platforms, and other stuff (see
 routines descriptions). See: http://sourceforge.net/projects/dclog/
 
 Usage :
 First call SetLog to initialize logging, choosing which stream to send the logging to. Then record logging events with error, warning, info and debug routines. Finally finish logging by calling CloseLog.
*/

#include "logging.h"
#include "dclog.h"


/// dclog log hanlde
static DCLog* _dclog = NULL;
/// switch to check we have one and only one log
static bool _LogInitialized = false;
/// string to append to log file name : time stamp
static char _TimeStamp[] = ".%d_%m_%y.%H-%M-%S";


//==============================================================================
/** Initialize logging.
 
 Create local handle, set log file to have a time stamp as extension, set level
 of logging and finally open log file.
 Use stdout as LogName to print out to console.
 */
//==============================================================================
void
SetLog(
const char*           LogName,  ///< name of the log file
const unsigned short  level )   ///< level of logging
{
  // to make sure we use only one log
  assert( _LogInitialized == false );

  _dclog = NewDCLog();
  assert( _dclog != NULL );
  assert( DCLogSetTimestampFormat( _dclog, _TimeStamp ) );
  
  // levels go from 0 to 5, see definitions in dclog.h
  assert( level <= 5 );
  assert( DCLogSetLevel( _dclog, level ) );
  assert( DCLogOpen( _dclog, LogName, "w" ) );
  
  _LogInitialized = true;
}
//==============================================================================
/** Write informations events to log.
 
 DCLogWrite has been modified in such a way it can be warpped with varying
 arguments, see www.c-faq.com/varargs/handoff.html for explanations.
 */
//==============================================================================
void
info(
const char *fmt, ///< format string
... )            ///< varying arguments
{
  va_list str_args;
  va_start( str_args, fmt );
  DCLogWrite(_dclog, DCLOG_INFO, fmt, str_args);
  va_end( str_args );
}
//==============================================================================
/** Write warning events to log.
 
 DCLogWrite has been modified in such a way it can be warpped with varying
 arguments, see www.c-faq.com/varargs/handoff.html for explanations.
 */
//==============================================================================
void
warning(
const char *fmt, ///< format string
... )            ///< varying arguments
{
  va_list str_args;
  va_start( str_args, fmt );
  DCLogWrite(_dclog, DCLOG_WARNING, fmt, str_args);
  va_end( str_args );
}
//==============================================================================
/** Write error events to log.
 
 DCLogWrite has been modified in such a way it can be warpped with varying
 arguments, see www.c-faq.com/varargs/handoff.html for explanations.
 */
//==============================================================================
void
error(
const char *fmt, ///< format string
... )            ///< varying arguments
{
  va_list str_args;
  va_start( str_args, fmt );
  DCLogWrite(_dclog, DCLOG_ERROR, fmt, str_args);
  va_end( str_args );

  abort();
}
//==============================================================================
/** Write debug events to log.
 
 DCLogWrite has been modified in such a way it can be warpped with varying
 arguments, see www.c-faq.com/varargs/handoff.html for explanations.
 */
//==============================================================================
void
debug(
const char *fmt, ///< format string
... )            ///< varying arguments
{
  va_list str_args;
  va_start( str_args, fmt );
  DCLogWrite(_dclog, DCLOG_DEBUG, fmt, str_args);
  va_end( str_args );
}
//==============================================================================
/** Close a log.
 
 Close log and destroy handle.
 */
//==============================================================================
void
CloseLog( void )
{
  assert( DCLogClose( _dclog ) );
  assert( DestroyDCLog( _dclog ) );
}
