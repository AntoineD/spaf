#ifndef FILE_SYSTEM_H
#define FILE_SYSTEM_H

#include "includes.h"

// macro to check result from a directory move
#define assert_Directory(filename,string) \
assert_error( (filename) != NULL, \
"Failed to go to %s directory",(string) )

bool
CheckFileExtension(
const char* extension,
const char* FileName );

bool
CheckFileExistence(
const char* FileName );

char*
ListFilesInDirectory(
const char* DirName );

bool
MakeAndGoToDirectory(
const char* DirName );

#define assert_MakeAndGoToDirectory(filename) \
assert_error( MakeAndGoToDirectory(filename), \
"Failed to make and go to directory %s", filename )

bool
GoToDirectoryUp( void );

#define assert_GoToDirectoryUp() \
do{ \
if ( GoToDirectoryUp() ) break; \
error("Failed to move one directory up from %s", GetWorkingDirectory() ); \
} while(0)

char*
GoToDirectory(
const char* DirectoryPath );

bool
GoBackToDirectory(
const char* DirectoryPath );

#define assert_GoBackToDirectory(DirectoryPath) \
assert_error( GoBackToDirectory(DirectoryPath), \
"Failed to go back to directory %s", DirectoryPath )

char*
GetWorkingDirectory( void );

char*
GetAbsolutePath(
char* path );

#endif
