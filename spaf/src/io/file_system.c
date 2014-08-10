/** \file
  Interactions with file system.
 
  This file contains the functions used to manage files and directories, such as:
  - check the extension of a file
  - check that a file exists
  - list all the files in a directory
  - make a directory
  - move to a directory ()
  - get the working directory
  - get the absolute path
*/

#include <dirent.h>
#include <sys/stat.h>
#include <unistd.h>
#include "logging.h"
#include "strings.h"
#include "file_system.h"

//==============================================================================
/** Check the extension of a file.
  Returns true if the extension matches, false otherwise.
*/
//==============================================================================
bool
CheckFileExtension(
const char* extension,  ///< file extension
const char* FileName )  ///< file name
{
  // first check extension
  if ( ! StringEndsWith(extension, FileName) ) return false;
  
  // then check separator
  return FileName[ strlen(FileName) - strlen(extension) - 1 ] == '.';
}
//==============================================================================
/** Check existence of a file.
  Returns true if the file exists, false otherwise.
 */
//==============================================================================
bool
CheckFileExistence(
const char* FileName )  ///< file name
{
  FILE* FileId = fopen( FileName, "r" );
  
  if ( ! FileId ) return false;
    
  fclose( FileId );
  
  return true;
}
//=============================================================================
/** Make a list of files in a directory.

  If the directory does not exists, NULL is returned.
  Entries . and .. are filtered and not included in the list.
  A string with the name of the files separated by a space is returned.
  If the directory is empty, a zero length string is returned.
  You'll have to manage memory by yourself afterwards.
*/
//=============================================================================
char*
ListFilesInDirectory(
const char* DirName ) ///< name of a directory
{
  if ( DirName == NULL ) return NULL;
    
  DIR *directory = opendir(DirName);
  
  if ( directory == NULL ) return NULL;
  
  // entry in the directory
  struct dirent *entry;
  
  // initialize file list to empty
  char* FileList = StringDuplicate("");

  // loop over entries
  while( (entry = readdir(directory)) != NULL )
  {
    // filter . and .. entries
    if ( StringCompare(entry->d_name,".") || 
         StringCompare(entry->d_name,"..") ) continue;
    
    // append file name and separator to list
    StringAppend(&FileList,entry->d_name);
    StringAppend(&FileList," ");
  }
  
  // close directory
  assert_error( closedir(directory) == 0, "Failed to close directory");
  
  return FileList;
}
//==============================================================================
/** Create a directory if it does not exists then go to it.

 Returns true or false whether all went right.
*/
//==============================================================================
bool
MakeAndGoToDirectory(
const char* DirName ) ///< directory name
{
  // check directory existence
  if ( ! CheckFileExistence(DirName) )
  {
    // if directory doesn't exist, create it
    if ( mkdir(DirName,S_IRWXU) != 0 ) return false;
  }

  // go to directory
  char* temp = GoToDirectory(DirName);
  
  bool retval = ( temp != NULL );
  
  free(temp);
  
  return retval;
}
//==============================================================================
/** Go one directory up.

 Change from current working directory to parent one, returns true or false
 whether the change was successful or not.
*/
//==============================================================================
bool
GoToDirectoryUp( void )
{
  if ( chdir("..") != 0 ) return false;
  
  // now working directory is DirName
  return true;
}
//==============================================================================
/** Get the working directory.
 
 Returns a string with the working directory full path, taken from glibc ,manual.
 Don't forget to free returned string when finished !
 */
//==============================================================================
char*
GetWorkingDirectory( void )
{
  size_t size = 100;
  
  while (1)
  {
    char *buffer = (char*) malloc(size);
    if ( getcwd(buffer, size) == buffer )
      return buffer;
    free(buffer);
    size *= 2;
  }
}
//==============================================================================
/** Go to a directory.
 
 This function set the current working directory to the passed directory, it
 also return the previous working directory. If something goes wrong, it
 returns NULL.
 */
//==============================================================================
char*
GoToDirectory(
const char* DirectoryPath ) ///< path to a directory
{
  // save current working directory
  char *WorkingDirectory = GetWorkingDirectory();
  
  // go to directory
  int retval = chdir(DirectoryPath);
  
  if ( retval != 0 )
  {
    free(WorkingDirectory);
    WorkingDirectory = NULL;
  }
  
  return WorkingDirectory;
}
//==============================================================================
/** Go to a directory but don't return a string, just bool.
 
 This function does the same as GoToDirectory, but don't return an array that will have to be freed, it just return a bool.
 */
//==============================================================================
bool
GoBackToDirectory(
const char* DirectoryPath ) ///< path to a directory
{
  // save current working directory
  char *WorkingDirectory = GoToDirectory( DirectoryPath );
  
  if ( WorkingDirectory == NULL )
    return false;
  else
  {
    free(WorkingDirectory);
    return true;
  }
}
//==============================================================================
/** Get the absolute path.

 Returns a string with the absolute path, given a path.
 Don't forget to free returned string when finished !
 */
//==============================================================================
char*
GetAbsolutePath(
char* path )  ///< a path
{
  // check we have a relative path, i.e. that starts with '/'
  // otherwise return the path itself
  if ( path[0] == '/' )
  {
    path = StringDuplicate(path);
    return path;
  }
        
  // get working directory
  char *AbsolutePath = GetWorkingDirectory();

  // append relative path to it
  StringAppend(&AbsolutePath, "/");
  StringAppend(&AbsolutePath, path);

  return AbsolutePath;
}
