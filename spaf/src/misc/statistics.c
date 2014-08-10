/** \file
  Keep track of statistics and write them to disc.
*/

#include "strings.h"
#include "logging.h"
#include "file_system.h"
#include "statistics.h"

//==============================================================================
/// statistics structure
//==============================================================================
typedef struct
{
  int
    number; ///< number of samples
    
  char
    name[STRING_LEN]; ///< name of sampled values

  double
    average,  ///< samples average
    min,      ///< minimum sample
    max,      ///< maximum sample
    last;     ///< last sample
}
stat_t;

// stat local variables
static stat_t* _stats = NULL;  // array of stat structures
static int _StatNb = 0;  // number of stat structures

//==============================================================================
/// create and allocate the stat structure if necessary, otherwise select stat
//==============================================================================
static stat_t*
StatInit(
const char* name )
{
  // check args
  assert_error( ( name != NULL ) && ( strlen(name) < STRING_LEN ),
           "invalid name string : %s",name);

  // check for stat existence
  // if yes then return the pointer to this structure
  for ( int iStat = 0 ; iStat < _StatNb ; iStat++ )
    if ( StringCompare( name, _stats[iStat].name ) )
      return &_stats[iStat];
  
  // otherwise update stat structure count
  _StatNb++; 

  // create a new element for those samples
  _stats = (stat_t*) realloc( _stats, _StatNb * sizeof(stat_t) );
  
  // set a pointer to the newly created element
  stat_t* stat = &_stats[_StatNb-1];
  
  // initialize this element
  stat->number = 0;
  stat->average = 0.;
  strcpy( stat->name, name );
  stat->max = DBL_MIN;
  stat->min = DBL_MAX;
  stat->last = 0.;

  return stat;
}
//==============================================================================
/// Create and add a sample to a stat
//==============================================================================
void
StatAdd(
const char*   name,   ///< name of stat
const double  sample) ///< sample to add to stat
{
  // select existing stat or create it
  stat_t* stat = StatInit( name );

  // update data structure
  stat->last = sample;
  stat->max = stat->last > stat->max ? stat->last : stat->max;
  stat->min = stat->last < stat->min ? stat->last : stat->min;
  stat->average = ( stat->average * stat->number + stat->last ) /
                  ( stat->number + 1 );

  // update count
  stat->number++;
}
//==============================================================================
/// Write all statistics.
/**
 Write all statistics, to file FileName if FileName is not NULL, otherwise to
 standard output.
*/
//==============================================================================
void
StatFlush(
const char* FileName )
{
  FILE* FilePtr = stdout;
  
  if ( FileName != NULL )
  {
    FilePtr = fopen(FileName, "w");
    assert_error(FilePtr, "Failed to open %s""\n", FileName);
  }
  
  for ( int iStat = 0 ; iStat < _StatNb ; iStat++ )
  {
    // pointer to the current stat element
    stat_t* stat = &_stats[iStat];

    fprintf(FilePtr,"%s :\n",stat->name);
    fprintf(FilePtr,"\t""number  : " INT_FMT "\n",stat->number);
    fprintf(FilePtr,"\t""last    : " DBL_FMT "\n",stat->last); \
    fprintf(FilePtr,"\t""minimum : " DBL_FMT "\n",stat->min); \
    fprintf(FilePtr,"\t""maximum : " DBL_FMT "\n",stat->max); \
    fprintf(FilePtr,"\t""average : " DBL_FMT "\n",stat->average);
    fprintf(FilePtr,"\n");
  }
  
  if ( FileName != NULL ) fclose(FilePtr);
}
