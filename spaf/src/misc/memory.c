/** \file
  Memory allocation for vectors and matrices, for different data type.
 */

#include "includes.h"
#include "memory.h"

/// Macro to expand all types 2d array memory allocation for 1-offset
// return NULL if failure
#define AllocM_(type) \
type** AllocM_##type( size_t NbOfRows, size_t NbOfColumns ){ \
NbOfRows++; NbOfColumns++; \
type** array = (type**) calloc( NbOfRows , sizeof(type*) ); \
if ( array == NULL ) return NULL; \
*array = (type*) calloc( NbOfRows * NbOfColumns , sizeof(type) ); \
if ( *array == NULL ) return NULL; \
for( size_t i = 1 ; i < NbOfRows ; i++ ) \
  array[i] = array[0] + i * NbOfColumns; \
return array; }

AllocM_(float)
AllocM_(double)
AllocM_(short)
AllocM_(int)
AllocM_(long)
AllocM_(bool)

/// Macro to expand all types 1d array memory allocation for 1-offset
// return NULL if failure
#define AllocV_(type) \
type* AllocV_##type( size_t NbOfEntries ){ \
NbOfEntries++; \
type* array = (type*) calloc( NbOfEntries , sizeof(type) ); \
return array; }

AllocV_(float)
AllocV_(double)
AllocV_(short)
AllocV_(int)
AllocV_(long)
AllocV_(bool)


/// Allocate a vector of matrices, with 1-offset
double***
AllocVM_double(
int m,
int n,
int z )
{
  m++;
  
  double*** myMatrix = (double***) calloc( m , sizeof(double**) );

  if ( myMatrix == NULL ) return NULL;

  for ( int i = 1 ; i < m ; i++ )
  {
    myMatrix[ i ] = AllocM_double( n, z );
    if ( myMatrix[ i ] == NULL ) return NULL;
  }
  
  return myMatrix;
}
/// Free a vector of matrices
void
FreeVM_double(
double*** myMat,
int  n )
{
  for ( int i = 1 ; i <= n ; i++ ) FreeM( myMat[ i ] );
  free( myMat );
}
