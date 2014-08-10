#ifndef MEMORY_H
#define MEMORY_H

#include "logging.h"
#include "includes.h"

// Chunk of memory allocated
#define MBLOCK 100


/// Macro to expand all types 2d array memory allocation for 1-offset
#define AllocM_proto(type) \
type** AllocM_##type( size_t NbOfRows, size_t NbOfColumns )

AllocM_proto(float);
AllocM_proto(double);
AllocM_proto(int);
AllocM_proto(short);
AllocM_proto(int);
AllocM_proto(long);
AllocM_proto(bool);

// these macros also check for successful allocation
#define AllocMdouble( NbOfRows, NbOfColumns, pointer ) \
do{ \
(pointer) = AllocM_double( (NbOfRows), (NbOfColumns) ); \
assert_error( (pointer) != NULL , "Memory allocation error" ); \
} while(0)

#define AllocMint( NbOfRows, NbOfColumns, pointer ) \
do{ \
(pointer) = AllocM_int( (NbOfRows), (NbOfColumns) ); \
assert_error( (pointer) != NULL , "Memory allocation error" ); \
} while(0)

/// Macro to expand all types 2d array memory freeing
#define FreeM( pointer ) \
do{ \
if ( (pointer) != NULL ) { free( *(pointer) ); free( (pointer) ); } \
} while(0)

/// Macro to expand all types 1d array memory allocation for 1-offset
// using function with void** arguments brings warning of type mismatch
#define AllocV_proto(type) \
type* AllocV_##type( size_t NbOfEntries )

AllocV_proto(float);
AllocV_proto(double);
AllocV_proto(int);
AllocV_proto(short);
AllocV_proto(int);
AllocV_proto(long);
AllocV_proto(bool);

// these macros also check for successful allocation
#define AllocVdouble( NbOfEntries, pointer ) \
do{ \
(pointer) = AllocV_double( (NbOfEntries) ); \
assert_error( (pointer) != NULL , "Memory allocation error" ); \
} while(0)

#define AllocVint( NbOfEntries, pointer ) \
do{ \
(pointer) = AllocV_int( (NbOfEntries) ); \
assert_error( (pointer) != NULL , "Memory allocation error" ); \
} while(0)

#define AllocVbool( NbOfEntries, pointer ) \
do{ \
(pointer) = AllocV_bool( (NbOfEntries) ); \
assert_error( (pointer) != NULL , "Memory allocation error" ); \
} while(0)


double***
AllocVM_double(
int m,
int n,
int z );

void
FreeVM_double(
double*** myMat,
int  n );

#endif
