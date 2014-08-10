/** \file
 Sparse matrix vector product routines.
 
 The functions in this file are used to compute the product with an add or substract operation. For the OpenMP version, the symetric case requires special attention as the transpose part of the product can update the output vector concurrently. This is why one vector per thread has to be used for this part of the operation.
 */

#include "includes.h"
#include "logging.h"
#include "memory.h"
#include "linalg.h"
#include "smvp.h"

#ifdef _OPENMP
//=============================================================================
/// File wise local variables
//=============================================================================
/// number of threads
static int _NbOfThreads = 0;
/// array of pointers of size _NbOfThreads, each thread have its own *_x vector to store the result. The first pointer, _x[0], is an alias of the array that acutally store the result, and thus is not allocated.
static double **_x = NULL;
/// maximum size of each *_x array, because all product do not involve the same vector length, we need arrays big enough for all vectors
static int _max_size = 0,
           _size = 0;

//=============================================================================
/// Allocate required memory for thread vectors.
//=============================================================================
static void
InitOMP(
const int     n,    ///< number of rows
      double* out ) ///< output vector
{
  // get current size
  _size = n;
  
#pragma omp parallel \
default(none) \
shared(_x,_NbOfThreads,_max_size)
  {
#pragma omp master
    {
      _NbOfThreads = omp_get_num_threads();
      
      // allocate _x array if not done already
      if ( _x == NULL )
        _x = (double**) calloc(_NbOfThreads, sizeof(double*));
      
      // allocate *_x arrays if there is more than one thread
      if ( _NbOfThreads > 1 )
      {
        // allocate or reallocate if the *_x array are not big enough
        if ( n > _max_size )
        {
          _max_size = n;
          
          // skip the first array _x[0] as it is just an alias
          for ( int thread = 1 ; thread < _NbOfThreads ; thread++ )
          {
            _x[thread] = (double*) realloc( _x[thread], _max_size * sizeof(double) );
            assert_error( _x[thread], "Allocation error");
          }
        }
      }
    }
  }

  // alias to the actual result vector
  _x[0] = out;
  
  // others are zeroed
  for ( int thread = 1 ; thread < _NbOfThreads ; thread++ )
    scal(_size, 0., _x[thread]);
}
//=============================================================================
/** Gather openmp vectors for matrix vector product.
 
 out = _x[0] = Sum over other threads of _x arrays
 */
//=============================================================================
static void
FinishOMP( void )
{
  for ( int thread = 1 ; thread < _NbOfThreads ; thread++ )
    axpby(_size, 1., _x[thread], 1., _x[0]);
}
#endif



//=============================================================================
/// Generic matrix multiply and add.
//=============================================================================
void
GenAdd(
const int     n,        ///< number of rows
const int*    offset,   ///< matrix offset
const int*    columns,  ///< matrix columns
const int*    rows,     ///< matrix rows
const double* entries,  ///< matrix entries
const double* in,       ///< input vector
      double* out )     ///< output vector
{
#pragma omp parallel for \
default(none) \
shared(offset,columns,in,out,entries,rows)
  for ( int row_id = 0 ; row_id < n ; row_id++ )
  {
    register int row = rows[ row_id ];    
    register double out_row = out[row];
    for ( int col = offset[row_id] ; col < offset[row_id+1] ; col++ )
      out_row += entries[col] * in[ columns[col] ];
    out[row] = out_row;
  }
}
//=============================================================================
/// Generic matrix multiply and substract.
//=============================================================================
void
GenSub(
const int     n,        ///< number of rows
const int*    offset,   ///< matrix offset
const int*    columns,  ///< matrix columns
const int*    rows,     ///< matrix rows
const double* entries,  ///< matrix entries
const double* in,       ///< input vector
      double* out )     ///< output vector
{
#pragma omp parallel for \
default(none) \
shared(offset,columns,in,out,entries,rows)
  for ( int row_id = 0 ; row_id < n ; row_id++ )
  {
    register int row = rows[ row_id ];    
    register double out_row = out[row];
    for ( int col = offset[row_id] ; col < offset[row_id+1] ; col++ )
      out_row -= entries[col] * in[ columns[col] ];
    out[row] = out_row;
  }
}
//=============================================================================
/// Generic matrix with successive rows multiply and add.
//==============================================================================
void
GenSuccAdd(
const int     n,        ///< number of rows
const int*    offset,   ///< matrix offset
const int*    columns,  ///< matrix columns
const double* entries,  ///< matrix entries
const double* in,       ///< input vector
      double* out )     ///< output vector
{
#pragma omp parallel for \
default(none) \
shared(offset,columns,in,out,entries)
  for ( int row = 0 ; row < n ; row++ )
  {
    register double out_row = out[row];
    for ( int col = offset[row] ; col < offset[row+1] ; col++ )
      out_row += entries[col] * in[ columns[col] ];
    out[row] = out_row;
  }
}
//=============================================================================
/** Generic matrix with successive rows multiply by tranpose and add.
 
 Here the number of column is required only for openmp as the output vector length is equal to the number of columns.
 */
//==============================================================================
void
GenSuccAdd_T(
const int     nc,       ///< number of columns
const int     n,        ///< number of rows
const int*    offset,   ///< matrix offset
const int*    columns,  ///< matrix columns
const double* entries,  ///< matrix entries
const double* in,       ///< input vector
      double* out )     ///< output vector
{
#ifdef _OPENMP
  InitOMP(nc, out);
#endif
  
#pragma omp parallel \
default(none) \
shared(offset,columns,in,_x,entries)
  {
#ifdef _OPENMP
    double *x = _x[omp_get_thread_num()];
#else
    double *x = out;
#endif
    
#pragma omp for
    for ( int row = 0 ; row < n ; row++ )
    {
      register double in_row = in[row];
      for ( int col = offset[row] ; col < offset[row+1] ; col++ )
        x[ columns[col] ] += entries[col] * in_row;
    }
  }

#ifdef _OPENMP
  FinishOMP();
#endif
}
//=============================================================================
/// Symmetric matrix with successive rows and diagonal first multiply and add.
//==============================================================================
void
SymDiagFirstSuccAdd(
const int     n,        ///< number of rows
const int*    offset,   ///< matrix offset
const int*    columns,  ///< matrix columns
const double* entries,  ///< matrix entries
const double* diagonal, ///< matrix entries
const double* in,       ///< input vector
      double* out )     ///< output vector
{  
#pragma omp parallel for \
default(none) \
shared(out,in,diagonal)
  for ( int row = 0 ; row < n ; row++ )
    out[row] += diagonal[row] * in[row];
  
#ifdef _OPENMP
  InitOMP(n, out);
#endif
  
#pragma omp parallel \
default(none) \
shared(offset,columns,in,_x,entries,diagonal,out)
  {
#ifdef _OPENMP
    double *x = _x[omp_get_thread_num()];
#else
    double *x = out;
#endif

#pragma omp for
    for ( int row = 0 ; row < n ; row++ )
    {
      register double x_row = x[row];
      register double in_row = in[row];
      for ( int col = offset[row]+1 ; col < offset[row+1] ; col++ )
      {
        register int column = columns[col];
        register double entry = entries[col];
        x_row     += entry * in[column];
        x[column] += entry * in_row;
      }
      x[row] = x_row;
    }
  }    
  
#ifdef _OPENMP
  FinishOMP();
#endif
  
// double *x = out;
//
//  for ( int row = 0 ; row < n ; row++ )
//  {
//    register double x_row = x[row];
//    register double in_row = in[row];
//    for ( int col = offset[row] ; col < offset[row+1] ; col++ )
//    {
//      register int column = columns[col];
//      register double entry = entries[col];
//      x_row     += entry * in[column];
//      if ( row != column ) x[column] += entry * in_row;
//    }
//    x[row] = x_row;
//  }
}
