/** \file
 Sparse matrix direct solvers.
 
 Diagonal or Choleski factored matrices are supported. 
 */

#include "sparse_matrix.h"
#include "sparse_solve.h"
#include "linalg.h"
#include "logging.h"

// this does not work, do not use or fix it
//#define USE_CHOLMOD
#ifdef USE_CHOLMOD
#include "cholmod_wrapper.h"
#endif

//=============================================================================
/** Direct solver for a Choleski factored matrix.
 */
//==============================================================================
static void
SolveLLt(
const SparseMatrix_t* m,  ///< sparse matrix
const double* b,          ///< right hand side
      double* x )         ///< unknown
{
#ifdef USE_CHOLMOD
  CholmodSolve(m, b, x);
#else
  // get matrix pointers
  int n = NbOfRowsStored(m),
      *offset = m->offset,
      *columns = m->columns;
  double *entries = m->entries;
  
  copy(n, b, x);

#if 0
  // get aux array
  static double *y = NULL;
  static int size = 0;
  
  if ( n > size )
  {
    size = n;
    y = (double*) realloc( y, size*sizeof(double) );
    assert_error( y, "Allocation error");
  }
  
  // Solve L y = b = x
  for ( int j = 0 ; j < n ; j++ )
  {
    int ip = offset[j];
//    assert(columns[ip]==j);
    register double y_j = y[j] = x[j] / entries[ip]; // diagonal
    
    for ( ip++ ; ip < offset[j+1] ; ip++ )
      x[columns[ip]] -= y_j * entries[ip];
  }
  
  // Solve L^T x = y
  for ( int i = n-1 ; i >= 0 ; i-- )
  {
    register double y_i = y[i];
    for ( int jp = offset[i] + 1 ; jp < offset[i+1] ; jp++)
      y_i -= x[columns[jp]] * entries[jp];
    
    x[i] = y_i / entries[offset[i]]; // diagonal
  }
#else
  // Solve L y = b = x
  for ( int j = 0 ; j < n ; j++ )
  {
    int ip = offset[j];
    register double y_j = x[j]  /= entries[ip];
    
    for ( ip++ ; ip < offset[j+1] ; ip++ )
      x[columns[ip]] -= y_j * entries[ip];
  }
  
  // Solve L^T x = y
  for ( int i = n-1 ; i >= 0 ; i-- )
  {
    register double y_i = x[i];
    int jp = offset[i+1] - 1;
    
    for ( ; jp >= offset[i]+1 ; jp--)
      y_i -= x[columns[jp]] * entries[jp];
    
    x[i] = y_i / entries[jp]; // diagonal
  }
#endif
#endif
}
//=============================================================================
/** Entry point to direct solver.
 */
//==============================================================================
void
Solve(
const SparseMatrix_t* m,  ///< sparse matrix
const double* b,          ///< right hand side
      double* x )         ///< unknown
{
  if ( IsDiagonal(m) )
    TensorProduct(NbOfRowsStored(m), m->diagonal, b, x);
  else
    SolveLLt(m, b, x);
}
