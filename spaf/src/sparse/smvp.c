/** \file
 Main entry to sparse matrix vector product.
 
 This file contains only one function that is used to dispatch the sparse matrix vector product according to the type of operation and type of sparse matrix.
 The operations can be add or substract the product to a vector. 
 */

#include "includes.h"
#include "logging.h"
#include "linalg.h"
#include "smvp.h"
#include "memory.h"

/// Use intel mkl sparse routines, was slower when tested
#ifdef USE_SPARSE_MKL
#include "mkl_spblas.h"
#else
#include "smvp_legacy.h"
#endif

//=============================================================================
/** Sparse matrix vector product.
 
 Computes the product of the sparse matrix 'm', if 'transpose' is false or its transpose otherwise, with the vector 'in', and put the result in 'out' if 'operation' is ' ' or add, resp. substract, the result to 'out' if 'operation is '+', resp. '-'.
 */
//=============================================================================
void
SparseMatrixVectorProduct(
      char  operation,
const bool  transpose,
const SparseMatrix_t* m,
      double* in,
      double* out )
{
  assert_error( operation == '+' || operation == '-' || operation == ' ',
    "Bad operation, hould be '+' or '-' or ' '");
  
  double *entries = m->entries,
         *diagonal = m->diagonal;
  
  int n = NbOfRowsStored(m),
      nc = NbOfColumns(m),
      *offset = m->offset,
      *columns = m->columns,
      *rows = m->rows;

#ifdef USE_SPARSE_MKL
  char transa = 'n';
  if ( transpose ) transa = 't';
  
  double alpha = 1.;
  
  char matdescra[4] = {'G','U','N','F'};
  if ( IsSymmetric(m) ) matdescra[0] = 'S';

  double beta = 1.;
  if ( operation == '-' ) beta = -1.;
  
#endif  
  
  // no update of the vector
  if ( operation == ' ' )
  {
    scal(n, 0., out);
    operation = '+';
  }

  
  if ( IsSymmetric(m)    == true &&
       operation         == '+'  &&
       RowsAreSuccessive(m) )
  {
#ifdef USE_SPARSE_MKL
    mkl_dcsrmv(
      &transa, &n, &nc, &alpha, matdescra, entries,
      columns, offset, &offset[1], in, &beta, out);
#else
    SymDiagFirstSuccAdd( n, offset, columns, entries, diagonal, in, out);
#endif
  }
  else if ( operation == '+' )
  {
    if ( transpose == true &&
         RowsAreSuccessive(m) )
    {
#ifdef USE_SPARSE_MKL
      mkl_dcsrmv(
        &transa, &n, &nc, &alpha, matdescra, entries,
        columns, offset, &offset[1], in, &beta, out);
#else
      GenSuccAdd_T( nc, n, offset, columns, entries, in, out);
#endif
    }
    else
    {
      if ( RowsAreSuccessive(m) )
      {
#ifdef USE_SPARSE_MKL
        mkl_dcsrmv(
          &transa, &n, &nc, &alpha, matdescra, entries,
          columns, offset, &offset[1], in, &beta, out);
#else
        GenSuccAdd( n, offset, columns, entries, in, out);
#endif
      }
      else
      {
        GenAdd( n, offset, columns, rows, entries, in, out);
      }
    }
  }
  
  else if ( operation == '-' &&
       transpose == false &&
       RowsAreSuccessive(m) == false )
  {
    GenSub( n, offset, columns, rows, entries, in, out);
  }  
}
