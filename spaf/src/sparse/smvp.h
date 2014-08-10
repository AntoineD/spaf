#ifndef SMVP_H
#define SMVP_H

#include "includes.h"
#include "sparse_matrix.h"

void
SparseMatrixVectorProduct(
      char  operation,
const bool  transpose,
const SparseMatrix_t*  myMat,
      double* in,
      double* out );

#endif
