#ifndef CHOLMOD_H
#define CHOLMOD_H

#include "sparse_matrix.h"

void
CholmodSolve(
const SparseMatrix_t* m,
const double* b,
      double* x );

#endif
