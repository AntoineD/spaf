#ifndef SOLVER_CG_H
#define SOLVER_CG_H

#include "mesh.h"
#include "fluid.h"

void
ReadConjugateGradientParameters(
FILE* FileId );

int
ConjugateGradient(
const bool* OutsideNodes,
const SparseMatrix_t* A,
const SparseMatrix_t* P,
const double* b,
      double* x );

#endif
