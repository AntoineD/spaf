#ifndef TAUCS_PRECONDITIONNER_H
#define TAUCS_PRECONDITIONNER_H

#include "sparse_matrix.h"

SparseMatrix_t*
TaucsPreconditioners(
const double          drop_tolerance,
const SparseMatrix_t* m );

void
TaucsReorder(
char* reordering,
const int NbOfMatrices,
SparseMatrix_t* m_ptr[] );

#endif
