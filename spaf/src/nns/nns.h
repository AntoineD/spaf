#ifndef NNS_H
#define NNS_H

#include "includes.h"
#include "sparse_matrix.h"

// just to give it a name
typedef void nns_t;

void
ReadNNSParameters(
FILE* FileId );

nns_t*
PrepareNNS(
const int       NbOfPoints,
      double**  points,
      int**    connectivity,
const SparseMatrix_t* InvConTab );

bool
IsInsideBoundingBox(
const nns_t*  tree,
const double  point[4] );

int
DoNNS(
const nns_t*  tree,
const double  point[4] );

#endif
