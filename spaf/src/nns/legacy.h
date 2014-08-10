#ifndef NNS_LEGACY_H
#define NNS_LEGACY_H

#include "includes.h"
#include "sparse_matrix.h"

typedef struct nns_struct nns_type;

nns_type*
PrepareNNS_Legacy(
const double    coef_min,
const double    coef_max,
const double    coef_round_off,
const int       NbOfNodes,
      double**  GeoTab,
      int**     ConTab,
const SparseMatrix_t* InvConTab );

void
WriteNNS(
const nns_type* nns );

bool
IsInsideBoundingBox_Legacy(
const nns_type* nns,
const double    point[4] );

int
DoNNS_Legacy(
const nns_type*  nns,
const double  point[4] );

#endif
