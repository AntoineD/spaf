#ifndef ENTRIES_H
#define ENTRIES_H

#include "includes.h"
#include "mesh.h"
#include "sparse_matrix.h"

void
GetOperatorEntries(
const mesh_t* mesh,
const int*    row_map,
const int     last_row,
const int*    col_map,
const int     first_col,
const int     last_col,
const bool    upper_half,
const char    type,
const int     dir,
      SparseMatrix_t* m );

#endif
