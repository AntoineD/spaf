#ifndef SPARSE_STORAGE_CREATE_H
#define SPARSE_STORAGE_CREATE_H

#include "includes.h"
#include "mesh.h"

void
GetSparsePattern(
const mesh_t* mesh,
const int*    row_map,
const int     last_row,
const int*    col_map,
const int     first_col,
const int     last_col,
const bool    upper_half,
      int*    NbOfRows,
      int*    NbOfEntries,
      int**   offset,
      int**   columns,
      int**   rows );

#endif
