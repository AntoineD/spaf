#ifndef SPARSE_STORAGE_H
#define SPARSE_STORAGE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "includes.h"

typedef struct SparseMatrix SparseMatrix_t;

struct SparseMatrix
{
  bool
  symmetric;      ///< whether the matrix is symmetric or not
  
  int
  NbOfEntries,    ///< number of entries
  NbOfRowsStored, ///< number of rows stored
  NbOfColumns,    ///< number of rows stored
  LastColumn,     ///< last column stored
  FirstColumn,    ///< first column stored
  *offset,        ///< offset list for each row
  *rows,          ///< index list to which row they correspond
  *columns,       ///< index list to which column they correspond
  *permutation,   ///< for reordering a vector
  *invpermutation;///< for reordering a vector

  double
  *entries,   ///< entries values if real
  *diagonal;  ///< diagonal entries
  
  int *entries_n;   ///< entries values if integer
};
  
void
FreeSparseMatrix(
SparseMatrix_t* m );
  
SparseMatrix_t*
CreateSparseMatrix(
const bool  symmetric,  
const int   NbOfRowsStored,
      int*  offset,
      int*  rows,
      int*  columns,
      double* diagonal,
      double* entries );

void
ProcessDiagonal(
const int NbOfMatrices,
      SparseMatrix_t* m[] );

int
Row(
const SparseMatrix_t*  storage,  
const int              RowIndex );

int
EntryOffset(
const SparseMatrix_t*  storage,
const int              row,    
const int              column );

int
NbOfEntriesInRow(
const SparseMatrix_t* m,
const int            row );      

int
NbOfColumns(
const SparseMatrix_t* m );

int
EntryNode(
const SparseMatrix_t* m,
const int             row,     
const int             index );      

double
EntryReal(
const SparseMatrix_t* m,
const int             row,     
const int             index );

double
DiagonalEntryReal(
const SparseMatrix_t* m,
const int index );
  
int
EntryColumn(
const SparseMatrix_t* m,
const int row,
const int index );

int
NbOfRowsStored(
const SparseMatrix_t* m );

int
NbOfEntries(
const SparseMatrix_t* m );

int
LastColumn(
const SparseMatrix_t* m );

int
FirstColumn(
const SparseMatrix_t* m );
  
bool
RowsAreSuccessive(
const SparseMatrix_t* storage);

bool
IsSymmetric(
const SparseMatrix_t* m);

bool
IsDiagonal(
const SparseMatrix_t* m);

#ifdef __cplusplus
}
#endif

#endif
