/** \file
 Manage sparse matrices.
 
 Here we can create and access data of sparse matrices.
 */

#include "sparse_matrix.h"
#include "includes.h"
#include "memory.h"
#include "output.h"
#include "data_io.h"
#include "logging.h"
#include "file_system.h"

//==============================================================================
/** Free a sparse matrix.
 */
//==============================================================================
void
FreeSparseMatrix(
SparseMatrix_t* m ) ///< sparse matrix
{
  free(m->offset); m->offset = NULL;
  free(m->rows); m->rows = NULL;
  free(m->columns); m->columns = NULL;
  free(m->entries); m->entries = NULL;
  free(m->diagonal); m->diagonal = NULL;
  free(m->entries_n); m->entries_n = NULL;
  free(m->permutation); m->permutation = NULL;
  free(m); 
}
//=============================================================================
/// Finish set up.
//=============================================================================
static void
ProcessSparseMatrix(
SparseMatrix_t* m ) ///< sparse matrix
{
  // comuting the number of entries
  if ( m->diagonal )
  {
    m->NbOfEntries = m->NbOfRowsStored;
    m->NbOfColumns = m->NbOfRowsStored;
    info("\t""\tdiagonal matrix with "INT_FMT" entries""\n", m->NbOfEntries);
    return;
  }
  else
  {
    int n = 0;
    for ( int row = 0 ; row < m->NbOfRowsStored ; row++ )
      n += NbOfEntriesInRow(m, row);
    
    m->NbOfEntries = n;    
  }
  
  info("\t""\t"INT_FMT" rows in""\n", m->NbOfRowsStored);
  info("\t""\t"INT_FMT" entries""\n", m->NbOfEntries);
  
  //----------------------------------------------------------------------------
  // get the first and last column indexes
  // this is to know which nodes interval the matrix acts on  
  //----------------------------------------------------------------------------
  if ( m->columns != NULL )
  {
    int max = 0, min = INT_MAX;
    
    for ( int i = 0 ; i < m->NbOfEntries ; i++ )
    {
      max = m->columns[i] > max ? m->columns[i] : max;
      min = m->columns[i] < min ? m->columns[i] : min;
    }
    
    m->FirstColumn = min;
    m->LastColumn = max;
    m->NbOfColumns = m->LastColumn - m->FirstColumn + 1;
    
    info("\t""\t"INT_FMT" columns""\n", m->NbOfColumns);
    info("\t""\t""columns range from "INT_FMT" to "INT_FMT"\n",
         m->FirstColumn, m->LastColumn);
  }
  //----------------------------------------------------------------------------
  // check if rows are successive
  // if it is the case then we don't have to look up for row number from row index, this speeds up matrix vector product
  // note this is for the left most, upper submatrix only !
  // this could be more general but it doesn't matter as this submatrix matters most for matrix vector product
  //----------------------------------------------------------------------------
  bool SuccessiveRows = true;
  
  if ( m->rows != NULL )
  {
    int max = 0, min = INT_MAX;
    
    for ( int i = 0 ; i < m->NbOfRowsStored  ; i++ )
    {
      max = m->rows[i] > max ? m->rows[i] : max;
      min = m->rows[i] < min ? m->rows[i] : min;
      
      if ( i != m->rows[ i ] ) SuccessiveRows = false;
    }
    info("\t""\t""rows range from "INT_FMT" to "INT_FMT"\n", min, max);
  }
  
  if ( SuccessiveRows == true )
  {
    // free rows array for successive storage, as it is not needed
    free(m->rows);
    // NULLify as it is how we know rows are successive
    m->rows = NULL;
    info("\t""\t""rows are successive""\n");    
  }
  else
    info("\t""\t""rows are not successive""\n");
  
  
  if ( m->symmetric ) info("\t""\t""symmetric""\n");
  else info("\t""\t""not symmetric""\n");

  //----------------------------------------------------------------------------
  // run some statistics on number of entries per row
  //----------------------------------------------------------------------------
  double average = 0.;
  int minimum = INT_MAX,
  maximum = 0;
  
  int counter = 0;
  
  for ( int row = 0 ; row < m->NbOfRowsStored ; row++ )
  {    
    int NbOfEntriesPerColumn = NbOfEntriesInRow(m, row);
    minimum = (NbOfEntriesPerColumn < minimum) ? NbOfEntriesPerColumn : minimum;
    maximum = (NbOfEntriesPerColumn > maximum) ? NbOfEntriesPerColumn : maximum;
    average = ( average * counter + NbOfEntriesPerColumn ) / ( counter + 1 );
    counter++;
  }
  
  info("\t""\t""entries per row : min , av , max : "
       INT_FMT" , "DBL_FMT" , "INT_FMT"\n",
       minimum, average, maximum);
}
//=============================================================================
/// Create a sparse matrix.
//=============================================================================
SparseMatrix_t*
CreateSparseMatrix(
const bool  symmetric,      ///< is the matrix symmetric ?
const int   NbOfRowsStored, ///< number of rows
      int*  offset,         ///< offset array
      int*  columns,        ///< columns array
      int*  rows,           ///< rows array
      double* diagonal,     ///< entries
      double* entries )     ///< entries
{
  info("\t""creating sparse matrix :""\n");
  
  SparseMatrix_t *m = (SparseMatrix_t*) calloc(1, sizeof(SparseMatrix_t));
  
  m->NbOfRowsStored = NbOfRowsStored;
  m->symmetric = symmetric;
  m->offset  = offset;
  m->rows    = rows;
  m->columns = columns;
  m->entries = entries;
  m->diagonal = diagonal;
  m->entries_n = NULL;
  m->permutation = NULL;

  ProcessSparseMatrix(m);
  
  return m;
}
//=============================================================================
/** Process and create matrix diagonal.
 
 If the matrix is symmetric and square, allocate and fill a vector with diagonal terms. The entries array is not modified at all and keeps those terms. Also put diagonal term in the first position in each rows.
 This function must be passed all the matrices that have the same sparsity pattern.
 */
//=============================================================================
void
ProcessDiagonal(
const int NbOfMatrices,     ///< number of matrices
      SparseMatrix_t* m[] ) ///< sparse matrix
{
#ifdef USE_SPARSE_MKL
  // if we use the mkl sparse matrix vector product, there is no need to get the diagonal in first position and we may better keep the reordering the way it was computed.
  return;
#endif
  
  info("\t""processing sparse matrix diagonal""\n");
  
  for ( int i = 0 ; i < NbOfMatrices ; i++ )
    assert( m[i]->entries && m[i]->symmetric && m[i]->NbOfRowsStored == m[i]->NbOfColumns );

  // check the matrices all have the same pattern
  for ( int i = 1 ; i < NbOfMatrices ; i++ )
    assert( NbOfRowsStored(m[0]) == NbOfRowsStored(m[i]) );
  
  for ( int row = 0 ; row < NbOfRowsStored(m[0]) ; row++ )
  {
    for ( int i = 1 ; i < NbOfMatrices ; i++ )
      assert( m[0]->offset[row] == m[i]->offset[row] );
    
    for ( int ip = m[0]->offset[row] ; ip < m[0]->offset[row+1] ; ip++ )
    for ( int i = 1 ; i < NbOfMatrices ; i++ )
      assert( m[0]->columns[ip] == m[i]->columns[ip] );
  }  
  
  // set diagonal term in first position in each row, this is required by the matrix vector product, which is faster that way
  assert_error( NbOfMatrices <= 2, "Increase entry array size");
  
  for ( int row = 0 ; row < NbOfRowsStored(m[0]) ; row++ )
  {
    // where is the diagonal term ?
    int offset = EntryOffset(m[0], row, row);
    
    // save what is presently in first position
    int column = m[0]->columns[m[0]->offset[row]];
    double entry[2];
    for ( int i = 0 ; i < NbOfMatrices ; i++ )
      entry[i] = m[i]->entries[m[i]->offset[row]];
    
    // set diagonal term in first position
    m[0]->columns[m[0]->offset[row]] = row;
    for ( int i = 0 ; i < NbOfMatrices ; i++ )
      m[i]->entries[m[0]->offset[row]] = m[i]->entries[offset];
    
    // switch other term
    m[0]->columns[offset] = column;
    for ( int i = 0 ; i < NbOfMatrices ; i++ )
      m[i]->entries[offset] = entry[i];
  }
  
  for ( int i = 0 ; i < NbOfMatrices ; i++ )
  {
    if ( ! m[i]->diagonal ) AllocVdouble(m[i]->NbOfRowsStored, m[i]->diagonal);
    
    for ( int row = 0 ; row < NbOfRowsStored(m[i]) ; row++ )
      m[i]->diagonal[row] = DiagonalEntryReal(m[i], row);
  }
}
//=============================================================================
/// Return the offset of a row.
//=============================================================================
static int
RowOffset(
const SparseMatrix_t* m,  ///< sparse matrix
const int row ) ///< row index
{
  return m->offset[row];
}
//=============================================================================
/// Return the integer entry with index of a row.
//=============================================================================
int
EntryNode(
const SparseMatrix_t* m,  ///< sparse matrix
const int row,     ///< row index
const int index )  ///< entry index
{
  return m->entries_n[ RowOffset(m,row) + index];
}
//=============================================================================
/// Return the double entry with index of a row.
//=============================================================================
double
EntryReal(
const SparseMatrix_t* m,  ///< sparse matrix
const int row,     ///< row index
const int index )  ///< entry index
{
  return m->entries[ RowOffset(m,row) + index ];
}
//=============================================================================
/// Return the indexth double diagonal entry.
//=============================================================================
double
DiagonalEntryReal(
const SparseMatrix_t* m,  ///< sparse matrix
const int index )  ///< index
{
  int offset = EntryOffset(m, index, index);
  return m->entries[ offset ];
}
//=============================================================================
/// Return the column of an entry with index of a row.
//=============================================================================
int
EntryColumn(
const SparseMatrix_t* m,  ///< sparse matrix
const int row,     ///< row index
const int index )  ///< entry index
{
  return m->columns[ RowOffset(m,row) + index ];
}
//==============================================================================
/// Get row in sparse storage form row index.
//==============================================================================
int
Row(
const SparseMatrix_t*  m, ///< sparse storage
const int RowIndex )  ///< row index
{
  if ( RowsAreSuccessive(m) ) return RowIndex;
  else return m->rows[ RowIndex ];
}
//==============================================================================
/// Return the number of row stored.
//==============================================================================
int
NbOfRowsStored(
const SparseMatrix_t* m ) ///< sparse storage
{
  return m->NbOfRowsStored;
}
//==============================================================================
/// Return the number of columns stored.
//==============================================================================
int
NbOfColumns(
const SparseMatrix_t* m ) ///< sparse storage
{
  return m->NbOfColumns;
}
//==============================================================================
/// Return the number of entries.
//==============================================================================
int
NbOfEntries(
const SparseMatrix_t* m ) ///< sparse storage
{
  return m->NbOfEntries;
}
//=============================================================================
/// Return the number of entries in a row.
//=============================================================================
int
NbOfEntriesInRow(
const SparseMatrix_t* m,  ///< sparse matrix
const int row ) ///< row index
{
  return m->offset[row+1] - m->offset[row];
}
//=============================================================================
/// Return last column index.
//=============================================================================
int
LastColumn(
const SparseMatrix_t* m ) ///< sparse matrix
{
  return m->LastColumn;
}
//=============================================================================
/// Return first column index.
//=============================================================================
int
FirstColumn(
const SparseMatrix_t* m ) ///< sparse matrix
{
  return m->FirstColumn;
}
//==============================================================================
/// Return true if rows are successive, false otherwise.
//==============================================================================
bool
RowsAreSuccessive(
const SparseMatrix_t* m)  ///< sparse storage
{
  return m->rows == NULL;
}
//==============================================================================
/// Return true if matrix is symmetric, false otherwise.
//==============================================================================
bool
IsSymmetric(
const SparseMatrix_t* m)  ///< sparse storage
{
  return m->symmetric;
}
//==============================================================================
/// Return true if matrix is diagonal, false otherwise.
//==============================================================================
bool
IsDiagonal(
const SparseMatrix_t* m)  ///< sparse storage
{
  return m->diagonal && !m->entries;
}
//==============================================================================
/// Get entry offset corresponding to a row and column.
//==============================================================================
int
EntryOffset(
const SparseMatrix_t*  m,  ///< sparse storage
const int row,      ///< row
const int column )  ///< column
{
  int rowIndex = 0;
  
  // if rows are successive then index = row
  if ( RowsAreSuccessive(m) )
    rowIndex = row;
  else
  { // otherwise we have to search in storage
    while ( m->rows[ rowIndex ] != row )
      rowIndex++;
    
    assert_error( rowIndex < m->NbOfRowsStored,
                 "row index beyond number of rows");
  }
  
  int
  offset = m->offset[ rowIndex ],
  rowlength = m->offset[ rowIndex + 1 ] - offset;
  
  // counter to check we stay in row
  int counter = 1;
  
  while ( m->columns[ offset ] != column )
  {
    counter++;
    offset++;
    
    // if there is no such entry in the matrix, return INT_MAX
    if ( counter > rowlength ) return INT_MAX;
  }
  
  return offset;
}
