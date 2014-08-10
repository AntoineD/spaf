/** Reorder sparse matrices, or compute preconditioners with taucs library.
 
 Reordering may speed up the sparse matrix vector product, a preconditioner may reduce the number of iterations of the conjugate gradient solver.
*/

#include "taucs_wrapper.h"
#include "taucs.h"
#include "logging.h"
#include "memory.h"
#include "linalg.h"
#include "strings.h"

//=============================================================================
/** Permute one or more symmetric matrices.
 
 This code comes from taucs, adapted to allow the permutation of more than one symmetric matrix, when the same sparse pattern is shared among multiple matrices.
 */
//=============================================================================
static void
PermuteSymmetrically(
const int*  invperm,        ///< inverse permutation
const int   NbOfMatrices,   ///< number of matrices
      SparseMatrix_t* m[] ) ///< array of matrices pointers
{
  int n = NbOfRowsStored(m[0]),
  *row_len = NULL;
  AllocVint(n-1, row_len);
  
  for ( int row = 0 ; row < n ; row++ )
  for ( int ip = m[0]->offset[row] ; ip < m[0]->offset[row+1] ; ip++ )
  {
    int column = m[0]->columns[ip],
    column_perm = invperm[column],
    row_perm = invperm[row];
    
    if ( column_perm < row_perm )
    {
      int T = column_perm; 
      column_perm = row_perm;
      row_perm = T;
    }
    
    row_len[row_perm] ++;
  }
  
  int *offset = NULL;
  AllocVint(n, offset);
  
  offset[0] = 0;
  for ( int j = 1 ; j <= n ; j++ )
    offset[j] = offset[j-1] + row_len[j-1];
  
  for ( int j = 0 ; j < n ; j++ )
    row_len[j] = offset[j];
  
  int *columns = NULL;
  AllocVint(NbOfEntries(m[0])-1, columns);
  
  // more than 2 matrices ?
  assert_error(NbOfMatrices <= 2, "Increase size of entries array");
  double *entries[2];
  for ( int i = 0 ; i < NbOfMatrices ; i++ )
    AllocVdouble(NbOfEntries(m[0])-1, entries[i]);
  
  for ( int row = 0 ; row < n ; row++ )
  for ( int ip = m[0]->offset[row] ; ip < m[0]->offset[row+1] ; ip++ )
  {
    int column = m[0]->columns[ip],
    column_perm = invperm[column],
    row_perm = invperm[row];
    
    if ( column_perm < row_perm )
    {
      int T = column_perm; 
      column_perm = row_perm;
      row_perm = T;
    }
    
    columns[ row_len[row_perm] ] = column_perm;
    
    for ( int i = 0 ; i < NbOfMatrices ; i++ )
      entries[i][ row_len[row_perm] ] = m[i]->entries[ip];
    
    row_len[row_perm] ++;
  }
  
  // free memory
  free(row_len);
  free(m[0]->offset);
  free(m[0]->columns);
  for ( int i = 0 ; i < NbOfMatrices ; i++ )
  {
    free(m[i]->entries);
    m[i]->entries = entries[i];
    m[i]->offset = offset;
    m[i]->columns = columns;
  }
}
//=============================================================================
/** Reorder one or more sparse matrices.
 
 The type of reordering is : identity, genmmd, md, mmd, amd.
 See taucs manual for more informations. Contrary to taucs, 'indentity' yields no permutation matrix, thus avoid to apply the permutation when required.
 */
//=============================================================================
void
TaucsReorder(
      char* reordering,     ///< type of reordering
const int   NbOfMatrices,   ///< number of matrices
      SparseMatrix_t* m[] ) ///< array of matrices pointers
{
  info("\t""reordering sparse matrix with %s""\n", reordering);
  
  // do nothing for identity reordering
  if ( StringCompare(reordering, "identity") ) return;
  
  // check the matrices have all the same sparse pattern
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
  
  // convert spaf matrix to taucs, shallow copy of attributes
  taucs_ccs_matrix* t = calloc(1, sizeof(taucs_ccs_matrix));
  t->m = NbOfRowsStored(m[0]);
  t->n = t->m;
  t->colptr = m[0]->offset;
  t->rowind = m[0]->columns;
  t->values.d = m[0]->entries;
  t->flags = TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER;
  
  // reorder
  int *perm = NULL, *invperm = NULL;
  
  taucs_ccs_order(t, &perm, &invperm, reordering);
  PermuteSymmetrically(invperm, NbOfMatrices, m);
  
//  free(invperm);

  // set permutation
  for ( int i = 0 ; i < NbOfMatrices ; i++ )
  {
    m[i]->permutation = perm;
    m[i]->invpermutation = invperm;
  }
  
  // shallow freeing
  free(t);
}
//=============================================================================
/** Compute a preconditioner.
 
 Returns the preconditioner as a sparse matrix.
 See taucs manaul for more informations.
 */
//=============================================================================
SparseMatrix_t*
TaucsPreconditioners(
const double          drop_tolerance, ///< drop tolerance
const SparseMatrix_t* m ) ///< sparse matrix
{
  info("\n""Building taucs preconditioner""\n");

  assert_error( IsSymmetric(m), "matrix is not symmetric");
  
  // create taucs matrix from spaf one
  taucs_ccs_matrix* t = calloc(1, sizeof(taucs_ccs_matrix));
  t->m = NbOfRowsStored(m);
  t->n = t->m;
  t->colptr = m->offset;
  t->rowind = m->columns;
  t->values.d = m->entries;
  t->flags = TAUCS_DOUBLE | TAUCS_SYMMETRIC | TAUCS_LOWER;

  // taucs preconditioner
  taucs_ccs_matrix *p = NULL;
  
  #if 0
  // Vaidya
  int modified = 0;
  srand(0);
  int rnd = rand();
  double subgraphs = 1.;
  int stretch_flag = 0;

  prec = taucs_amwb_preconditioner_create(t, rnd, subgraphs, stretch_flag);  
  p = taucs_ccs_factor_llt(prec, drop_tolerance, modified);
  p_f = taucs_ccs_solve_llt;
  #elif 0
  // Vaidya multi-frontal
  srand(0);
  int rnd = rand();
  double subgraphs = 1.;
  int stretch_flag = 0;

  prec = taucs_amwb_preconditioner_create(t, rnd, subgraphs, stretch_flag);  
  void* temp = taucs_ccs_factor_llt_mf(prec);
  p = taucs_supernodal_factor_to_ccs(temp);
  taucs_supernodal_factor_free(temp);
  p_f = taucs_ccs_solve_llt;
  #elif 0
  // supernodal multi-frontal
  void* temp = taucs_ccs_factor_llt_mf(t);
  prec_l = taucs_supernodal_factor_to_ccs(temp);
  taucs_supernodal_factor_free(temp);
  p_f = taucs_supernodal_solve_llt;
  #elif 0
  // supernodal left-looking
  p = taucs_ccs_factor_llt_ll(t);
  p_f = taucs_supernodal_solve_llt;
  
  #elif 1
  // icc
  int modified = 0; // 1 not working
  p = taucs_ccs_factor_llt(t, drop_tolerance, modified);
  #endif

  // shallow freeing
  free(t);

  // convert taucs matrix to spaf matrix
  SparseMatrix_t* spm = CreateSparseMatrix(
    true, p->n, p->colptr, p->rowind, NULL, NULL, p->values.d );
  
  // shallow freeing
  free(p);
  
  return spm;
}
