/** \file
 Compute the entries of the operators.
 */

#include "memory.h"
#include "logging.h"
#include "linalg.h"
#include "shape_functions.h"
#include "entries.h"
//==============================================================================
/// Compute local mass matrix.
//==============================================================================
static inline void
LocalMassMatrix(
const int       size,
const double    coef,
const double*   phi,
      double**  local_op )
{
  for ( int row = 1; row <= size; row++ )
  for ( int col = 1; col <= size; col++ )
    local_op[ row ][ col ] += coef * phi[ row ] * phi[ col ];
}
//==============================================================================
/// Compute local stiffness matrix.
//==============================================================================
static inline void
LocalStiffnessMatrix(
const int       size,
const double    coef,
      double**  dphi,
      double**  local_op )
{
  for ( int row = 1; row <= size; row++ )
  for ( int col = 1; col <= size; col++ )
  for ( int x_index = 1; x_index <= 3; x_index++ )
    local_op[ row ][ col ] +=
      coef * dphi[ row ][ x_index ] * dphi[ col ][ x_index ];
}
//==============================================================================
/// Compute local gradient matrix.
//==============================================================================
static inline void
LocalGradientMatrix(
const int       size1,
const int       size2,
const double    coef,
const double*   integrand1,
      double**  integrand2,
      int       dir,
      double**  local_op )
{
  for ( int row = 1 ; row <= size1 ; row++ )
  for ( int col = 1 ; col <= size2 ; col++ )
    local_op[ row ][ col ] += coef * integrand1[ col ] * integrand2[ row ][dir];
}
//==============================================================================
/** Assemble matrix or check there is possible assembly.

  This routine exist because its looping structure is used twice, but with a different intend. It can check if there is an local entry that will contribute to the global one if m is NULL, otherwise it updates the global entry.
 */
//==============================================================================
static bool
Assemble(
const mesh_t* mesh,       ///< mesh structure
const int*    row_map,    ///< row mapping
const int     last_row,   ///< last row to be processed
const int*    col_map,    ///< column mapping
const int     first_col,  ///< first column to process
const int     last_col,   ///< last column to process
const bool    upper_half, ///< process upper half only ?
const int     element,    ///< current element
const int     nb_of_local_rows, ///< number of local rows
const int     nb_of_local_cols, ///< number of local columns
      double**  local_operator, ///< local operator
      SparseMatrix_t* m )       ///< sparse matrix
{
  bool process = false;
  
  for ( int local_row = 1; local_row <= nb_of_local_rows; local_row++ )
  {
    // get global row
    int global_row = mesh->ConTab[ element ][ local_row ];
    
    // skip if node is mapped to nothing, otherwise get mapped node index
    if ( row_map )
    {
      global_row = row_map[ global_row ];
      if ( ! global_row ) continue;
    }
    
    // skip if we are out of the sub-matrix
    if ( global_row > last_row ) continue;
    
    for ( int local_col = 1 ; local_col <= nb_of_local_cols ; local_col++ )
    {
      // get global column
      int global_col = mesh->ConTab[ element ][ local_col ];
      
      // skip if node is mapped to nothing, otherwise get mapped node index
      if ( col_map )
      {
        global_col = col_map[ global_col ];
        if ( ! global_col ) continue;
      }
      
      // skip if we are out of the sub-matrix
      if ( global_col > last_col || global_col < first_col ) continue;
      if ( upper_half && global_row > global_col ) continue;
      
      // if no matrix is provided, just report there is an entry
      if ( ! m )
      {
        process = true;
        break;
      }
      
      // get offset
      int offset = EntryOffset(m, global_row-1, global_col-1);
      
      // update entry value
      m->entries[ offset ] += local_operator[ local_row ][ local_col ];
    }

    // there is an entry
    if ( ! m && process ) break;
  }
  
  return process;
}
//==============================================================================
/** Compute the entries of an operator.

  Only a sub-matrix of the complete operator may be processed, which can be selected with the first and last column, and the last row. If the operator is symmetric, only the upper half part can be processed. Two maps are used to choose between processing velocity or pressure nodes.
 
 On return, the entries array of the sparse matrix m has been allocated and filled.
 */
//==============================================================================
void
GetOperatorEntries(
const mesh_t* mesh,       ///< mesh structure
const int*    row_map,    ///< row mapping
const int     last_row,   ///< last row to be processed
const int*    col_map,    ///< column mapping
const int     first_col,  ///< first column to process
const int     last_col,   ///< last column to process
const bool    upper_half, ///< process upper half only ?
const char    type,       ///< type of operator
const int     dir,        ///< direction for gradient operator
      SparseMatrix_t* m ) ///< sparse matrix
{
  info("\t""computing entries for ");

  // number of local matrix elements in each directions
  int nb_of_local_cols = col_map ? PRE_NODES_PER_EL : NODES_PER_EL,
      nb_of_local_rows = row_map ? PRE_NODES_PER_EL : NODES_PER_EL;
   
  // quadrature data
  gauss_t* gauss = NULL;

  // get quadrature according to operator type
  switch (type)
  {
    case 'm':
      if ( row_map )
        gauss = GaussCreate( nb_of_local_rows, 1 );
      else
        gauss = GaussCreate( nb_of_local_rows, 5 );
      
      info("mass matrix""\n");
      break;
      
    case 's':
      if ( row_map )
        gauss = GaussCreate( nb_of_local_rows, 1 );
      else
        gauss = GaussCreate( nb_of_local_rows, 5 );

      info("stiffness matrix""\n");
      break;
      
    case 'g':
      gauss = GaussCreate( nb_of_local_rows, 5 );
      
      info("gradient matrix in direction %d""\n", dir);
      break;
      
    default:
      error("bad operator type : %c",type);
      break;
  }
  
  AllocVdouble(m->NbOfEntries-1, m->entries);

  // shape functions derivatives at Gauss pts
  double **global_dphi = NULL,
  // inverse of jacobian
  **JacInv = NULL,
  // local operator
  **local_operator = NULL;
  
  AllocMdouble( nb_of_local_rows, 3, global_dphi );
  AllocMdouble( 3, 3, JacInv );
  AllocMdouble( nb_of_local_rows, nb_of_local_cols, local_operator);
  
  // loop over all elements
  for ( int element = 1 ; element <= mesh->NbOfElements ; element++ )
  {
    // check if we have entry to process in this element
    bool process = Assemble(
      mesh, row_map, last_row, col_map, first_col, last_col, upper_half,
      element, nb_of_local_rows, nb_of_local_cols, local_operator, NULL );
    
    // skip if there is no entry to process
    if ( ! process ) continue;
    
    
    // initialize arrays
    scal( (nb_of_local_rows+1)*(nb_of_local_cols+1), 0., local_operator[0]);
    
    // compute local matrix
    for ( int gauss_pt = 1 ; gauss_pt <= gauss->NbOfPoints ; gauss_pt++ )
    {
      double J_times_wt = 0.;
      
      GetJacobianAtPoint(
        nb_of_local_rows, mesh->points, mesh->ConTab[element],
        gauss->local_dphi[gauss_pt], gauss->weight[gauss_pt],
        JacInv, &J_times_wt);
            
      switch (type)
      {
        case 'm':
          LocalMassMatrix(nb_of_local_rows, J_times_wt, gauss->phi[ gauss_pt ],
                          local_operator);
          break;
          
        case 's':
          GetGradientAtPoint(nb_of_local_rows, gauss->local_dphi[gauss_pt],
                             JacInv, global_dphi );
          LocalStiffnessMatrix(nb_of_local_rows, J_times_wt, global_dphi,
                               local_operator);
          break;
          
        case 'g':
          {
            // pressure shape functions
            double psi[PRE_NODES_PER_EL+1] = { 0., 0., 0., 0., 0. };
            GetShapeFunctionsAtPoint( 1, gauss->point[ gauss_pt ], &psi[1] );
            GetGradientAtPoint(nb_of_local_rows, gauss->local_dphi[gauss_pt],
                               JacInv, global_dphi );
            LocalGradientMatrix(nb_of_local_rows, nb_of_local_cols, J_times_wt,
                                psi, global_dphi, dir, local_operator);
          }
          break;

        default:
          error("bad operator type : %c",type);
          break;
      }
    }

    // Assemble matrix
    Assemble(mesh, row_map, last_row, col_map, first_col, last_col, upper_half,
             element, nb_of_local_rows, nb_of_local_cols, local_operator, m );
  }
  
  FreeM( global_dphi );
  FreeM( JacInv );
  FreeM( local_operator );
  FreeGauss(gauss);
}
