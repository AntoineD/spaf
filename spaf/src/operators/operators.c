/** \file
 Operators managment.
 
 Here the mass and stiffness operators for the velocity and pressure, as well as the gradient operator, are computed, stored and used. 
*/

#include "sparse_pattern.h"
#include "entries.h"
#include "strings.h"
#include "logging.h"
#include "memory.h"
#include "output.h"
#include "linalg.h"
#include "data_io.h"
#include "mesh.h"
#include "operators.h"
#include "smvp.h"
#include "statistics.h"
#include "pcg.h"
#include "version_z.h"
#include "parse.h"
#include "taucs_wrapper.h"

/*
 This is how the velocity operators look like. They are symetric, but only submatrix 0 is. We store upper row submatrices only, lower are just transposed. This is how those submatrices should be, most properties are set by inspecting the submatrix, except symmetry which has to be prescribed explicitely.
 
 submatrix 0: symetric, successive row, diagonal first
 submatrix 1: not symetric, no successive row, no diagonal first
 
 |-------|---| 
 |       |   |
 |   0   | 1 |
 |       |   |
 |-------|---| 
 |           |
 |-----------|
 */

/*
 This is how the pressure operators look like. The operators are symetric, but only submatrix 0 and 1 are. 
 
 submatrix 0: symetric, successive row, diagonal first
 submatrix 1: symetric, no successive row, no diagonal first
 
 Only submatrices 0 and 1 are required for solving the pressure poisson equation.
 
 |-------|---| 
 |       |   |
 |   0   | 1 |
 |       |   |
 |-------|   | 
 |   1     \ |
 |-----------| 
 */

// those variables are local to this file
/// sparse matrices
// the following are just pointers that will point to either the _m or _M variables
// for the 1 way version it only points to the _m ones
static SparseMatrix_t
  *v_stiffness = NULL, ///< velocity stiffness matrix
  *v_mass = NULL,      ///< velocity mass matrix
  *gradient = NULL,  ///< velocity-pressure nabla operator
  *p_stiffness = NULL, ///< Pressure stiffness matrix  
/// Preconditioners
  **v_stiffness_prec = NULL,  ///< for the velocity stiffness matrix
  **v_mass_prec = NULL,       ///< for the velocity mass matrix
  **p_stiffness_prec = NULL;  ///< for the pressure stiffness matrix

// operators for the micro grid
static SparseMatrix_t
v_stiffness_m[2], ///< velocity stiffness matrix
v_mass_m[2],      ///< velocity mass matrix
gradient_m[3+1],  ///< velocity-pressure nabla operator
p_stiffness_m[2], ///< Pressure stiffness matrix  
/// Preconditioners
*v_stiffness_prec_m = NULL,  ///< for the velocity stiffness matrix
*v_mass_prec_m = NULL,       ///< for the velocity mass matrix
*p_stiffness_prec_m = NULL;  ///< for the pressure stiffness matrix

// operators for the macro grid
static SparseMatrix_t
v_stiffness_M[2], ///< velocity stiffness matrix
v_mass_M[2],      ///< velocity mass matrix
gradient_M[3+1],  ///< velocity-pressure nabla operator
p_stiffness_M[2], ///< Pressure stiffness matrix  
/// Preconditioners
*v_stiffness_prec_M = NULL,  ///< for the velocity stiffness matrix
*v_mass_prec_M = NULL,       ///< for the velocity mass matrix
*p_stiffness_prec_M = NULL;  ///< for the pressure stiffness matrix

// type of matrix reordering
static char* _reordering = NULL;
// drop tolerance, 0 means use the default preconditioner, > 0 means use incomplete Choleski factorization.
static double _drop_tolerance = 0.;

#ifdef PRE_INCREMENTAL_ROTATIONAL
static SparseMatrix_t
  p_mass[2],            ///< Pressure mass matrix
  *p_mass_prec = NULL;  ///< for the pressure mass matrix
#endif

//=============================================================================
/** Switch current operators to macro grid ones.
 */
//=============================================================================
void
UseMacroGrid(void)
{
  v_stiffness = v_stiffness_M;
  v_mass = v_mass_M;
  gradient = gradient_M;
  p_stiffness = p_stiffness_M;
  v_stiffness_prec = &v_stiffness_prec_M; 
  v_mass_prec = &v_mass_prec_M;
  p_stiffness_prec = &p_stiffness_prec_M;
}
//=============================================================================
/** Switch current operators to micro grid ones.
 */
//=============================================================================
void
UseMicroGrid(void)
{
  v_stiffness = v_stiffness_m;
  v_mass = v_mass_m;
  gradient = gradient_m;
  p_stiffness = p_stiffness_m;
  v_stiffness_prec = &v_stiffness_prec_m; 
  v_mass_prec = &v_mass_prec_m;
  p_stiffness_prec = &p_stiffness_prec_m;
}
//=============================================================================
/** Switch current operators to micro grid ones.
 */
//=============================================================================
bool
UsingMicroGrid(void)
{
  return v_stiffness == v_stiffness_m;
}
//=============================================================================
/** Read preconditioner parameters.
 */
//=============================================================================
void
ReadPreconditionerParameters(
FILE* FileId ) ///< file id
{
  // read reordering
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "reordering" ),
               "Missing preconditioner parameter : ""reordering");
  _reordering = StringDuplicate(Token(3));
  
  // read drop tolerance
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "drop_tolerance" ),
               "Missing preconditioner parameter : ""drop_tolerance");
  _drop_tolerance = Token2double( 3, '[', 0., DBL_MAX,']' );
  
  // check end of section
  TokenizeFileLine( FileId );
  assert_error( SectionIsEnd(),
               "Preconditioner parameters : missing section = end");
  
  // print infos
  info("\nPreconditioner parameters :\n" ); 
  info("\treordering = %s""\n", _reordering);
  
  if (_drop_tolerance == 0.)
    info("\tUsing matrix diagonal.\n" );
  else
    info("\tUsing incomplete Choleski factorization.""\n"
         "\t""drop tolerance = %g\n", _drop_tolerance);  
}
//=============================================================================
/** Apply an operator to a vector and add/substract it to another one.
 
 First make 0 offset pointers to all vectors, then make the matrix vector product and addition or substraction.
 For free nodes submatrix of the velocity mass operator, we have to take care of the permutation.
 */
//=============================================================================
void
ApplyOperator(
const char*   operator, ///< string that select the operation
      double* x_[],     ///< right hand side
      double* rhs_[] )  ///< operand
{
  // 0 offsetting
  double *x[3+1] = {&x_[0][1], &x_[1][1], &x_[2][1], &x_[3][1]},
         *rhs[3+1] = {&rhs_[0][1], &rhs_[1][1], &rhs_[2][1], &rhs_[3][1]};

  //---------------------------------------------------------------------------
  // substract product with free sub matrix of velocity stiffness
  if ( StringCompare(operator,"v_bc") )
  {
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      SparseMatrixVectorProduct( '-', false, &v_stiffness[1], x[dir], rhs[dir] );
  }
  //---------------------------------------------------------------------------
  // multiply with velocity mass and add to rhs
  else if ( StringCompare(operator,"convection") )
  {   
    double *temp = NULL;

    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      if ( v_mass[0].permutation )
      {
        if ( ! temp ) AllocVdouble(NbOfRowsStored(&v_mass[0])-1, temp);
        
        copy(NbOfRowsStored(&v_mass[0]), rhs[dir], temp);
        permute(NbOfRowsStored(&v_mass[0]), v_mass[0].permutation, temp, rhs[dir]);
        permute(NbOfRowsStored(&v_mass[0]), v_mass[0].permutation, x[dir], temp);

        SparseMatrixVectorProduct( '+', false, &v_mass[0], temp, rhs[dir] );

        copy(v_mass[0].NbOfRowsStored, rhs[dir], temp);
        permuteinv(NbOfRowsStored(&v_mass[0]), v_mass[0].permutation, temp, rhs[dir]);
      }
      else
        SparseMatrixVectorProduct( '+', false, &v_mass[0], x[dir], rhs[dir] );
      
      SparseMatrixVectorProduct( '+', false, &v_mass[1], x[dir], rhs[dir] );
    }
    free(temp);
  }
  //---------------------------------------------------------------------------
  // take gradient and add to rhs
  else if ( StringCompare(operator,"gradient") )
  {
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      SparseMatrixVectorProduct( '+', false, &gradient[dir], *x, rhs[dir] );
  }
  //---------------------------------------------------------------------------
  // take divergence and add to rhs
  else if ( StringCompare(operator,"divergence") )
  {
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      SparseMatrixVectorProduct( '+', true, &gradient[dir], x[dir], *rhs );
  }
  //---------------------------------------------------------------------------
  else
    error("Bad solver : %s", operator);
}
//=============================================================================
/** Solve a linear system with an operator.
 
  First make 0 offset pointers to rhs and unknown vectors, then invert system with conjugate gradient and get inversion statistics.
 */
//=============================================================================
void
SolveOperator(
const char*   operator,     ///< string that specify the operation
const bool*   OutsideNodes, ///< list of outside nodes
      double* rhs_[],       ///< right hand side
      double* x_[] )        ///< operand
{
  // 0 offsetting
  double *x[3+1] = {&x_[0][1], &x_[1][1], &x_[2][1], &x_[3][1]},
         *rhs[3+1] = {&rhs_[0][1], &rhs_[1][1], &rhs_[2][1], &rhs_[3][1]};
  
  int NbOfIterationsMax = 0;
  char StatName[30] = "";

  //---------------------------------------------------------------------------
  // velocity mass
  if ( StringCompare(operator,"v_mass") )
  {
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      int NbOfIterations = ConjugateGradient(OutsideNodes, &v_mass[0],
                                             *v_mass_prec, rhs[dir], x[dir]);

      debug( "\tPCG dir "INT_FMT" : %d iterations\n", dir, NbOfIterations );
      
      if ( NbOfIterations > NbOfIterationsMax )
        NbOfIterationsMax = NbOfIterations;
    }
    
    strcpy(StatName, "Velocity projection : PCG");
  }
  //---------------------------------------------------------------------------
  // velocity stiffness
  else if ( StringCompare(operator,"v_stiffness") )
  {
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      int NbOfIterations = ConjugateGradient(OutsideNodes, &v_stiffness[0],
                                             *v_stiffness_prec, rhs[dir], x[dir]);
      
      debug( "\tPCG dir "INT_FMT" : %d iterations\n", dir, NbOfIterations );
      
      if ( NbOfIterations > NbOfIterationsMax )
        NbOfIterationsMax = NbOfIterations;
    }
    
    strcpy(StatName, "Velocity diffusion : PCG");
  }
  //---------------------------------------------------------------------------
  // pressure stiffness
  else if ( StringCompare(operator,"p_stiffness") )
  {
    NbOfIterationsMax = ConjugateGradient(NULL, &p_stiffness[0],
                                          *p_stiffness_prec, *rhs, *x);
    
    debug( "\tPCG : %d iterations\n", NbOfIterationsMax );

    strcpy(StatName, "Pressure prediction : PCG");
  }
  //---------------------------------------------------------------------------
  else
    error("Bad solver : %s", operator);

  // update stats
  StatAdd(StatName, NbOfIterationsMax);
}
//==============================================================================
/** Compute a preconditioner.
 
 if drop tolerance is 0 then the preconditioner is just the inverse of the matrix diagonal. Otherwise the preconditioner is computed with Taucs.
 The preconditioner is returned as a sparse matrix.
 */
//==============================================================================
static SparseMatrix_t*
GetPreconditioner(
const double          drop_tolerance, ///< drop tolerance
const SparseMatrix_t* m )             ///< sparse matrix
{
  if ( drop_tolerance == 0. )
  {
    // use diagonal
    int size = NbOfColumns(m);
    
    double* p = NULL;
    AllocVdouble(size-1, p);
    
    for ( int i = 0 ; i < size ; i++ )
      p[i] = 1. / DiagonalEntryReal(m, i);
    
    return CreateSparseMatrix(false, NbOfRowsStored(m), NULL, NULL, NULL, p,
                              NULL );    
  }
  else
    // use taucs
    return TaucsPreconditioners(drop_tolerance, m);
}
//=============================================================================
/** Compute all operators.
 
 For each operator, compute the sparse pattern, create sparse matrix and get the entries. Then for the operators that are used to solve linear system, get diagonal, reorder and compute preconditioner.
*/
//=============================================================================
void
GetOperators(
const double  tau0,   ///< time discretisation coefficient
const double  Re,     ///< Reynolds number
const mesh_t* mesh )  ///< mesh structure
{
  PrintTitle("Preparing operators");

  // dummy variable used for all operators
  SparseMatrix_t *storage = NULL;

  int first_col = 0,
      last_col = 0,
      last_row = 0,
      nb_of_entries = 0,
      nb_of_rows = 0,
      dir = 0;
  
  bool upper_half = false;
  
  int *rows = NULL,
      *columns = NULL,
      *offset = NULL,
      *row_map = NULL,
      *col_map = NULL;

  //---------------------------------------------------------------------------
  info("\n""Velocity free nodes :""\n");
  
  upper_half = true;
  first_col = 1;
  last_col = mesh->NbOfFreeVelocityNodes;
  last_row = mesh->NbOfFreeVelocityNodes;
  row_map = NULL;
  col_map = NULL;
  dir = 0;

  GetSparsePattern(mesh, row_map, last_row, col_map, first_col, last_col,
                   upper_half, &nb_of_rows, &nb_of_entries, &offset,
                   &columns, &rows);
  
  storage = CreateSparseMatrix(true, nb_of_rows, offset, columns, rows,
                               NULL, NULL);
  
  GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                     upper_half, 's', dir, storage);
  
  v_stiffness[0] = *storage;
    
  GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                     upper_half, 'm', dir, storage);

  v_mass[0] = *storage;
  
  //---------------------------------------------------------------------------
  info("\n""Velocity bc nodes :""\n");
  
  upper_half = true;
  first_col = mesh->NbOfFreeVelocityNodes + 1;
  last_col = mesh->NbOfNodes;
  last_row = mesh->NbOfFreeVelocityNodes;
  row_map = NULL;
  col_map = NULL;
  dir = 0;
  
  GetSparsePattern(mesh, row_map, last_row, col_map, first_col, last_col,
                   upper_half, &nb_of_rows, &nb_of_entries, &offset,
                   &columns, &rows);
  
  storage = CreateSparseMatrix(false, nb_of_rows, offset, columns, rows,
                               NULL, NULL);
  
  GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                     upper_half, 's', dir, storage);
  
  v_stiffness[1] = *storage;

  GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                     upper_half, 'm', dir, storage);
  
  v_mass[1] = *storage;
  
  //---------------------------------------------------------------------------
  info("\n""Pressure free nodes :""\n");
  
  upper_half = true;
  first_col = 1;
  last_col = mesh->NbOfFreePressureNodes;
  last_row = mesh->NbOfFreePressureNodes;
  row_map = mesh->VelToPreNodeMap;
  col_map = mesh->VelToPreNodeMap;
  dir = 0;
  
  GetSparsePattern(mesh, row_map, last_row, col_map, first_col, last_col,
                   upper_half, &nb_of_rows, &nb_of_entries, &offset,
                   &columns, &rows);
  
  storage = CreateSparseMatrix(true, nb_of_rows, offset, columns, rows, NULL,
                               NULL);
  
  GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                     upper_half, 's', dir, storage);

  p_stiffness[0] = *storage;

  //---------------------------------------------------------------------------
  info("\n""Gradient :""\n");

  upper_half = false;
  first_col = 1;
  last_col = mesh->NbOfPressureNodes;
  last_row = mesh->NbOfNodes;
  row_map = NULL;
  col_map = mesh->VelToPreNodeMap;
  
  GetSparsePattern(mesh, row_map, last_row, col_map, first_col, last_col,
                   upper_half, &nb_of_rows, &nb_of_entries, &offset,
                   &columns, &rows);

  storage = CreateSparseMatrix(false, nb_of_rows, offset, columns, rows, NULL,
                               NULL);
    
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                       upper_half, 'g', dir, storage);
    
    gradient[dir] = *storage;
  }

#ifdef PRE_INCREMENTAL_ROTATIONAL
  //---------------------------------------------------------------------------
  info("\n""Pressure free nodes :""\n");
  
  GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                     upper_half, 'm', dir, storage );
  
  p_mass[0] = *storage;
  
  //---------------------------------------------------------------------------
  info("\n""Pressure bc nodes :""\n");
  
  upper_half = true;
  first_col = mesh->NbOfFreePressureNodes + 1;
  last_col = mesh->NbOfFreePressureNodes;
  last_row = mesh->NbOfFreePressureNodes;
  row_map = mesh->VelToPreNodeMap;
  col_map = mesh->VelToPreNodeMap;
  dir = 0;
  
  GetSparsePattern(mesh, row_map, last_row, col_map, first_col, last_col,
                   upper_half, &nb_of_rows, &nb_of_entries, &offset,
                   &columns, &rows);
  
  storage = CreateSparseMatrix(true, nb_of_rows, offset, columns, rows,
                               NULL, NULL );
  
  GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                     upper_half, 's', dir, storage );
  
  p_stiffness[1] = *storage;
  
  GetOperatorEntries(mesh, row_map, last_row, col_map, first_col, last_col,
                     upper_half, 'm', dir, storage);
  
  p_mass[1] = *storage;

  // build pressure correction operator, this will multiply tau*div(vel) by the
  // inverse of ( tau0 * Re ) * Mass matrix in order to get 1 / Re * div(vel)
  for ( int i = 0 ; i < 2 ; i++ )
    scal( NbOfEntries(&p_mass[i]), tau0 * Re, p_mass[i].entries );
  
  // for pressure correction in rotational form we use the lumped modified mass matrix
  // so we store only its inverse
  double *p_mass_prec_[2] = {NULL,NULL};
  
  // first get inverse
  p_mass_prec[0] = GetPreconditioner( &p_mass[0] );
  p_mass_prec[1] = GetPreconditioner( &p_mass[1] );
  
  // then concatenate entries
  AllocVdouble( mesh->NbOfPressureNodes, p_mass_prec );
  
  for ( int i = FirstColumn(p_storage[0]) ; i <= LastColumn(p_storage[0]) ; i++ ) 
    p_mass_prec[i] = p_mass_prec_[0][i];
  
  int index = 1;
  for ( int i = FirstColumn(p_storage[1]) ; i <= LastColumn(p_storage[1]) ; i++ ) 
  {
    p_mass_prec[i] = p_mass_prec_[1][index];
    index++;
  }
  
  free(p_mass_prec_[0]);
  free(p_mass_prec_[1]);
#endif  

  //---------------------------------------------------------------------------
  // Build Stokes operator
  for ( int submatrix = 0 ; submatrix < 2 ; submatrix++ )
    axpby(NbOfEntries(&v_stiffness[ submatrix ]),
          tau0, v_mass[ submatrix ].entries,
          1./Re, v_stiffness[ submatrix ].entries);

  SparseMatrix_t* temp[] = {&v_stiffness[0], &v_mass[0]};
  ProcessDiagonal(2, temp);
  
  //---------------------------------------------------------------------------
  info( "\n""Reordering matrices""\n" );

#ifndef PRE_INCREMENTAL_ROTATIONAL
  SparseMatrix_t* spm[] = { &p_stiffness[0] };
  int NbOfMatrices = 1;
#else
  SparseMatrix_t* spm[] = { &p_stiffness[0], &p_mass[0] };
  int NbOfMatrices = 2;
#endif
  
  TaucsReorder(_reordering, NbOfMatrices, spm);
  ProcessDiagonal(NbOfMatrices, spm);

  //---------------------------------------------------------------------------
  info( "\n""Computing preconditioners""\n" );
  
  *p_stiffness_prec = GetPreconditioner(_drop_tolerance, &p_stiffness[0]);
  *v_mass_prec      = GetPreconditioner(_drop_tolerance, &v_mass[0]);
  *v_stiffness_prec = GetPreconditioner(_drop_tolerance, &v_stiffness[0]);
  
//#ifdef VERSION_Z && TWO_WAY
//  // copy back handlers
//  if (UsingMicroGrid())
//  {
//    p_stiffness_prec_m = p_stiffness_prec;
//    v_stiffness_prec_m = v_stiffness_prec;
//    v_mass_prec_m = v_mass_prec;
//  }
//  else
//  {
//    p_stiffness_prec_M = p_stiffness_prec;
//    v_stiffness_prec_M = v_stiffness_prec;
//    v_mass_prec_M = v_mass_prec;
//  }
//#endif
}
//=============================================================================
/** Update operators and preconditioners.
 
 Update the velocity stiffness operator from first to second order time discretization. Then reorder it as well as the velocity mass matrix as they share the same sparse pattern, and compute their preconditioners.
 */
//=============================================================================
void
UpdateOperators(
const double dt ) ///< time discretisation coefficient update
{
  PrintTitle("Updating operators for second order time discretization.");
  
  //---------------------------------------------------------------------------
  // Update velocity stiffness matrix and reorder
  for ( int submatrix = 0 ; submatrix < 2 ; submatrix++ )
    axpby(v_stiffness[ submatrix ].NbOfEntries, 1. / ( 2. * dt ),
      v_mass[ submatrix ].entries, 1., v_stiffness[ submatrix ].entries);

  //---------------------------------------------------------------------------
  info( "\n""Reordering matrices""\n" );
  
  SparseMatrix_t* spm[] = { &v_stiffness[0], &v_mass[0] };
  TaucsReorder(_reordering, 2, spm);
  ProcessDiagonal(2, spm);
  
  //---------------------------------------------------------------------------
  info( "Computing preconditioners\n" );
  
  FreeSparseMatrix(*v_stiffness_prec);
  FreeSparseMatrix(*v_mass_prec);
  
  *v_stiffness_prec = GetPreconditioner(_drop_tolerance, &v_stiffness[0]);
  *v_mass_prec      = GetPreconditioner(_drop_tolerance, &v_mass[0]);

  //---------------------------------------------------------------------------
#ifdef PRE_INCREMENTAL_ROTATIONAL
  // update pressure correction modified mass matrix
  // we use the same arrays as for the velocity
  for ( int submatrix = 0 ; submatrix < 2 ; submatrix++ )
    scal( p_mass[submatrix].NbOfEntries, 3./2., p_mass[submatrix].entries);
  
  // for pressure correction in rotational form we use the lumped operator
  
  // first get inverse
  p_mass_prec[0] = GetPreconditioner( &p_mass[0] );
  p_mass_prec[1] = GetPreconditioner( &p_mass[1] );
  
  for ( int i = FirstColumn(storage[0]) ; i <= LastColumn(storage[0]) ; i++ ) 
    p_mass_prec[i] = p_mass_prec[0][i];

  int index = 1;
  for ( int i = FirstColumn(storage[1]) ; i <= LastColumn(storage[1]) ; i++ ) 
  {
    p_mass_prec[i] = p_mass_prec[1][index];
    index++;
  }
#endif
}
