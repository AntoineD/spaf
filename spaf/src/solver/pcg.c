/** \file
 Conjugate gradient solver.
 
 Here the solver parameters can be read and a linear system can be solved for a sparse matrix.
 */

#include "includes.h"
#include "memory.h"
#include "strings.h"
#include "logging.h"
#include "linalg.h"
#include "smvp.h"
#include "parse.h"
#include "pcg.h"
#include "sparse_solve.h"
#include "smvp_legacy.h"

/// maximum number of iterations of the conjugate gradient
static int _IterMax;
/// convergence tolerance of the conjugate gradient \f$ \epsilon \f$ 
static double _tolerance;

//==============================================================================
/// Read conjugate gradient parameters.
//==============================================================================
void
ReadConjugateGradientParameters(
FILE* FileId )  ///< file identifier
{
  // read max number of iterations
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "iteration_max" ),
               "Missing solver parameter : ""iteration_max");
  _IterMax = Token2int( 3, '[', 0, USHRT_MAX,']' );
  
  // read convergence criterion \f$ \epsilon \f$ 
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "tolerance" ),
               "Missing solver parameter : ""tolerance");
  _tolerance = Token2double( 3, '[', DBL_EPSILON, 1.,']' );
  
  // check end of section
  TokenizeFileLine( FileId );
  assert_error( SectionIsEnd(),
               "Solver parameters : missing section = end");
  
  // print infos
  info("\nConjugate Gradient parameters :\n" ); 
  info("\tMaximum number of iterations = " INT_FMT "\n", _IterMax);
  info("\tTolerance = %g\n", _tolerance);
}
//==============================================================================
/** Process outside nodes.
 
 If a node is registered as outside in the OutsideNodes array, set the corresponding element of vector to 0.
 Note that OutsideNodes is 1 offset, whereas vector is 0 offset.
 */
//==============================================================================
static inline void
ProcessOutsideNodes(
const bool*   OutsideNodes,   ///< list of nodes not solved
const int*    invpermutation, ///< inverse permutation
const int     NbOfNodes,      ///< number of node to process
      double* vector )        ///< vector to be processed
{
  if ( !OutsideNodes ) return;
  
  if (invpermutation)
  {
#pragma omp parallel for \
default(none) \
shared(OutsideNodes,vector,invpermutation)
    for ( int i = 1 ; i <= NbOfNodes ; i++ )
      if ( OutsideNodes[i] ) vector[invpermutation[i-1]] = 0.;
  }
  else
  {
#pragma omp parallel for \
default(none) \
shared(OutsideNodes,vector,invpermutation)
    for ( int i = 1 ; i <= NbOfNodes ; i++ )
      if ( OutsideNodes[i] ) vector[i-1] = 0.;
  }
}
//==============================================================================
/** Preconditioned conjugate gradient solver.
 
 The convergence criterion is 
 \f$ \frac{< r_{k+1} , P^{-1} r_{k+1} >}{< r_0 , P^{-1} r_0 >} < \epsilon^2 \f$
 */
//==============================================================================
int
ConjugateGradient(
const bool* OutsideNodes, ///< list of nodes not solved
const SparseMatrix_t* A,  ///< matrix
const SparseMatrix_t* P,  ///< preconditioner
const double*         b,  ///< rhs vector
      double*         x ) ///< unknown vector, first guess on entry
{
  // number of free nodes
  int NbOfNodes = NbOfRowsStored(A);
  
  // some temporary vectors
  double *r = NULL,     // residual vector
         *d = NULL,     // direction of the update vector
         *temp = NULL;  // for temporary computations
  
  AllocVdouble(NbOfNodes-1, r);
  AllocVdouble(NbOfNodes-1, d);
  AllocVdouble(NbOfNodes-1, temp);
  
  // if the matrix has been permuted, we need to permute the rhs
  if ( A->permutation )
  {
    // use temp to permute \f$ x_0 \f$ in place
    copy(NbOfNodes, x, temp);
    permute(NbOfNodes, A->permutation, temp, x);
    
    // permute rhs in temp
    permute(NbOfNodes, A->permutation, b, temp);
  }
  else
    // temp holds b
    copy(NbOfNodes, b, temp);
  
  // \f$ r_0 = b - A x_0 \f$
  SparseMatrixVectorProduct('+', false, A, x, r);
  axpby(NbOfNodes, 1., temp, -1., r);
  
  ProcessOutsideNodes(OutsideNodes, A->invpermutation, NbOfNodes, r);
  
  // \f$ d_0 = P^{-1} r_0 \f$ 
  if (P) Solve(P, r, d);
  else   copy(NbOfNodes, r, d);

  ProcessOutsideNodes(OutsideNodes, A->invpermutation, NbOfNodes, d);
  
  // \f$ previous_scal_prod = < r_0 , P r_0 > \f$ 
  double previous_scal_prod = dot(NbOfNodes, r, d),
  
  // convergence_value = \f$ \epsilon^2 \f$
  convergence_value = previous_scal_prod * pow(_tolerance, 2);
  
  int iteration = 0;
  
  double test_value = previous_scal_prod;
  
  while ( test_value > convergence_value )
  {
    iteration++;
    
    // \f$ temp = A d_k \f$ 
    SparseMatrixVectorProduct(' ',false, A, d, temp);    
    
    // \f$ \alpha = < r_k , P r_k > / < d_k , A d_k > \f$ 
    double alpha = previous_scal_prod / dot(NbOfNodes, d, temp);
    
    // \f$ x_{k+1} = x_k + \alpha d_k \f$ 
    axpby(NbOfNodes, alpha, d, 1., x);
    
    // \f$ r_{k+1} = r_k - \alpha A d_k \f$ 
    axpby(NbOfNodes, -alpha, temp, 1., r);
    
    ProcessOutsideNodes(OutsideNodes, A->invpermutation, NbOfNodes, r);
    
    // \f$ temp = P^{-1} r_{k+1} \f$ 
    if (P) Solve(P, r, temp);
    else   copy(NbOfNodes, r, temp);

    ProcessOutsideNodes(OutsideNodes, A->invpermutation, NbOfNodes, temp);
    
    // \f$ < r_{k+1} , P^{-1} r_{k+1} > \f$ 
    double new_scal_prod = dot( NbOfNodes, r, temp );
    
    // \f$ \beta = < r_{k+1} , P^{-1} r_{k+1} > / < r_k , P^{-1} r_k > \f$ 
    double beta = new_scal_prod / previous_scal_prod;
    
    previous_scal_prod = new_scal_prod;
    
    // \f$ d_{k+1} = \beta d_k + P^{-1} r_{k+1} \f$ 
    axpby(NbOfNodes, 1., temp, beta, d);
    
    // the following is not necessary as d and temp have been processed already
    // ProcessOutsideNodes(OutsideNodes, NbOfNodes, d);
    
    assert_error(iteration <= _IterMax,
                 "Conjugate Gradient DID NOT CONVERGE, residual = "DBL_FMT,
                 test_value);
    
    test_value = new_scal_prod;
  }
  
  if ( A->permutation )
  {
    // use temp to permute back x in place
    copy(NbOfNodes, x, temp);
    permuteinv(NbOfNodes, A->permutation, temp, x);
  }
  
  free( r );
  free( d );
  free( temp );
  
  return iteration;
}
