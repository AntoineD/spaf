/** \file
 Choleski decomposition solver wrapper.
 
 This is a wrapper to a library found on the web. This library is supposed to have a fast Choleski solver implementation but I cannpot get it to work, so be sure to fix it before going further.
 */

#include "cholmod_wrapper.h"
#include "memory.h"
#include "linalg.h"

// defined a minimal cholmod sparse matrix structure
typedef struct
{
  int n, *nz, *i, *p;
  double *x;
} cholmod_factor;

#define REAL
#define LL
#define Int int
#define cholmod_dense double
#define ASSERT assert
#define PREFIX r_
#include "t_cholmod_lsolve.h"
#define LL
#include "t_cholmod_ltsolve.h"

void
CholmodSolve(
const SparseMatrix_t* m,
const double* b,
      double* x )
{
  cholmod_factor t;
  t.p = m->offset;
  t.i = m->columns;
  t.x = m->entries;
  t.n = m->NbOfRowsStored;

  AllocVint(t.n-1, t.nz);

  for ( int i = 0 ; i < t.n ; i++ )
  {
    t.nz[i] = t.p[i+1] - t.p[i];
    
    for (int j = t.p[i]; j < t.p[i+1]-1; j++)
    {
      assert( t.i[j+1] > t.i[j] );
    }
  }
  
  copy(t.n, b, x);
  
  r_ll_lsolve_1(&t, x);
  r_ll_ltsolve_1(&t, x);
  
  free(t.nz);
}

