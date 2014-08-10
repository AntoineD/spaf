/** \file
 Find the element that contains a given point, compute the point coordinates local to this element.
 */

#include "logging.h"
#include "parse.h"
#include "memory.h"
#include "mesh_search.h"

static int
// Maximum number of iterations for computing local coordinates of a point when
// the element is non linear.
_max_nb_of_iterations = 0,
// Maximum number of element to search in, when searching for an element
// which contains a point.
_max_nb_of_elements_to_search = 0;

static double
// Cutoff criterion for divergence of the iteration for computing
// the local coordinates of a point.
_iterations_tolerance_max = 0.,
// Convergence criterion of the iteration for computing the local coordinates of a point.
// Relative to the components of the non-linear transform.
_iterations_tolerance_min = 0.,
// Tolerances on the coordinates of a point when checking if it is in an element.
_coordinates_tolerance = 0.,
// Minimum determinant value under which an element is considered as having
// a too small volume.
_min_determinant_value = 0.;

static int
// table of the elements already checked
*_ElemChecked = NULL,
// tables to store the neighbours of the nodes of
// the elements already checked
*_NeighborNodesCurrentLayer = NULL,
*_NeighborNodesNextLayer = NULL;

#ifdef _OPENMP
/// File wise local variables for open mp
/// thread local private arrays for mesh search
static int **_ElemChecked_OMP = NULL,
            **_NeighborNodesCurrentLayer_OMP = NULL,
            **_NeighborNodesNextLayer_OMP = NULL;
//=============================================================================
/// Allocate required memory for openmp vectors used for mesh search.
//=============================================================================
static void
AllocOMP( void )
{
#pragma omp parallel \
default(none) \
shared(_max_nb_of_elements_to_search) \
shared(_ElemChecked_OMP,_NeighborNodesCurrentLayer_OMP,_NeighborNodesNextLayer_OMP)
{
#pragma omp master
{
  // get the number of threads
  int NbOfThreads = omp_get_num_threads();
  
  // allocate arrays
  AllocMint(
    NbOfThreads,_max_nb_of_elements_to_search,_ElemChecked_OMP);

  AllocMint(
    NbOfThreads,_max_nb_of_elements_to_search,_NeighborNodesCurrentLayer_OMP);

  AllocMint(
    NbOfThreads,_max_nb_of_elements_to_search,_NeighborNodesNextLayer_OMP);
}
}
}
#endif
//=============================================================================
/// Read mesh search parameters
//=============================================================================
void
ReadMeshSearchParameters(
FILE* FileId )  ///< file identifier
{
  // read max_nb_of_iterations
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("max_nb_of_iterations"),
               "Missing mesh search parameter : ""max_nb_of_iterations");
  _max_nb_of_iterations = Token2int(3,'[',0,USHRT_MAX,']');

  // read max_nb_of_elements_to_search
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("max_nb_of_elements_to_search"),
               "Missing mesh search parameter : ""max_nb_of_elements_to_search");
  _max_nb_of_elements_to_search = Token2int(3,'[',0,USHRT_MAX,']');

  // read iterations_tolerance_max
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("iterations_tolerance_max"),
               "Missing mesh search parameter : ""iterations_tolerance_max");
  _iterations_tolerance_max = Token2double(3,'[',0.,DBL_MAX,']');

  // read iterations_tolerance_min
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("iterations_tolerance_min"),
               "Missing mesh search parameter : ""iterations_tolerance_min");
  _iterations_tolerance_min = Token2double(3,'[',0.,DBL_MAX,']');

  // read coordinates_tolerance
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("coordinates_tolerance"),
               "Missing mesh search parameter : ""coordinates_tolerance");
  _coordinates_tolerance = Token2double(3,'[',0.,DBL_MAX,']');

  // read _min_determinant_value
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("min_determinant_value"),
               "Missing mesh search parameter : ""min_determinant_value");
  _min_determinant_value = Token2double(3,'[',0.,DBL_MAX,']');
  
  // check end of section
  TokenizeFileLine(FileId);
  assert_error( SectionIsEnd(),
               "Mesh search parameters : missing section = end");

  // print parameters
  info("\nMesh search parameters :\n" ); 
  info("\tmax_nb_of_iterations         = " INT_FMT "\n",
       _max_nb_of_iterations);
  info("\tmax_nb_of_elements_to_search = " INT_FMT "\n",
       _max_nb_of_elements_to_search);
  info("\titerations_tolerance_max     = " "%g" "\n",
       _iterations_tolerance_max);
  info("\titerations_tolerance_min     = " "%g" "\n",
       _iterations_tolerance_min);
  info("\tcoordinates_tolerance        = " "%g" "\n",
       _coordinates_tolerance);
  info("\tmin_determinant_value        = " "%g" "\n",
       _min_determinant_value);
  
  // allocate for filewise static arrays
#ifndef _OPENMP
  AllocVint(_max_nb_of_elements_to_search, _ElemChecked);
  AllocVint(_max_nb_of_elements_to_search, _NeighborNodesCurrentLayer);
  AllocVint(_max_nb_of_elements_to_search, _NeighborNodesNextLayer);
#else
  AllocOMP();
#endif
}
//=============================================================================
/** Find local position of a point in an element.
 
 Computes local coordinates in the reference element of the point
 which coordinates are in point. Then it check whether or not those
 coordinates refer to a point which is inside, by a tolerance given in the
 variables _coordinates_tolerance_*.
 Returns true and the local coordinates if the point is inside, return false
 otherwise.
 
 A linearization of the mapping between physical element and reference element
 is computed, Newton iterations are used if necessary, and the local
 coordinates are computed. The iterations are controlled by the variable 
 _max_nb_of_iterations which is the maximum number of iterations and the
 variable _iterations_tolerance_min whc=ich measure the convergence of the iteration
 process. 
 
 The volume of the element, or the determinant of the mapping is check against
 the variable _min_determinant_value.
*/
//=============================================================================
bool
MeshGetLocalPosition(
const mesh_t* mesh,         ///< mesh structure
const int    el,           ///< element index
const double  point[4],     ///< point global coordiantes
      double  LocalPos[4] ) ///< point local coordinates
{
  double x, y, z, x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4, x5, y5, z5,
  x6, y6 , z6, x7, y7, z7, x8, y8, z8, x9, y9, z9, x10, y10, z10;  /* coordinates of the nodes */
  double j11, j12, j13, j21, j22, j23, j31, j32, j33;    /* the components of the jacobean */
  double f1, f2, f3;                                     /* the components of the non-linear transform  */
  double l1, l2, l3, l4, determ, eps11, eps12, eps21,
  eps22;                                          /* l is the linear (current) guess;
   determ is working variable*/
  double norm1;                                          /* norm of the correction */
  int n = 1;                                             /* counter of the iterations */
  
  eps11 = 1.;
  eps12 = 2.;  /* Cutoff value for quick & dirty element check */
  eps21 = _coordinates_tolerance;
  eps22 = 1. + _coordinates_tolerance;
  
  /***  Load the coordinates of the nodal points of the element and the Lagrange point  ***/
  
  x = point[1];
  y = point[2];
  z = point[3];
  // -------------------------------------------------------------------
  // Note:
  // I CHANGED NUMBERING HERE TO CONSIDER FIRST 4 NODES PRESSURE NODES
  // Carolina May 00
  // -------------------------------------------------------------------
  x1 = mesh->points[ mesh->ConTab[ el ][ 1 ] ][ 1 ];
  x3 = mesh->points[ mesh->ConTab[ el ][ 2 ] ][ 1 ];
  x5 = mesh->points[ mesh->ConTab[ el ][ 3 ] ][ 1 ];
  x10 = mesh->points[ mesh->ConTab[ el ][ 4 ] ][ 1 ];
  x2 = mesh->points[ mesh->ConTab[ el ][ 5 ] ][ 1 ];
  x4 = mesh->points[ mesh->ConTab[ el ][ 6 ] ][ 1 ];
  x6 = mesh->points[ mesh->ConTab[ el ][ 7 ] ][ 1 ];
  x7 = mesh->points[ mesh->ConTab[ el ][ 8 ] ][ 1 ];
  x8 = mesh->points[ mesh->ConTab[ el ][ 9 ] ][ 1 ];
  x9 = mesh->points[ mesh->ConTab[ el ][ 10 ] ][ 1 ];
  
  y1 = mesh->points[ mesh->ConTab[ el ][ 1 ] ][ 2 ];
  y3 = mesh->points[ mesh->ConTab[ el ][ 2 ] ][ 2 ];
  y5 = mesh->points[ mesh->ConTab[ el ][ 3 ] ][ 2 ];
  y10 = mesh->points[ mesh->ConTab[ el ][ 4 ] ][ 2 ];
  y2 = mesh->points[ mesh->ConTab[ el ][ 5 ] ][ 2 ];
  y4 = mesh->points[ mesh->ConTab[ el ][ 6 ] ][ 2 ];
  y6 = mesh->points[ mesh->ConTab[ el ][ 7 ] ][ 2 ];
  y7 = mesh->points[ mesh->ConTab[ el ][ 8 ] ][ 2 ];
  y8 = mesh->points[ mesh->ConTab[ el ][ 9 ] ][ 2 ];
  y9 = mesh->points[ mesh->ConTab[ el ][ 10 ] ][ 2 ];
  
  z1 = mesh->points[ mesh->ConTab[ el ][ 1 ] ][ 3 ];
  z3 = mesh->points[ mesh->ConTab[ el ][ 2 ] ][ 3 ];
  z5 = mesh->points[ mesh->ConTab[ el ][ 3 ] ][ 3 ];
  z10 = mesh->points[ mesh->ConTab[ el ][ 4 ] ][ 3 ];
  z2 = mesh->points[ mesh->ConTab[ el ][ 5 ] ][ 3 ];
  z4 = mesh->points[ mesh->ConTab[ el ][ 6 ] ][ 3 ];
  z6 = mesh->points[ mesh->ConTab[ el ][ 7 ] ][ 3 ];
  z7 = mesh->points[ mesh->ConTab[ el ][ 8 ] ][ 3 ];
  z8 = mesh->points[ mesh->ConTab[ el ][ 9 ] ][ 3 ];
  z9 = mesh->points[ mesh->ConTab[ el ][ 10 ] ][ 3 ];
  
  /*
   info("\n x1=%lf\t y1=%lf \t z1=%lf \n x2=%lf\t y2=%lf  \t z2=%lf \n x3=%lf \t y3=%lf  \t z3=%lf \n x4=%lf \t y4=%lf \t z4=%lf \n x5=%lf\t y5=%lf \t z5=%lf \n x6=%lf\t y6=%lf  \t z6=%lf \n x7=%lf\t y7=%lf \t z7=%lf \n x8=%lf\t y8=%lf  \t z8=%lf  \n x9=%lf\t y9=%lf \t z9=%lf \n x10=%lf\t y10=%lf  \t z10=%lf \n x=%lf \t y=%lf \t z=%lf\n",x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x,y,z);
   */
  
  /*
   *  Compute the barycentric coordinates by Cramers Rule.
   */
  
  determ = ( x1 - x10 ) * ( y5 - y10 ) * ( z3 - z10 ) + ( x3 - x10 ) * ( y1 - y10 ) * ( z5 - z10 ) +
  ( x5 - x10 ) * ( y3 - y10 ) * ( z1 - z10 ) - ( x3 - x10 ) * ( y5 - y10 ) * ( z1 - z10 ) -
  ( x5 - x10 ) * ( y1 - y10 ) * ( z3 - z10 ) - ( x1 - x10 ) * ( y3 - y10 ) * ( z5 - z10 );
  
  assert_error( fabs( determ ) > _min_determinant_value,
          "volume of linear element "INT_FMT" is almost 0.", el );
  
  l1 = ( ( x - x10 ) * ( y5 - y10 ) * ( z3 - z10 ) + ( x3 - x10 ) * ( y - y10 ) * ( z5 - z10 ) +
        ( x5 - x10 ) * ( y3 - y10 ) * ( z - z10 ) - ( x3 - x10 ) * ( y5 - y10 ) * ( z - z10 ) -
        ( x5 - x10 ) * ( y - y10 ) * ( z3 - z10 ) - ( x - x10 ) * ( y3 - y10 ) * ( z5 - z10 ) ) / determ;
  
  l2 = ( ( x1 - x10 ) * ( y - y10 ) * ( z3 - z10 ) + ( x3 - x10 ) * ( y1 - y10 ) * ( z - z10 ) +
        ( x - x10 ) * ( y3 - y10 ) * ( z1 - z10 ) - ( x3 - x10 ) * ( y - y10 ) * ( z1 - z10 ) -
        ( x - x10 ) * ( y1 - y10 ) * ( z3 - z10 ) - ( x1 - x10 ) * ( y3 - y10 ) * ( z - z10 ) ) / determ;
  
  l3 = ( ( x1 - x10 ) * ( y5 - y10 ) * ( z - z10 ) + ( x - x10 ) * ( y1 - y10 ) * ( z5 - z10 ) +
        ( x5 - x10 ) * ( y - y10 ) * ( z1 - z10 ) - ( x - x10 ) * ( y5 - y10 ) * ( z1 - z10 ) -
        ( x5 - x10 ) * ( y1 - y10 ) * ( z - z10 ) - ( x1 - x10 ) * ( y - y10 ) * ( z5 - z10 ) ) / determ;
  
  l4 = l1 + l2 + l3;
  
  
  /***  Check if the "linear" local coordinates are too far away from the interval [-1,1] ***/
  
  if ( l1 > eps12 || l2 > eps12 ||
      l3 >  eps12 || ( 1. - l4 ) > eps12 ||
      l1 < -eps11 || l2 < -eps11 ||
      l3 < -eps11 || ( 1. - l4 ) < -eps11 )
  {
    /*info("\n Wrong initial guess for %lf %lf %lf - l1=%lf \t l2=%lf \t l3=%lf \t l4=%lf",Lag_node_x,
     Lag_node_y, Lag_node_z,l1,l2,l3,l4);  */ 
    return false;
  }
  
  /* Functions f1, f2 & f3 are based on isoparametric expressions for x,y,z in terms of barycentric co-ords, i.e.
   f1 = sum (i=1...10) X_i*phi_i(xi, eta, zeta) - x */
  f1 = (double) ( x1 * l1 * ( 2. * l1 - 1. ) + 4. * x2 * l1 * l3 + x3 * l3 * ( 2. * l3 - 1. ) + 4. * x4 * l2 * l3 + x5 * l2 * ( 2. * l2 - 1. ) + 4. * x6 * l1 * l2 + 4. * x7 * l1 * ( 1. - l4 ) + 4. * x8 * l3 * ( 1. - l4 ) + 4. * x9 * l2 * ( 1. - l4 ) + x10 * ( 1. - 2. * l4 ) * ( 1. - l4 ) - x );
  
  f2 = (double) ( y1 * l1 * ( 2. * l1 - 1. ) + 4. * y2 * l1 * l3 + y3 * l3 * ( 2. * l3 - 1. ) + 4. * y4 * l2 * l3 + y5 * l2 * ( 2. * l2 - 1. ) + 4. * y6 * l1 * l2 + 4. * y7 * l1 * ( 1. - l4 ) + 4. * y8 * l3 * ( 1. - l4 ) + 4. * y9 * l2 * ( 1. - l4 ) + y10 * ( 1. - 2. * l4 ) * ( 1. - l4 ) - y );
  
  f3 = (double) ( z1 * l1 * ( 2. * l1 - 1. ) + 4. * z2 * l1 * l3 + z3 * l3 * ( 2. * l3 - 1. ) + 4. * z4 * l2 * l3 + z5 * l2 * ( 2. * l2 - 1. ) + 4. * z6 * l1 * l2 + 4. * z7 * l1 * ( 1. - l4 ) + 4. * z8 * l3 * ( 1. - l4 ) + 4. * z9 * l2 * ( 1. - l4 ) + z10 * ( 1. - 2. * l4 ) * ( 1. - l4 ) - z );
  
  norm1 = (double) sqrt( f1 * f1 + f2 * f2 + f3 * f3 );
  
  /***  Loop for non-linear iterations  ***/
  
  while ( norm1 > _iterations_tolerance_min  && n <= _iterations_tolerance_max )
  {
    
    /***  Compute the Jacobean ***/
    
    j11 =  ( x1 * ( 4. * l1 - 1. ) + 4. * x2 * l3 + 4. * x6 * l2 + 4. * x7 * ( 1 - l4 - l1 ) - 4. * x8 * l3 - 4. * x9 * l2 +
            x10 * ( 4. * l4 - 3. ) );
    j12 =  ( 4. * x4 * l3 + x5 * ( 4. * l2 - 1. ) + 4. * x6 * l1 - 4. * x7 * l1 - 4. * x8 * l3 + 4. * x9 * ( 1 - l4 - l2 ) +
            x10 * ( 4. * l4 - 3. ) );
    j13 =  ( 4. * x2 * l1 + x3 * ( 4. * l3 - 1. ) + 4. * x4 * l2 - 4. * x7 * l1 + 4. * x8 * ( 1 - l4 - l3 ) - 4. * x9 * l2 +
            x10 * ( 4. * l4 - 3. ) );
    j21 =  ( y1 * ( 4. * l1 - 1. ) + 4. * y2 * l3 + 4. * y6 * l2 + 4. * y7 * ( 1 - l4 - l1 ) - 4. * y8 * l3 - 4. * y9 * l2 +
            y10 * ( 4. * l4 - 3. ) );
    j22 =  ( 4. * y4 * l3 + y5 * ( 4. * l2 - 1. ) + 4. * y6 * l1 - 4. * y7 * l1 - 4. * y8 * l3 + 4. * y9 * ( 1 - l4 - l2 ) +
            y10 * ( 4. * l4 - 3. ) );
    j23 =  ( 4. * y2 * l1 + y3 * ( 4. * l3 - 1. ) + 4. * y4 * l2 - 4. * y7 * l1 + 4. * y8 * ( 1 - l4 - l3 ) - 4. * y9 * l2 +
            y10 * ( 4. * l4 - 3. ) );
    j31 =  ( z1 * ( 4. * l1 - 1. ) + 4. * z2 * l3 + 4. * z6 * l2 + 4. * z7 * ( 1 - l4 - l1 ) - 4. * z8 * l3 - 4. * z9 * l2 +
            z10 * ( 4. * l4 - 3. ) );
    j32 =  ( 4. * z4 * l3 + z5 * ( 4. * l2 - 1. ) + 4. * z6 * l1 - 4. * z7 * l1 - 4. * z8 * l3 + 4. * z9 * ( 1 - l4 - l2 ) +
            z10 * ( 4. * l4 - 3. ) );
    j33 =  ( 4. * z2 * l1 + z3 * ( 4. * l3 - 1. ) + 4. * z4 * l2 - 4. * z7 * l1 + 4. * z8 * ( 1 - l4 - l3 ) - 4. * z9 * l2 +
                    z10 * ( 4. * l4 - 3. ) );
    
    /***  Calculate the components of the transform - f1,f2,f3   ***/
    
    f1 =  ( x1 * l1 * ( 2. * l1 - 1. ) + 4. * x2 * l1 * l3 + x3 * l3 * ( 2. * l3 - 1. ) + 4. * x4 * l2 * l3 + x5 * l2 * ( 2. * l2 - 1. ) + 4. * x6 * l1 * l2 +
                   4. * x7 * l1 * ( 1. - l4 ) + 4. * x8 * l3 * ( 1. - l4 ) + 4. * x9 * l2 * ( 1. - l4 ) + x10 * ( 1. - 2. * l4 ) * ( 1. - l4 ) - x );
    
    f2 =  ( y1 * l1 * ( 2. * l1 - 1. ) + 4. * y2 * l1 * l3 + y3 * l3 * ( 2. * l3 - 1. ) + 4. * y4 * l2 * l3 + y5 * l2 * ( 2. * l2 - 1. ) + 4. * y6 * l1 * l2 +
                   4. * y7 * l1 * ( 1. - l4 ) + 4. * y8 * l3 * ( 1. - l4 ) + 4. * y9 * l2 * ( 1. - l4 ) + y10 * ( 1. - 2. * l4 ) * ( 1. - l4 ) - y );
    
    f3 =  ( z1 * l1 * ( 2. * l1 - 1. ) + 4. * z2 * l1 * l3 + z3 * l3 * ( 2. * l3 - 1. ) + 4. * z4 * l2 * l3 + z5 * l2 * ( 2. * l2 - 1. ) + 4. * z6 * l1 * l2 +
                   4. * z7 * l1 * ( 1. - l4 ) + 4. * z8 * l3 * ( 1. - l4 ) + 4. * z9 * l2 * ( 1. - l4 ) + z10 * ( 1. - 2. * l4 ) * ( 1. - l4 ) - z );
    
    /***  Solve the linearized system        ***/
    
    determ = j11 * j22 * j33 + j31 * j12 * j23 + j13 * j21 * j32 - j13 * j22 * j31 - j11 * j23 * j32 - j12 * j21 * j33;
    
    assert_error( fabs( determ ) > _min_determinant_value,
            "jacobian is almost 0" );
    
    l1 =  ( -( f1 * j22 * j33 + f3 * j12 * j23 + j13 * f2 * j32 - j13 * j22 * f3 - f1 * j23 * j32 - j12 * f2 * j33 ) /
           determ + l1 );
    
    l2 =  ( -( j11 * f2 * j33 + j31 * f1 * j23 + j13 * j21 * f3 - j13 * f2 * j31 - j11 * j23 * f3 - f1 * j21 * j33 ) /
           determ + l2 );
    
    l3 =  ( -( j11 * j22 * f3 + j31 * j12 * f2 + f1 * j21 * j32 - f1 * j22 * j31 - j11 * f2 * j32 - j12 * j21 * f3 ) /
           determ + l3 );
    
    l4 =  ( l1 + l2 + l3 );
    
    /***   Check if the current guess is far from the allowed coordinates  ***/
    
    if ( l1 > eps12 || l2 > eps12 ||
         l3 > eps12 || ( 1. - l4 ) > eps12 ||
         l1 < -eps11 || l2 < -eps11 ||
         l3 < -eps11 || ( 1. - l4 ) < -eps11 )
      return false;
    
    norm1 = (double) sqrt( f1 * f1 + f2 * f2 + f3 * f3 );
    
    /*** Check if the current norm is larger than the maximum allowed; if yes - divergence is supposed ***/
    
    if ( norm1 > _iterations_tolerance_max ) return false;
    else n++;
  }
  
  if ( n > _max_nb_of_iterations ) return false;
  
  /***  Check if the nonlinear local coordinates are between 0 and 1  ***/
  
  LocalPos[1] = l1;
  LocalPos[2] = l2;
  LocalPos[3] = l3;
  
  if ( l1 > eps22 ||
       l2 > eps22 ||
       l3 > eps22 ||
       ( 1. - l4 ) > eps22 ||
       l1 < -eps21 ||
       l2 < -eps21 ||
       l3 < -eps21 ||
       ( 1. - l4 ) < -eps21 )
    return false;
  
  return true;
}

//=============================================================================
/** Finds local coordinates of a point and the element that contains it.

 Search in all elements that contains the closest node guess, then extend to
 neighbor elements of those and so on until the container element is found
 or the search stop criterium is reached.
 Returns the index of the element which contains the point,
 otherwise returns 0.
*/
//=============================================================================
static int
FindElemContainingPoint(
const mesh_t* mesh,             ///< mesh structure
const double  point[4],         ///< point
const int    StartingNode,     ///< node to start the search from
      int*   NbOfElemScanned ) ///< number of element scanned
{ 
  assert_error(StartingNode > 0,"bad node");
  
  // not used but required by MeshGetLocalPosition
  double LocalPos[4] = {0.,0.,0.,0.};

#ifdef _OPENMP
  // get thread index
  int thread = omp_get_thread_num()+1;
  
  // we have to hide the definition of the local filewise arrays and use function local arrays that are local to the thread
  int
  *_ElemChecked = _ElemChecked_OMP[thread],
  *_NeighborNodesCurrentLayer = _NeighborNodesCurrentLayer_OMP[thread],
  *_NeighborNodesNextLayer = _NeighborNodesNextLayer_OMP[thread];
#endif
  
  int
    NbOfElemChecked = 0,
    NbOfNeighborNodesNextLayer = 0,
    // we start with only one neighbor node
    NbOfNeighborNodesCurrentLayer = 1;
  
  // and this node is the starting node
  _NeighborNodesCurrentLayer[1] = StartingNode;
  
  // switch to stop search when the max number of elements has been reached
  bool StopSearch = false;

  // make sure we stay within the range of elements in the mesh
  while ( NbOfElemChecked <= _max_nb_of_elements_to_search )
  {
    // go over the neighbor nodes
    for ( int iNeighborNode = 1 ;
          iNeighborNode <= NbOfNeighborNodesCurrentLayer ;
          iNeighborNode++ )
    {
      int NeighborNode = _NeighborNodesCurrentLayer[ iNeighborNode ];      
      
      // Look over the elements containing the current neighbor node
      for ( int iElemWithNode = 1 ;
            iElemWithNode <= NbOfEntriesInRow(mesh->InvConTab, NeighborNode-1) ;
            iElemWithNode++ )
      {
        int element = EntryNode(mesh->InvConTab, NeighborNode-1, iElemWithNode-1);

        // Has this element already been checked ?
        // Search its entry in checked element list.
        int iElemChecked;
        
        for ( iElemChecked = 1 ; iElemChecked <= NbOfElemChecked ; iElemChecked++ )
          if ( element == _ElemChecked[ iElemChecked ] ) break;
            
        // If already been checked then jump to next element containing node
        if ( iElemChecked != NbOfElemChecked+1 ) continue;

        NbOfElemChecked++;

        // Does this element contains the point ?
        if ( MeshGetLocalPosition( mesh, element, point, LocalPos ) == true )
        {
          if ( NbOfElemScanned != NULL ) *NbOfElemScanned += NbOfElemChecked;
          return element;
        }

        // have we reached the max number of elements to search in ?
        if ( NbOfElemChecked > _max_nb_of_elements_to_search )
        {
          StopSearch = true;
          break;
        }
        
        // mark this element as checked and add its vertices to the list of neighbor nodes.
        _ElemChecked[ NbOfElemChecked ] = element;

        // Add the element vertices to the neighbour table for the next
        // search if they are not already there.
        for ( int PreNode = 1 ; PreNode <= PRE_NODES_PER_EL ; PreNode++ )
        {
          int node = mesh->ConTab[ element ][ PreNode ];

          if ( node == NeighborNode ) continue;

          int i_NeighborNode;
          
          for ( i_NeighborNode = 1 ; i_NeighborNode <= NbOfNeighborNodesNextLayer ; i_NeighborNode++ )
            if ( node == _NeighborNodesNextLayer[ i_NeighborNode ] ) break;
          
          if ( i_NeighborNode != NbOfNeighborNodesNextLayer+1 ) continue;

          NbOfNeighborNodesNextLayer++;
          assert_error( NbOfNeighborNodesNextLayer < _max_nb_of_elements_to_search,
                       "Increase _max_nb_of_elements_to_search");
          _NeighborNodesNextLayer[ NbOfNeighborNodesNextLayer ] = node;
        }
      }
      // get out if searching is over
      if ( StopSearch ) break;
    }
    // get out if searching is over
    if ( StopSearch ) break;

    // copy next layer into current one
    int* temp = _NeighborNodesCurrentLayer;
    _NeighborNodesCurrentLayer = _NeighborNodesNextLayer;
    _NeighborNodesNextLayer = temp;
    
    NbOfNeighborNodesCurrentLayer = NbOfNeighborNodesNextLayer;
    NbOfNeighborNodesNextLayer = 0;
  }

// for version z, this is not necessary as nodes outside outer domain are taken care of
#ifndef VERSION_Z
  warning("max number of elements to search reached at node "INT_FMT", "
          INT_FMT" elements have been checked\n",
          StartingNode, NbOfElemChecked );
#endif
  
  // Update scanned element count, clean memory up and exit.
  if ( NbOfElemScanned != NULL ) *NbOfElemScanned += NbOfElemChecked;
  
  return 0;
}
//==============================================================================
/** Search for an element containing a point.
 With BB search and spiral search.
 ClosestNode_ptr is the previous closest node, it is returned with a new value
 only if the point is found in an element.
 */
//==============================================================================
int
SearchElemContainingPoint(
const mesh_t* mesh,             ///< pointer to mesh structure
      bool    UseNNS,           ///< switch to search in bb or not
const double  point[4],         ///< the point of interest
      int*   ClosestNode_ptr,  ///< potential node to start search from
      int*   ElemSearchedNb )  ///< total number of element searched
{
  // is the point inside bb ?
  if ( IsInsideBoundingBox(mesh->nns, point) ) return 0;

  // use the provided nearest node
  int ClosestNode = *ClosestNode_ptr;
  
  // unless we have to search for a new one
  if ( UseNNS == true ) ClosestNode = DoNNS(mesh->nns, point);
  
  // if the point has not been found we exit
  if ( ClosestNode == -1 ) return 0;
  
  // search for the element containing the point
  int ElemContainer =
  FindElemContainingPoint( mesh, point, ClosestNode, ElemSearchedNb );
  
  // if the point is inside, we keep closest point found by the bb search
  if ( ElemContainer != 0 ) *ClosestNode_ptr = ClosestNode;
  
  return ElemContainer;
}
