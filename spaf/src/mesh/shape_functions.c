/** \file
 Shape functions computations.
 
 Here you can:
 - compute the values of the shape functions or the derivatives at a local point
 - interpolate a field at a local point
 - compute the weighted jacobian and the inverse jacobian at a local point
 - same but for a sub element inside a particle
 - compute the gradient at a local point
 */

#include "includes.h"
#include "logging.h"
#include "linalg.h"
#include "memory.h"
#include "shape_functions.h"

// Warning printed if Jacobian det <= this value
#define JAC_DET_WARN DBL_EPSILON

//==============================================================================
/// Get the values of the shape function at a local point.
//==============================================================================
void
GetShapeFunctionsAtPoint(
const int     order,          ///< order of the shape function
const double  point[3],       ///< local coordinates of a point
      double* ShapeFunction ) ///< shape functions values at this point
{
  double
  xi   = point[0],
  eta  = point[1],
  zeta = point[2],
  sum  = xi + eta + zeta;
  
  if ( order == 1 )
  {
    ShapeFunction[ 0 ] = xi;
    ShapeFunction[ 1 ] = zeta;
    ShapeFunction[ 2 ] = eta;
    ShapeFunction[ 3 ] = 1. - sum;
  }
  else if ( order == 5 )
  {
    ShapeFunction[ 0 ] = xi   * ( 2. * xi   - 1. );
    ShapeFunction[ 1 ] = zeta * ( 2. * zeta - 1. );
    ShapeFunction[ 2 ] = eta  * ( 2. * eta  - 1. );
    ShapeFunction[ 3 ] = ( 1. - 2. * sum ) * ( 1. - sum );
    ShapeFunction[ 4 ] = 4. * xi   * zeta;
    ShapeFunction[ 5 ] = 4. * zeta * eta;
    ShapeFunction[ 6 ] = 4. * xi   * eta;
    ShapeFunction[ 7 ] = 4. * xi   * ( 1. - sum );
    ShapeFunction[ 8 ] = 4. * zeta * ( 1. - sum );
    ShapeFunction[ 9 ] = 4. * eta  * ( 1. - sum );    
  }
  else
    error("Bad order");
}
//==============================================================================
/** Get the values of the derivatives of the shape functions at a local point.

 This function computes the 3 components.
 */
//==============================================================================
void
GetDerivativeOfShapeFunctionsAtPoint(
const int       order,                    ///< order of the shape function
const double    point[3],                 ///< local coordinates of a point
      double**  ShapeFunctionDerivative ) ///< derivative shape functions values at this point
{
  double
  xi   = point[0],
  eta  = point[1],
  zeta = point[2],
  sum  = xi + eta + zeta;
  
  if ( order == 1 )
  {
    ShapeFunctionDerivative[ 1 ][ 1 ] = 1.;
    ShapeFunctionDerivative[ 2 ][ 1 ] = 0.;
    ShapeFunctionDerivative[ 3 ][ 1 ] = 0.;
    ShapeFunctionDerivative[ 4 ][ 1 ] = -1.;
    
    ShapeFunctionDerivative[ 1 ][ 2 ] = 0.;
    ShapeFunctionDerivative[ 2 ][ 2 ] = 0.;
    ShapeFunctionDerivative[ 3 ][ 2 ] = 1.;
    ShapeFunctionDerivative[ 4 ][ 2 ] = -1.;
    
    ShapeFunctionDerivative[ 1 ][ 3 ] = 0.;
    ShapeFunctionDerivative[ 2 ][ 3 ] = 1.;
    ShapeFunctionDerivative[ 3 ][ 3 ] = 0.;
    ShapeFunctionDerivative[ 4 ][ 3 ] = -1.;          
  }
  else if ( order == 5 )
  {
    ShapeFunctionDerivative[ 1 ][ 1 ] = 4. * xi - 1.;
    ShapeFunctionDerivative[ 7 ][ 1 ] = 4. * eta;
    ShapeFunctionDerivative[ 3 ][ 1 ] = 0.;
    ShapeFunctionDerivative[ 6 ][ 1 ] = 0.;
    ShapeFunctionDerivative[ 2 ][ 1 ] = 0.;
    ShapeFunctionDerivative[ 5 ][ 1 ] = 4. * zeta;
    ShapeFunctionDerivative[ 8 ][ 1 ] = 4. * ( 1. - sum - xi );
    ShapeFunctionDerivative[ 10 ][ 1 ]= -4. * eta;
    ShapeFunctionDerivative[ 9 ][ 1 ] = -4. * zeta;
    ShapeFunctionDerivative[ 4 ][ 1 ] = 4. * sum - 3.;
    
    ShapeFunctionDerivative[ 1 ][ 2 ] = 0.;
    ShapeFunctionDerivative[ 7 ][ 2 ] = 4. * xi;
    ShapeFunctionDerivative[ 3 ][ 2 ] = 4. * eta - 1.;
    ShapeFunctionDerivative[ 6 ][ 2 ] = 4. * zeta;
    ShapeFunctionDerivative[ 2 ][ 2 ] = 0.;
    ShapeFunctionDerivative[ 5 ][ 2 ] = 0.;
    ShapeFunctionDerivative[ 8 ][ 2 ] = -4. * xi ;
    ShapeFunctionDerivative[ 10 ][ 2 ]= 4. * ( 1. - sum - eta );
    ShapeFunctionDerivative[ 9 ][ 2 ] = -4. * zeta;
    ShapeFunctionDerivative[ 4 ][ 2 ] = 4. * sum - 3.;
    
    ShapeFunctionDerivative[ 1 ][ 3 ] = 0.;
    ShapeFunctionDerivative[ 7 ][ 3 ] = 0.;
    ShapeFunctionDerivative[ 3 ][ 3 ] = 0.;
    ShapeFunctionDerivative[ 6 ][ 3 ] = 4. * eta;
    ShapeFunctionDerivative[ 2 ][ 3 ] = 4. * zeta - 1.;
    ShapeFunctionDerivative[ 5 ][ 3 ] = 4. * xi;
    ShapeFunctionDerivative[ 8 ][ 3 ] = -4. * xi;
    ShapeFunctionDerivative[ 10 ][ 3 ]= -4. * eta;
    ShapeFunctionDerivative[ 9 ][ 3 ] = 4. * ( 1. - sum - zeta );
    ShapeFunctionDerivative[ 4 ][ 3 ] = 4. * sum - 3.;          
  }
  else
    error("Bad order");
}
//==============================================================================
/// Second order interpolation in 3D.
/**
 Interpolate at second order from the variable input at a point which is in an
 element. The local connectivity of this element is given in ConTab, the local
 position of the point relatively to the element is in LocalPos.
 The result is returned to the caller.
 */
//==============================================================================
double
InterpolateOrder2(
const int     ConTab[11],   ///< element local connectivity
const double  LocalPos[4],  ///< local position of a point
const double* input )       ///< interpolate from this variable
{
  double s = LocalPos[1] + LocalPos[2] + LocalPos[3];
  
  double u1 = input[ ConTab[ 1 ] ],
  u3 = input[ ConTab[ 2 ] ],
  u5 = input[ ConTab[ 3 ] ],
  u10 =input[ ConTab[ 4 ] ],
  u2 = input[ ConTab[ 5 ] ],
  u4 = input[ ConTab[ 6 ] ],
  u6 = input[ ConTab[ 7 ] ],
  u7 = input[ ConTab[ 8 ] ],
  u8 = input[ ConTab[ 9 ] ],
  u9 = input[ ConTab[ 10 ] ];
  
  return
  u1 * LocalPos[1] * ( 2. * LocalPos[1] - 1. ) +
  4. * u2 * LocalPos[1] * LocalPos[3] +
  u3 * LocalPos[3] * ( 2. * LocalPos[3] - 1. ) +
  4. * u4 * LocalPos[2] * LocalPos[3] +
  u5 * LocalPos[2] * ( 2. * LocalPos[2] - 1. ) +
  4. * u6 * LocalPos[1] * LocalPos[2] +
  4. * u7 * LocalPos[1] * ( 1. - s ) +
  4. * u8 * LocalPos[3] * ( 1. - s ) +
  4. * u9 * LocalPos[2] * ( 1. - s ) +
  u10 * ( 1. - 2. * s ) * ( 1. - s );
}
/*------------------------------------------------------------------------------
 PURPOSE:
 Computes and stores inverse Jacobian and determinant for the new elements.
 
 INPUT:
 gauss: structure containing all gaussian quadrature information
 mesh->ConTab: connectivity table
 mesh->points: nodal co-ordinates
 local_dphi: derivatives of shape fcns in local co-ord system
 JacInv: inverse of Jacobian (empty at call time) * Jacobian determinant * gauss pt weight
 JacWeight: Jacobian magnitude times gauss pt weight (empty at call time)
 element: current new element # being processed.
 ConTab: connectivity table of newly formed elements
 GeoTab: nodal coordinates of new nodes
 
 OUTPUT:
 Filled JacInv and JacWeight, return true if the element volume is not negligible, false otherwise.
 
 By:
 Updated: 3/may/04 (Veeramani)
 ------------------------------------------------------------------------------*/
bool
GetJacobianPar(
const int       NbOfNodesPerElement,  ///< number of nodes per elements
const int       NbOfNodes,            ///< number of nodes in the mesh
const int       NbOfElements,         ///< number of elements in the mesh
      double**  GeoTab,               ///< geometric table
const gauss_t*  gauss,                ///< gauss quadrature
      double*   JacWeight,            ///< weighted jacobian
      int       element,              ///< index of an element
      int**     ParConTab,            ///< connectivity table of a particle
      double**  ParGeoTab )           ///< geometric table of a particle
{
  int ParElem = element - NbOfElements;
  
  for ( int gauss_pt = 1 ; gauss_pt <= gauss->NbOfPoints ; gauss_pt++ )
  {
    /* Calculate the Jacobian ; first zeroing all entries.
     Recall that first index in Jacobian cycles through xi, eta, zeta;
     second index cycles x, y, z.  Schematically, J is:
     dx/dxi   dy/dxi   dz/dxi
     dx/deta  dy/deta  dz/deta
     dx/dzeta dy/dzeta dz/dzeta     */
    
    double **Jacobian = NULL;
    AllocMdouble(3,3,Jacobian);
    
    for ( int local_node = 1 ; local_node <= NbOfNodesPerElement ; local_node++ )
    {
      int node = ParConTab[ local_node ][ ParElem ];
      
      double coordinate[4];

      if ( node <= NbOfNodes )
      {
        coordinate[1] = GeoTab[ node ][1];
        coordinate[2] = GeoTab[ node ][2];
        coordinate[3] = GeoTab[ node ][3];
      }
      else
      {
        coordinate[1] = ParGeoTab[1][ node - NbOfNodes ];
        coordinate[2] = ParGeoTab[2][ node - NbOfNodes ];
        coordinate[3] = ParGeoTab[3][ node - NbOfNodes ];
      }

      double *dphi = gauss->local_dphi[gauss_pt][ local_node ];
      
      for ( int xi_index = 1 ; xi_index <= 3 ; xi_index++ )
      for ( int x_index  = 1 ; x_index  <= 3 ; x_index++  )
        Jacobian[ xi_index ][x_index] += coordinate[x_index] * dphi[xi_index];
    }

    /*  Calculate the determinant of the Jacobian */
    double J_det = det3by3(Jacobian);
    FreeM(Jacobian);

    // check the value of the determinant, 
    if ( J_det <= JAC_DET_WARN )
    {
      assert_error( J_det > 0.,
        "Jacobian is negative = "DBL_FMT"\n",J_det);
      
      // the element is very small
      return false;
    }
    
    JacWeight[ gauss_pt ] = J_det * gauss->weight[ gauss_pt ];
  }
  
  return true;
}
//==============================================================================
/** Get weighted jacobian and/or the inverse jacobian.
 
 */
//==============================================================================
void
GetJacobianAtPoint(
const int       NbOfNodesPerElement,  ///< number of nodes per elements
      double**  GeoTab,               ///< global geometric table of the mesh
const int*      ElemCon,              ///< element connectivity table
      double**  local_dphi,           ///< local derivative of shape functions 
const double    weight,               ///< a weight
      double**  JacInv,               ///< inverse of jacobian
      double*   JacWeight )           ///< weighted jacobian
{  
  /* Calculate the Jacobian.
   Recall that first index in Jacobian cycles through xi, eta, zeta;
   second index cycles x, y, z. Schematically, J is
   
   dx/dxi    dy/dxi    dz/dxi
   dx/deta   dy/deta   dz/deta
   dx/dzeta  dy/dzeta  dz/dzeta
  */
  
  double **Jacobian = NULL;
  AllocMdouble(3, 3, Jacobian);
  
  for ( int local_node = 1; local_node <= NbOfNodesPerElement; local_node++ )
  {
    int node = ElemCon[local_node];
    double *coordinate = GeoTab[node],
         *dphi = local_dphi[ local_node ];
    
    for ( int xi_index = 1 ; xi_index <= 3 ; xi_index++ )
    for ( int x_index  = 1 ; x_index  <= 3 ; x_index++  )
      Jacobian[ xi_index ][x_index] += coordinate[x_index] * dphi[xi_index];
  }
  
  // get determinant
  double J_det = det3by3(Jacobian);
  
  // check its value
  assert_error(J_det > JAC_DET_WARN,
               "|Jacobian determinant| <= " DBL_FMT, JAC_DET_WARN  );
  
  /* Compute inverse times determinant times weight.
   The indexing key is shown below :
   
   dxi/dx   deta/dx   dzeta/dx
   dxi/dy   deta/dy   dzeta/dy
   dxi/dz   deta/dz   dzeta/dz
  */
  if ( JacInv ) inv3by3(Jacobian, JacInv);
  
  FreeM(Jacobian);
  
  // compute weighted jacobian
  *JacWeight = J_det * weight;
}
//==============================================================================
/** Compute the gradient at a point.
 */
//==============================================================================
void
GetGradientAtPoint(
const int  NbOfNodesPerElement, ///< number of nodes per elements
      double**  LocalDerPhi,    ///< local derivative at current point
      double**  JacInv,         ///< inverse of jacobian of current element
      double**  GlobalDerPhi )  ///< global derivative at current point
{  
  // Compute dphi/dx dphi/dy, dphi/dz at all Gauss points within this element
  for ( int local_node = 1 ; local_node <= NbOfNodesPerElement ; local_node++ )
    gemv3(JacInv, LocalDerPhi[ local_node ], GlobalDerPhi[ local_node ]);
}
