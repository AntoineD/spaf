/** \file
 Gauss quadratures.
 
 Create and compute the points and weights for gauss quadrature, as well as values of the local shape functions and their derivatives at those points.
 */

#include "includes.h"
#include "logging.h"
#include "memory.h"
#include "shape_functions.h"
#include "gauss.h"

//==============================================================================
/** Quadrature rule of degree 5 with 15 points.
 */
//==============================================================================
static void
TetrahedronOrder5(
int*      NbOfPoints, ///< number of gauss points
double**  weight,     ///< weights at those points
double*** point )     ///< coordinates of the points
{
  double w[4] = {
    6.02678571428571428571428571428571e-3,
    0.0302836780970891758063769725577305,
    0.0116452490860289694108936091443380,
    0.0109491415613864593456430191124068};
  
  double coord[7];
  
  coord[1] = 0.;
  coord[5] = 0.0665501535736642982398804642263025;
  coord[3] = 1. / 11.;
  coord[2] = 0.25;
  coord[0] = 1. / 3.;
  coord[6] = 0.433449846426335701760119535773697;
  coord[4] = 8. / 11.;
  
  *NbOfPoints = 15;
  
  AllocVdouble((*NbOfPoints), (*weight));
  AllocMdouble((*NbOfPoints), 3-1, (*point));
    
  (*weight)[1]  = w[0];
  (*weight)[2]  = w[0];
  (*weight)[3]  = w[0];
  (*weight)[4]  = w[0];
  (*weight)[5]  = w[1];
  (*weight)[6]  = w[2];
  (*weight)[7]  = w[2];
  (*weight)[8]  = w[2];
  (*weight)[9]  = w[2];
  (*weight)[10] = w[3];
  (*weight)[11] = w[3];
  (*weight)[12] = w[3];
  (*weight)[13] = w[3];
  (*weight)[14] = w[3];
  (*weight)[15] = w[3];

  (*point)[1] [0] = coord[0];
  (*point)[2] [0] = coord[1];
  (*point)[3] [0] = coord[0];
  (*point)[4] [0] = coord[0];
  (*point)[5] [0] = coord[2];
  (*point)[6] [0] = coord[3];
  (*point)[7] [0] = coord[4];
  (*point)[8] [0] = coord[3];
  (*point)[9] [0] = coord[3];
  (*point)[10][0] = coord[5];
  (*point)[11][0] = coord[6];
  (*point)[12][0] = coord[5];
  (*point)[13][0] = coord[5];
  (*point)[14][0] = coord[6];
  (*point)[15][0] = coord[6];

  (*point)[1] [1] = coord[0];
  (*point)[2] [1] = coord[0];
  (*point)[3] [1] = coord[1];
  (*point)[4] [1] = coord[0];
  (*point)[5] [1] = coord[2];
  (*point)[6] [1] = coord[3];
  (*point)[7] [1] = coord[3];
  (*point)[8] [1] = coord[4];
  (*point)[9] [1] = coord[3];
  (*point)[10][1] = coord[6];
  (*point)[11][1] = coord[5];
  (*point)[12][1] = coord[6];
  (*point)[13][1] = coord[5];
  (*point)[14][1] = coord[5];
  (*point)[15][1] = coord[6];

  (*point)[1] [2] = coord[0];
  (*point)[2] [2] = coord[0];
  (*point)[3] [2] = coord[0];
  (*point)[4] [2] = coord[1];
  (*point)[5] [2] = coord[2];
  (*point)[6] [2] = coord[3];
  (*point)[7] [2] = coord[3];
  (*point)[8] [2] = coord[3];
  (*point)[9] [2] = coord[4];
  (*point)[10][2] = coord[6];
  (*point)[11][2] = coord[5];
  (*point)[12][2] = coord[5];
  (*point)[13][2] = coord[6];
  (*point)[14][2] = coord[6];
  (*point)[15][2] = coord[5];
}
//==============================================================================
/** Quadrature rule of degree 1 with 4 points.
 */
//==============================================================================
static void
TetrahedronOrder1(
int*      NbOfPoints, ///< number of gauss points
double**  weight,     ///< weights at those points
double*** point )     ///< coordinates of the points
{
  *NbOfPoints = 4;
  
  AllocVdouble((*NbOfPoints),(*weight));
  AllocMdouble((*NbOfPoints),3-1,(*point));
  
  for ( int i = 1 ; i <= *NbOfPoints ; i++ )
    (*weight)[i] = 1. / 24.;
  
  (*point)[1][0] = 1.;
  (*point)[1][1] = 0.;
  (*point)[1][2] = 0.;
  
  (*point)[2][0] = 0.;
  (*point)[2][1] = 0.;
  (*point)[2][2] = 1.;
  
  (*point)[3][0] = 0.;
  (*point)[3][1] = 1.;
  (*point)[3][2] = 0.;
  
  (*point)[4][0] = 0.;
  (*point)[4][1] = 0.;
  (*point)[4][2] = 0.;
}
//==============================================================================
/** Quadrature rule of degree 3 with 5 points.
 */
//==============================================================================
static void
TetrahedronOrder3(
int*      NbOfPoints, ///< number of gauss points
double**  weight,     ///< weights at those points
double*** point )     ///< coordinates of the points
{
  *NbOfPoints = 5;
  
  AllocVdouble((*NbOfPoints),(*weight));
  AllocMdouble((*NbOfPoints),3-1,(*point));
  
  (*weight)[1] = -0.13333333333333333333333;
  
  (*point)[1][0] = 0.25;
  (*point)[1][1] = 0.25;
  (*point)[1][2] = 0.25;

  for ( int i = 2 ; i <= *NbOfPoints ; i++ )
    (*weight)[i] = 0.075;
          
  (*point)[2][0] = 0.166666666666666666666666;
  (*point)[2][1] = 0.166666666666666666666666;
  (*point)[2][2] = 0.166666666666666666666666;
          
  (*point)[3][0] = 0.166666666666666666666666;
  (*point)[3][1] = 0.166666666666666666666666;
  (*point)[3][2] = 0.5;
          
  (*point)[4][0] = 0.166666666666666666666666;
  (*point)[4][1] = 0.5;
  (*point)[4][2] = 0.166666666666666666666666;
          
  (*point)[5][0] = 0.5;
  (*point)[5][1] = 0.166666666666666666666666;
  (*point)[5][2] = 0.166666666666666666666666;
}
//==============================================================================
/** Create and fill a gauss quadrature structure.
 */
//==============================================================================
gauss_t*
GaussCreate(
const int NbOfNodesPerElement,  ///< number of nodes per element
const int order )               ///< order of the quadrature
{
  gauss_t* gauss = (gauss_t*) calloc( 1, sizeof(gauss_t) );
  
  gauss->order = order;

  // get gaussian quadrature data
  if ( gauss->order == 1 )
    TetrahedronOrder1( &gauss->NbOfPoints, &gauss->weight, &gauss->point);
  
  else if ( gauss->order == 3 )
    TetrahedronOrder3( &gauss->NbOfPoints, &gauss->weight, &gauss->point);
  
  else if ( gauss->order == 5 )
    TetrahedronOrder5( &gauss->NbOfPoints, &gauss->weight, &gauss->point);

  else
    error("Bad order");

  // Allocate shape functions stuff
  gauss->local_dphi = AllocVM_double(gauss->NbOfPoints, NbOfNodesPerElement, 3);
  
  AllocMdouble( gauss->NbOfPoints, NbOfNodesPerElement, gauss->phi );
  
  // get shape function stuff at all points
  for ( int gauss_pt = 1 ; gauss_pt <= gauss->NbOfPoints ; gauss_pt++ )
  {
    double
    *point   = gauss->point[gauss_pt],
    *phi     = gauss->phi[ gauss_pt ],
    **phider = gauss->local_dphi[ gauss_pt ]; 
        
    GetShapeFunctionsAtPoint( gauss->order, point, &phi[1] );
    GetDerivativeOfShapeFunctionsAtPoint( gauss->order, point, phider );
  }
    
  return gauss;
}
//==============================================================================
/** Free a gauss quadrature structure memory.
 */
//==============================================================================
void
FreeGauss(
gauss_t* gauss )  ///< gauss structure
{
  free(gauss->weight);
  FreeM(gauss->point);
  FreeM(gauss->phi);
  FreeVM_double(gauss->local_dphi,gauss->NbOfPoints);
  free(gauss);
}
