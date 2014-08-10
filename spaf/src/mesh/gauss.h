#ifndef GAUSS_H
#define GAUSS_H


typedef struct
{
  /// order of quadrature
  int order;
  
  int NbOfPoints; ///< number of points
  
  double
    **point,       ///< coordinates of the points
    *weight,       ///< weights at points
    ***local_dphi, ///< derivatives of shape functions at points for each direction
    **phi;         ///< shape functions at points
}
gauss_t;

gauss_t*
GaussCreate(
const int  NbOfNodesPerElement,
const int  QuadratureOrder );

void
FreeGauss(
gauss_t* gauss );

#endif
