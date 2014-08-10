#ifndef SHAPE_FUNCTIONS_H
#define SHAPE_FUNCTIONS_H

#include "gauss.h"

void
GetShapeFunctionsAtPoint(
const int  order,
const double    point[3],
      double*   PhiAtPoint );

void
GetDerivativeOfShapeFunctionsAtPoint(
const int  order,
const double    point[3],
      double**  DerPhi );

double
InterpolateOrder2(
const int  ConTab[11],
const double    LocalPos[4],
const double*   input );

void
GetJacobianAtPoint(
                   const int  NbOfNodesPerElement,
                   double**  GeoTab,                   
                   const int* ElemCon,            
                   double**  local_dphi,                
                   const double    weight,             
                   double**  JacInv,                   
                   double*   JacWeight );              

bool
GetJacobianPar(
               const int    NbOfNodesPerElement,
               const int    NbOfNodes,
               const int    NbOfElements,
               double**    GeoTab,
               const gauss_t*  Gauss,
               double*     J_times_wt,
               int    element,
               int**  ParConTab,
               double**    ParGeoTab );

void
GetGradientAtPoint(
                    const int  NbOfNodesPerElement,
                   double**  LocalDerPhi,
                    double** JacInv,
                   double** GlobalDerPhi);

#endif
