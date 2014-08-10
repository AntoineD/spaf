#ifndef MESH_H
#define MESH_H

#include "planes.h"
#include "gauss.h"
#include "nns.h"
#include "sparse_matrix.h"

/// Nodes per element, velocity nodes
#define NODES_PER_EL 10 
/// Pressure nodes per element
#define PRE_NODES_PER_EL 4

/// defines an io for checking characteristic that crosses it
typedef struct
{
  char type;
  double radius, center[3], normal[3];
  equation_t plane;
}
io_t;

/// defines a mesh structure
typedef struct
{
  int
    NbOfNodes,              ///< number of nodes (velocity)
    NbOfPressureNodes,      ///< number of pressure nodes
    NbOfFreeVelocityNodes,  ///< number of free velocity nodes
    NbOfFreePressureNodes,  ///< number of free pressure nodes
    NbOfElements,           ///< number of elements
    **ConTab,               ///< Connectivity table
    **ConTabLocal,          ///< Local connectivity within an element
    *VelToPreNodeMap;       ///< map between velocity and pressure nodes
  
  /// list of inner domain nodes that are outside the outer domain (version z only)
  bool *OutsideNodes;
  
  /// Inverse Connectivity table
  SparseMatrix_t* InvConTab;
  
  /// coordinates of the nodes
  double **points;

  /// Geometry i/o parameter for checking if a characteristic exits the computational domain
  int NbOfOpenings;
  io_t io[3];
  
  /// gauss quadrature structure for integration
  gauss_t* gauss;
  
  /// nearest neighbor search structure
  nns_t* nns;
  
  /// nodes Courant number (min internodal distance around each point)
  double *courant;

}
mesh_t;

mesh_t*
ReadMesh(
const char* FileName );

SparseMatrix_t*
GetInverseConnectivity(
const int    NbOfNodes,
const int    NbOfElements,
const int    NbOfNodesPerElement,
      int**  ConTab );

#endif
