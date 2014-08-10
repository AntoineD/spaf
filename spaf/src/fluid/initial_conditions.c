/** \file
  Read initial conditions of the fluid.
*/

#include "logging.h"
#include "strings.h"
#include "mesh_search.h"
#include "output.h"
#include "convection.h"
#include "file_system.h"
#include "shape_functions.h"
#include "initial_conditions.h"
#include "version_z.h"
#include "fluid_io.h"

//==============================================================================
/** Interpolate the pressure at a point.
 Given the local position of the point in the current element, computes the first order interpolation of a field at this point.
 Returns the value of the field at the point.
 */
//==============================================================================
static double
PreInterpolate(
const mesh_t* mesh,         ///< mesh structure
const int     ElemId,       ///< index of an element
const double  LocalPos[4],  ///< local position of a point in the element
const double* Pre )         ///< a field
{
  return
  Pre[ mesh->VelToPreNodeMap[ mesh->ConTab[ ElemId ][ 1 ] ] ] * LocalPos[1] +
  Pre[ mesh->VelToPreNodeMap[ mesh->ConTab[ ElemId ][ 2 ] ] ] * LocalPos[3] +
  Pre[ mesh->VelToPreNodeMap[ mesh->ConTab[ ElemId ][ 3 ] ] ] * LocalPos[2] +
  Pre[ mesh->VelToPreNodeMap[ mesh->ConTab[ ElemId ][ 4 ] ] ] *
  ( 1. - LocalPos[1] - LocalPos[3] - LocalPos[2] );
}
//=============================================================================
/** Read fluid initial conditions from disc.
 
  Initial conditions must be written as a field.
 */
//=============================================================================
void
ReadFluidInitialConditions(
const char*     PathToCase,  ///< path to mesh directory
const mesh_t*   mesh,        ///< mesh structure
      double**  Vel,         ///< velocity field
      double*   Pre )        ///< pressure field
{
  PrintTitle("Reading fluid initial conditions");
  
  // get the path to the mesh directory inside the case one
  char *PathToIC = StringConcat(PathToCase, "/initial_conditions");
  
  info("from %s\n", PathToIC );
  
  // read the field
  ReadField( PathToIC, mesh->NbOfPressureNodes, Pre, mesh->NbOfNodes, Vel);
  
  free(PathToIC);
}
//=============================================================================
/** Compute initial conditions by interpolation from outer domain, version Z only.
 
  FluidVar->VelOld1 and FluidVar->Pre are filled here, pressure is shifted in such a way it has homogeneous Dirichlet bc.
*/
//=============================================================================
void
GetFluidInitialConditions_z(
const double      position[4],///< particle position
const mesh_t*     mesh,       ///< inner mesh structure
const mesh_t*     mesh_o,     ///< inner mesh structure
const FluidVar_t* FluidVar_o, ///< outer fluid structure
      FluidVar_t* FluidVar )  ///< inner fluid structure
{
  info("Applying fluid initial conditions : ");

  // loop over all inner nodes
  for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
  {
    double
    *point = mesh->points[NodeId],
    NodeCoord[4] = {0., position[1]+point[1], position[2]+point[2], position[3]+point[3]}, // node coordinates in mesh_o frame
    node_LocalCoord_o[4] = {0.,0.,0.,0.};  // local node coordinates in outer domain
    
    // find outer element containing inner node, we do full bb search
    int temp = 0,  // not used, but required by SearchElemContainingPoint
          ElemContainer_o = SearchElemContainingPoint(
                            mesh_o, true, NodeCoord, &temp, &temp );
        
    // check that the node is not outside outer domain
    assert_error( ElemContainer_o != 0,
      "node %u is not in global mesh", NodeId );

    // get the local position of the inner node in the outer element
    bool ElemContainerFound = MeshGetLocalPosition(
      mesh_o, ElemContainer_o, NodeCoord, node_LocalCoord_o );

    assert_error( ElemContainerFound == true ,
             "node %u is not in the element %u", NodeId, ElemContainer_o );

    // interpolate the velocity from the outer element at the inner node
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      FluidVar->VelOld1[dir][NodeId] = InterpolateOrder2(
        mesh_o->ConTab[ ElemContainer_o ], node_LocalCoord_o, FluidVar_o->Vel[dir]);
    
    // interpolate the pressure from the outer element at the inner node
    int PreNode = mesh->VelToPreNodeMap[ NodeId ]; // pressure node index
    
    // is it a pressure node ?
    if ( PreNode != 0 )
      FluidVar->Pre[ PreNode ] = PreInterpolate(
        mesh_o, ElemContainer_o, node_LocalCoord_o, FluidVar_o->Pre );
  }

  // find a prescribed pressure node
  int PrescribedNodeId = 0;
  for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
  {
    int PreNodeId = mesh->VelToPreNodeMap[ NodeId ];
    
    if ( PreNodeId > mesh->NbOfFreePressureNodes )
    {
      PrescribedNodeId = NodeId;
      break;
    }
  }
  
  double PreValue = FluidVar->Pre[mesh->VelToPreNodeMap[PrescribedNodeId]];
  
  // loop over all inner nodes and substract pressure value
  for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
  {
    int PreNode = mesh->VelToPreNodeMap[ NodeId ]; // pressure node index
    
    if ( PreNode != 0 )
      FluidVar->Pre[ PreNode ] -= PreValue;        
  }
}
