/** This file is used to search points and elements of a mesh that are inside, at least partly, a particle.
 */

#include "includes.h"
#include "logging.h"
#include "linalg.h"
#include "mesh_search.h"
#include "memory.h"
#include "particle_mesh_find_elem.h"

//=============================================================================
/** Find whether a point is in a particle.

 We get the point coordinates in the frame of the particle and then check whether it is in the particle according to its kind.
 Return true if the point is in the particle, false otherwise.
 */
//=============================================================================
bool
IsPointInPar(
const particle_t* particle,     ///< pointer to a particle
const double        point[ 4 ] )  ///< point coordinates
{
  // copy point because it will be overwritten
  double PointInParFrame[4] = {0.,point[1],point[2],point[3]};

  // get point in particle frame
  ChangeFrame( particle->RotMat, particle->Pos0, PointInParFrame );
  
  // select the particle kind
  switch ( particle->kind )
  {
  case KIND_ELLIPSOID :
    return PointIsInEllipsoid( particle->ellipsoid, PointInParFrame );
    break;
  
  default :
    error( "In IsPointInPar : bad particle kind\n" );
  }
  
  error( "In IsPointInPar : bad particle kind\n" );
  return false;
}
//=============================================================================
/** Find whether a node is in a particle.
 
 Just call IsPointInPar with node coordinates.
 Return true if the point is in the particle, false otherwise.
 */
//=============================================================================
bool
IsNodeInPar(
const mesh_t*     mesh,
const particle_t* particle,
const int      node )
{
  return IsPointInPar( particle, mesh->points[ node ] );
}
//=============================================================================
/** Give the number of pressure nodes of an element that are inside a particle.
 
 Used to know if an element crosses a particle surface.
 */
//=============================================================================
static int
ElemNodeNbInPar(
const mesh_t*     mesh,     ///< mesh structure
const particle_t* particle, ///< pointer to a particle
const int      theElem ) ///< an element
{
  int numNodes = 0;

  for ( int i = 1 ; i <= PRE_NODES_PER_EL ; i++ )
    if ( IsNodeInPar( mesh, particle, mesh->ConTab[ theElem ][ i ] ) )
        numNodes++;

  return numNodes;
}
/*--------------------------------------------------------------------
Purpose:
  Get the elements and nodes on the velocity mesh that are touched by the particle.
  The elements may be partly inside, but only nodes inside are registered

Input:
  d: data structure
  particle: the particle
  element: current element checked
  NodeVisited: node visited is marked in the list of all nodes
  ElemVisited: element visited is marked in list of all elements

Output:
  ElemTouched: list of elements touched by the particles,
  ElemTouchedNodeIn: list of elements partly/fully inside the particles,
  nodesIS: list of nodes touched by the particles

By: Veeramani
Update: 29/Sep/06 Antoine
--------------------------------------------------------------------*/
static void
ParElemNodeTouchWrapped(
const mesh_t*     mesh,
const particle_t* particle,
const int      element,
      int**    ElemTouchedPtr,
      int**    ElemTouchedNodeInPtr,
      int**    NodeInPtr,
      bool*       NodeVisited,
      bool*       ElemVisited )
{
  int node, theNode, iElem, nextEl,
        NodeInPar,
        count; // Temporary variable for keeping count

  ElemVisited[ element ] = true;
  
  // get the number of nodes of element in particle
  NodeInPar = ElemNodeNbInPar( mesh, particle, element );

  if ( NodeInPar )
  {
    // increment the number of elements touched
    count = ++(*ElemTouchedPtr)[ 0 ];
    
    // Reallocate only when count exceeds a multiple of memory block
    if ( ( count % MBLOCK ) == 1 )
    {
      *ElemTouchedPtr = (int*) realloc( (*ElemTouchedPtr), ( count + MBLOCK ) * sizeof(int) );
      assert_error( *ElemTouchedPtr != NULL, "Memory allocation error");
    }
    
    (*ElemTouchedPtr)[ count ] = element;

    count = ++(*ElemTouchedNodeInPtr)[ 0 ];
    
    // Reallocate only when count exceeds a multiple of memory block
    if ( ( count % MBLOCK ) == 1 )
    {
      *ElemTouchedNodeInPtr = (int*) realloc( (*ElemTouchedNodeInPtr), ( count + MBLOCK ) * sizeof(int) );
      assert_error( *ElemTouchedNodeInPtr != NULL, "Memory allocation error");
    }
    (*ElemTouchedNodeInPtr)[ count ] = NodeInPar;

    // check neighbouring elements recursively
    for ( node = 1 ; node <= NODES_PER_EL ; node++ )
    {
      theNode = mesh->ConTab[ element ][ node ];
      
      if ( NodeVisited[ theNode ] == false )
      {
        NodeVisited[ theNode ] = true;

        // Check whether node is inside the particle
        if ( IsNodeInPar( mesh, particle, theNode ) )
        {
          count = ++(*NodeInPtr)[ 0 ];
          
          // Reallocate only when count exceeds a multiple of memory block
          if ( ( count % MBLOCK ) == 1 )
          {
            *NodeInPtr = ( int * ) realloc( (*NodeInPtr), ( count + MBLOCK ) * sizeof(int) );
            assert_error( *NodeInPtr != NULL, "Memory allocation error");
          }
          
          (*NodeInPtr)[ count ] = theNode;
        }

        // check in inverse ConTab
        for ( iElem = 1 ; iElem <= NbOfEntriesInRow(mesh->InvConTab, theNode-1);  iElem++ )
        {
          nextEl = EntryNode(mesh->InvConTab, theNode-1, iElem-1);
          
          if ( ElemVisited[ nextEl ] == false )
            ParElemNodeTouchWrapped(
              mesh, particle, nextEl,
              ElemTouchedPtr, ElemTouchedNodeInPtr, NodeInPtr, NodeVisited, ElemVisited );
        }
      }
    }
  }
}

void
ParElemNodeTouch(
const mesh_t*       mesh,
const particle_t*   particle,
const int        element,
      int**      ElemTouchedPtr,
      int**      ElemTouchedNodeInPtr,
      int**      NodesIn )
{
  bool *NodeVisited = NULL, // nodes and
       *ElemVisited = NULL; // elements visited
  
  AllocVbool( mesh->NbOfNodes,NodeVisited );
  AllocVbool( mesh->NbOfElements,ElemVisited );

  for ( int i = 0 ; i <= mesh->NbOfElements ; i++ )
    ElemVisited[i] = false;
    
  for ( int i = 0 ; i <= mesh->NbOfNodes ; i++ )
    NodeVisited[i] = false;

  (*NodesIn) = (int*) realloc( (*NodesIn), sizeof(int) );
  
  assert_error( (*NodesIn) != NULL, "Memory allocation error");

  (*NodesIn)[ 0 ] = 0; // Number of nodes inside particle

  ParElemNodeTouchWrapped(
    mesh, particle, element, ElemTouchedPtr, ElemTouchedNodeInPtr,
    NodesIn, NodeVisited, ElemVisited );

  debug( "\t" INT_FMT " nodes inside\n", (*NodesIn)[ 0 ] );

  free(ElemVisited);
  free(NodeVisited);
}
