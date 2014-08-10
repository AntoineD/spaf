/** \file
 Manage particle mesh.
 */

#include "includes.h"
#include "memory.h"
#include "linalg.h"
#include "logging.h"
#include "statistics.h"
#include "particle.h"
#include "particle_mesh_subdivide.h"
#include "particle_mesh_find_elem.h"
#include "particle_io.h"
#include "mesh_search.h"
#include "shape_functions.h"
#include "particle_mesh.h"

//=============================================================================
/** Sort elements touched by the particle into full or partial elements.
*/
//=============================================================================
static void
SortElements(
int*  elemT,
int*  elemP,
int** ElemFullPtr,
int** ElemPartlyPtr )
{
  int elemCount = 0; // Counter for number of new elements

  // Get the original elements from elemT and elemP and distribute
  for ( int iElem = 1 ; iElem <= elemT[ 0 ] ; iElem++ )
  {
    // Check whether it is fully or partly inside
    if ( elemP[ iElem ] == PRE_NODES_PER_EL )
    {
      // Element is fully inside, register in ElemFull
      elemCount = ++(*ElemFullPtr)[ 0 ];

      if ( ( elemCount % MBLOCK ) == 1 )
      {
        *ElemFullPtr = (int*) realloc( (*ElemFullPtr), ( elemCount + MBLOCK ) * sizeof(int) );
         assert_error( *ElemFullPtr != NULL, "Memory allocation error");
      }
      
      (*ElemFullPtr)[ elemCount ] = elemT[ iElem ];
    }
    else
    {
      // Element partly inside, register in ElemPartly
      elemCount = ++(*ElemPartlyPtr)[ 0 ];

      if ( ( elemCount % MBLOCK ) == 1 )
      {
       (*ElemPartlyPtr) = (int*) realloc( (*ElemPartlyPtr), ( elemCount + MBLOCK ) * sizeof(int) );
       assert_error( *ElemPartlyPtr != NULL, "Memory allocation error");
      }
           
      (*ElemPartlyPtr)[ elemCount ] = elemT[ iElem ];
    }
  }
}
//=============================================================================
/** Initialize particle mesh.
 
 Reallocate all arrays such that we are in the same memory state as when this
 structure was created. Each array has one element which is 0.
 */
//=============================================================================
static void
InitializeMeshStructure(
Parmesh_t* ParMesh )  ///< particle mesh structure
{
  ParMesh->ElemFull = (int*) realloc( ParMesh->ElemFull, 1 * sizeof(int) );
  assert_error( ParMesh->ElemFull != NULL, "Memory allocation error");
  ParMesh->ElemFull[0] = 0;
    
  ParMesh->ElemPartly= (int*) realloc( ParMesh->ElemPartly, 1 * sizeof(int) );
  assert_error( ParMesh->ElemPartly != NULL, "Memory allocation error");
  ParMesh->ElemPartly[0] = 0;

  ParMesh->points[0][0] = 0.;
  
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    ParMesh->points[dir] = (double*)
      realloc(ParMesh->points[dir],1*sizeof(double));
    assert_error( ParMesh->points[dir] != NULL, "Memory allocation error");
  }
  
  ParMesh->ConTab[0][0] = 0;

  for ( int i = 1 ; i <= 10 ; i++ )
  {
    ParMesh->ConTab[ i ] = (int*)
      realloc( ParMesh->ConTab[ i ], 1*sizeof(int) );
    assert_error( ParMesh->ConTab[ i ] != NULL, "Memory allocation error");
  }
}
//=============================================================================
/** Initialize particle mesh.
 
 Reallocate all arrays such that we are in the same memory state as when this
 structure was created. Each array has one element which is 0
 */
//=============================================================================
void
GetParticleMesh(
const mesh_t*     mesh,
      particle_t  particle[] )
{
  if ( particle == NULL ) return;
  
  debug( "\nParticles : meshing\n");
  
  // Elements touched by the particles
  int *ElemTouched = NULL,
         *ElemTouchedNodesIn = NULL,
         // store the surface nodes inside particle domain
         *NodesConnectedToOutside = NULL,
         // store the surface innode-outnode connectivity
         ***InOutSurfConTab = NULL;
  
  AllocVint( NODELMT, NodesConnectedToOutside );
  
  InOutSurfConTab = (int***) calloc( mesh->NbOfNodes + 1, sizeof(int**) );
  assert_error( InOutSurfConTab != NULL, "Memory allocation error");
  
  for ( int node = 1 ; node <= mesh->NbOfNodes ; node++ )
  {
    InOutSurfConTab[ node ] = (int**) calloc( 3, sizeof(int*) );
    assert_error( InOutSurfConTab[ node ] != NULL, "Memory allocation error");
    
    for ( int i = 0 ; i <= 2 ; i++ )
      AllocVint(0,InOutSurfConTab[ node ][ i ]);
  }

  // fill in the ParMesh structure
  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    // pointer to particle mesh
    Parmesh_t *ParMesh = particle[iPar].mesh;
    
    // Initialization
    InitializeMeshStructure( ParMesh );

    // Number of surface nodes in particle domain
    NodesConnectedToOutside[ 0 ] = 0;
  
    ElemTouched = (int*) realloc( ElemTouched, sizeof(int) );
    ElemTouched[ 0 ] = 0; // Number of element touched
    
    ElemTouchedNodesIn = (int*) realloc( ElemTouchedNodesIn, sizeof(int) );
    ElemTouchedNodesIn[ 0 ] = 0; // Number of element partly inside

    // Find the element which contains the center of mass of the particle
    int temp = 0, // for unused variables

    ElemMassCenter = SearchElemContainingPoint(
                      mesh, true, particle[iPar].Pos0, &temp, &temp);

    // Find the elements and nodes inside the current particle
    ParElemNodeTouch(
      mesh, &particle[ iPar ], ElemMassCenter,
      &ElemTouched, &ElemTouchedNodesIn, &(ParMesh->NodesIn) );

    // Distributes elements touched into elements fully and partly inside
    SortElements(
      ElemTouched, ElemTouchedNodesIn,
      &(ParMesh->ElemFull), &(ParMesh->ElemPartly) );

    debug( "\t" INT_FMT " elements fully inside" "\n",
          ParMesh->ElemFull[ 0 ] );
    debug( "\t" INT_FMT " elements partly inside" "\n",
          ParMesh->ElemPartly[ 0 ] );

    // Find the node-node connections cut
    // find edge node and fill InOutSurfConTab structure
    FillInOutSurfConTab(
      mesh, &particle[ iPar ],
      ParMesh->ElemPartly, ParMesh->points,
      InOutSurfConTab, NodesConnectedToOutside );
    
    // Surface elements are processed on the basis of InOutSurfConTab
    // and added to ElemFull
    ParConnGeo(
      mesh, &particle[ iPar ], ParMesh->points, ParMesh->ConTab,
      &(ParMesh->ElemFull), ParMesh->ElemPartly,
      InOutSurfConTab, NodesConnectedToOutside );
    
    debug( "\t" INT_FMT " nodes on surface\n",
          NodesConnectedToOutside[ 0 ] );
    debug( "\t" INT_FMT " nodes created\n",
          (int) ParMesh->points[ 0 ][ 0 ] );
    debug( "\t" INT_FMT " elements inside after sub-meshing\n",
          ParMesh->ElemFull[ 0 ] );
    
    // write particle mesh to file
//    WriteParMeshVTK( iPar, mesh, ParMesh );

    // statistics
    StatAdd("Particle meshing : elements fully inside",
            ParMesh->ElemFull[ 0 ]);
  }
  
  free(ElemTouched);
  free(ElemTouchedNodesIn);

  free(NodesConnectedToOutside);
  
  for ( int node = 1 ; node <= mesh->NbOfNodes ; node++ )
  {
    for ( int i = 0 ; i <= 2 ; i++ )
      free(InOutSurfConTab[ node ][i]);
    
    free(InOutSurfConTab[ node ]);
  }
  
  free(InOutSurfConTab);
}
//=============================================================================
/// 
//=============================================================================
bool
GetParElementJacobian(
const mesh_t*     mesh,
const Parmesh_t*  ParMesh,
const int      element,
      double*       J_times_wt )
{
  return GetJacobianPar( NODES_PER_EL,
    mesh->NbOfNodes, mesh->NbOfElements, mesh->points,
    mesh->gauss, J_times_wt,
    element, ParMesh->ConTab, ParMesh->points );
}

