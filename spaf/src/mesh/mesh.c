/** \file
  Create a mesh structure, compute inverse connectivity table.
  */

#include "includes.h"
#include "memory.h"
#include "data_io.h"
#include "logging.h"
#include "output.h"
#include "strings.h"
#include "parse.h"
#include "collision.h"
#include "shape_functions.h"
#include "file_system.h"
#include "sparse_pattern.h"
#include "mesh_search.h"
#include "mesh.h"

//==============================================================================
/// Read openings data.
//==============================================================================
static void
ReadOpenings(
const char*   FileName, ///< name of the file with openings data
      mesh_t* mesh )    ///< mesh structure
{
  FILE* FileId = fopen(FileName, "r");
  
  // return if no openings file
  if ( FileId == NULL ) return;
  
  // read sections
  while ( TokenizeFileLine(FileId) != 0 )
  {    
    char type = ' ';
    
    if ( SectionIs("inlet") )
      type = 'i';
    else if ( SectionIs("outlet") )
      type = 'o';
    else
      error( "Bad section, should be inlet or outlet");
    
    mesh->NbOfOpenings++;
    assert_error( mesh->NbOfOpenings <= 3,"Increase NbOfOpenings");
    int i = mesh->NbOfOpenings - 1;

    mesh->io[i].type = type;

    // read radius
    TokenizeFileLine(FileId);
    assert_error( TokenNameIs("radius"), "keyword should be radius");
    mesh->io[i].radius = Token2double(3,'[',-DBL_MAX,DBL_MAX,']');
    
    // read center
    TokenizeFileLine(FileId);
    assert_error( TokenNameIs("center"), "keyword should be center");
    for ( int dir = 0 ; dir < 3 ; dir++ )
      mesh->io[i].center[dir] = Token2double(dir+3,'[',-DBL_MAX,DBL_MAX,']');
    
    // read normal
    TokenizeFileLine(FileId);
    assert_error( TokenNameIs("normal"), "keyword should be normal");
    for ( int dir = 0 ; dir < 3 ; dir++ )
      mesh->io[i].normal[dir] = Token2double(dir+3,'[',-DBL_MAX,DBL_MAX,']');
    
    // check end of current section
    TokenizeFileLine(FileId);
    assert_error( SectionIsEnd(), "No end section");
  };
  
  // print infos
  for ( int i = 0 ; i < mesh->NbOfOpenings ; i++ )
  {
    if ( mesh->io[i].type == 'i' )
      info("Inlet "INT_FMT" :""\n",i);
    else
      info("Outlet "INT_FMT" :""\n",i);
    
    info("\t""radius = "DBL_FMT"\n",mesh->io[i].radius);
    
    info("\t""center = ");
    for ( int dir = 0 ; dir < 3 ; dir++ )
      info("%g ",mesh->io[i].center[dir]);
    
    info("\n""\t""normal = ");
    for ( int dir = 0 ; dir < 3 ; dir++ )
      info("%g ",mesh->io[i].normal[dir]);
    info("\n");
  }
  
  // fill plane data
  for ( int i = 0 ; i < mesh->NbOfOpenings ; i++ )
  {
    io_t* io = &mesh->io[i];
    
    double
    A = io->normal[0],
    B = io->normal[1],
    C = io->normal[2],
    D = - A * io->center[0] - B * io->center[1] - C * io->center[2];
    
    io->plane = (equation_t){{D,A,B,C}};      
  }
}
//==============================================================================
/// Create an inverse connectivity table.
//==============================================================================
SparseMatrix_t*
GetInverseConnectivity(
const int   NbOfNodes,            ///< number of nodes
const int   NbOfElements,         ///< number of elements
const int   NbOfNodesPerElement,  ///< number of nodes per elements
      int** ConTab )              ///< connectivity table
{  
  int *NbOfElemPerNode = NULL;
  AllocVint( NbOfNodes,NbOfElemPerNode );
  
  // register the number of elements a node is in
  for ( int element = 1; element <= NbOfElements; element++ )
  for ( int LocalNodeId = 1; LocalNodeId <= NbOfNodesPerElement; LocalNodeId++ )
  {
    int node = ConTab[element][LocalNodeId];
    NbOfElemPerNode[ node ] ++;
  }
  
  int NbOfEntries = 0;
  
  for ( int NodeId = 1 ; NodeId <= NbOfNodes ; NodeId++ )
    NbOfEntries += NbOfElemPerNode[ NodeId ] ;
  
  int *offset = NULL;
  AllocVint(NbOfNodes, offset);
  
  // Allocate memory for inverse of connectivity table
  int *entries = NULL;
  AllocVint(NbOfEntries-1, entries);
  
  // Fill in the offsets and copy value into the global structure
  offset[0] = 0;
  
  for ( int iNode = 1; iNode <= NbOfNodes; iNode++ )
  {
    offset[iNode] = offset[iNode-1] + NbOfElemPerNode[ iNode ];
    NbOfElemPerNode[ iNode ] = 0;
  }
  
  for ( int iElem = 1; iElem <= NbOfElements; iElem++ )
  for ( int j = 1; j <= NbOfNodesPerElement; j++ )
  {
    int node = ConTab[ iElem ][ j ];
    NbOfElemPerNode[ node ]++;
    entries[ offset[ node-1 ] + NbOfElemPerNode[ node ] - 1 ] = iElem; // ?
  }
  
  free( NbOfElemPerNode );

  info("Inverse connectivity table :" "\n");
  SparseMatrix_t* InvConTab = CreateSparseMatrix(
    false, NbOfNodes, offset, NULL, NULL, NULL, NULL );
  
  InvConTab->entries_n = entries;
  
  return InvConTab;
}
//==============================================================================
/** Compute the minimum node to node distance for all nodes.
 */
//==============================================================================
static double*
GetCourant(
const int             NbOfNodes,  ///< number of nodes
      int**           ConTab,     ///< connectivity
const SparseMatrix_t* InvConTab,  ///< inverse connectivity
      double**        points )    ///< nodes coordinates
{
  double* courant = NULL;
  AllocVdouble(NbOfNodes,courant);
  
  for ( int node = 1 ; node <= NbOfNodes ; node++ )
  {
    double *point = points[ node ];
    double dist_min = DBL_MAX;
    
    // go over elements which contain the node
    for ( int j = 1 ; j <= NbOfEntriesInRow(InvConTab, node-1) ; j++ )
    {
      int element = EntryNode(InvConTab, node-1, j-1);
      
      // go over nodes of this element
      for ( int local_node = 1 ; local_node <= NODES_PER_EL ; local_node++ )
      {
        int neighbor = ConTab[ element ][ local_node ];
        
        // compute node to node distance
        if ( neighbor != node )
        {
          double* neighbor_point = points[ neighbor ];
          
          double dist = sqrt( pow( neighbor_point[ 1 ] - point[ 1 ] ,2) +
                             pow( neighbor_point[ 2 ] - point[ 2 ] ,2) +
                             pow( neighbor_point[ 3 ] - point[ 3 ] ,2) );
          
          if ( dist < dist_min ) dist_min = dist;
        }
      }
    }
    
    // for each point load the minimum distance to any other point
    courant[ node ] = dist_min;
  }
  
  return courant;
}
//==============================================================================
/** Read the mesh and geometric data.
 
 This function reads the mesh, openings data if any and planes data if any.
*/
//==============================================================================
mesh_t*
ReadMesh(
const char* PathToCase ) ///< path to case directory
{
  PrintTitle("Reading mesh");
  
  // get the path to the mesh directory inside the case one
  char *PathToMesh = StringConcat(PathToCase, "/mesh");

  info("from %s\n", PathToMesh );
  
  // move to mesh directory
  char *WorkingDirectory = GoToDirectory(PathToMesh);  

  if ( WorkingDirectory == NULL ) error("Cannot open mesh""\n");
  
  free(PathToMesh);
  
  mesh_t* mesh = (mesh_t*) calloc( 1, sizeof(mesh_t) );
  
  //---------------------------------------------------------------------------
  // read points
  //---------------------------------------------------------------------------
  int number = 0;
 
  double *array_double = NULL;

  ReadArray_double( "points", &number, &array_double );
  
  mesh->NbOfNodes = number / 3;
  AllocMdouble( mesh->NbOfNodes, 3, mesh->points );

  int id = 0;
  for ( int node = 1 ; node <= mesh->NbOfNodes ; node++ )
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    mesh->points[node][dir] = array_double[id++];
  
  free(array_double);
  
  //---------------------------------------------------------------------------
  // read connectivity
  //---------------------------------------------------------------------------
  int *array_int = NULL;

  ReadArray_int( "connectivity", &number, &array_int );
  
  mesh->NbOfElements = number / NODES_PER_EL;
  AllocMint( mesh->NbOfElements, NODES_PER_EL, mesh->ConTab );
  
  id = 0;
  for ( int element = 1 ; element <= mesh->NbOfElements ; element++ )
  for ( int node = 1 ; node <= NODES_PER_EL ; node++ )
    mesh->ConTab[element][node] = array_int[id++];
  
  free(array_int);

  //---------------------------------------------------------------------------
  // read velocity to pressure map
  //---------------------------------------------------------------------------
  array_int = NULL;
  
  ReadArray_int( "velocity_pressure_map", &number, &array_int );
  
  assert_error( number == mesh->NbOfNodes,
    "velocity_pressure_map doesn't match the number of nodes of the mesh");
  
  AllocVint( mesh->NbOfNodes,mesh->VelToPreNodeMap );
  
  for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
  {
    mesh->VelToPreNodeMap[NodeId] = array_int[NodeId-1];
    if ( mesh->VelToPreNodeMap[NodeId] != 0 ) mesh->NbOfPressureNodes++;
  }  
  
  free(array_int);
  
  //---------------------------------------------------------------------------
  // read opening geometric data
  //---------------------------------------------------------------------------
  ReadOpenings( "openings", mesh );
  
  //---------------------------------------------------------------------------
  // read number of free nodes
  //---------------------------------------------------------------------------
  number = 0;
  int* NbOfNodesFreeNodes = NULL;
  
  // for velocity
  ReadArray_int("free_nodes_velocity", &number, &NbOfNodesFreeNodes);
  mesh->NbOfFreeVelocityNodes = NbOfNodesFreeNodes[0];

  free(NbOfNodesFreeNodes);
  NbOfNodesFreeNodes = NULL;
  
  // for pressure
  ReadArray_int("free_nodes_pressure", &number, &NbOfNodesFreeNodes);
  mesh->NbOfFreePressureNodes = NbOfNodesFreeNodes[0];
  
  free(NbOfNodesFreeNodes);  

  //---------------------------------------------------------------------------
  // read planes definition, if any
  //---------------------------------------------------------------------------
  GetCollisionPlanes("collision_planes");
  
  // finished reading
  // go back to previous workind directory  
  assert_GoBackToDirectory(WorkingDirectory);
  free(WorkingDirectory);

  //---------------------------------------------------------------------------
  // Print infos about the mesh
  //---------------------------------------------------------------------------
  info("\n");
  info(INT_FMT" velocity nodes, "INT_FMT" free\n",
       mesh->NbOfNodes, mesh->NbOfFreeVelocityNodes );
  info(INT_FMT" pressure nodes, "INT_FMT" free\n",
       mesh->NbOfPressureNodes, mesh->NbOfFreePressureNodes);
  info("%u elements\n", mesh->NbOfElements);
  info("\n");  
  
  //---------------------------------------------------------------------------
  // get inverse connectivity table
  //---------------------------------------------------------------------------
  mesh->InvConTab = GetInverseConnectivity(mesh->NbOfNodes, mesh->NbOfElements,
                                           NODES_PER_EL, mesh->ConTab);
  
  //---------------------------------------------------------------------------
  // get quadrature
  //---------------------------------------------------------------------------
  mesh->gauss = GaussCreate( NODES_PER_EL, 5 );
  
  //---------------------------------------------------------------------------
  // get NNS structure
  //---------------------------------------------------------------------------
  mesh->nns = PrepareNNS(mesh->NbOfNodes, mesh->points, mesh->ConTab,
                         mesh->InvConTab );

  //---------------------------------------------------------------------------
  // get minimum distance in mesh around each node
  //---------------------------------------------------------------------------
  mesh->courant = GetCourant(mesh->NbOfNodes, mesh->ConTab, mesh->InvConTab,
                             mesh->points );

  //---------------------------------------------------------------------------
  // This is usefull for particle meshing
  //---------------------------------------------------------------------------
  AllocMint(4,4,mesh->ConTabLocal);
  
  mesh->ConTabLocal[ 1 ][ 2 ] = mesh->ConTabLocal[ 2 ][ 1 ] = 5;
  mesh->ConTabLocal[ 1 ][ 3 ] = mesh->ConTabLocal[ 3 ][ 1 ] = 7;
  mesh->ConTabLocal[ 1 ][ 4 ] = mesh->ConTabLocal[ 4 ][ 1 ] = 8;
  mesh->ConTabLocal[ 2 ][ 3 ] = mesh->ConTabLocal[ 3 ][ 2 ] = 6;
  mesh->ConTabLocal[ 2 ][ 4 ] = mesh->ConTabLocal[ 4 ][ 2 ] = 9;
  mesh->ConTabLocal[ 3 ][ 4 ] = mesh->ConTabLocal[ 4 ][ 3 ] = 10;
  
  //---------------------------------------------------------------------------
  // version z stuff
  //---------------------------------------------------------------------------
#ifdef VERSION_Z
  // allocate outside nodes list
  mesh->OutsideNodes = (bool*) calloc( mesh->NbOfNodes+1 , sizeof(bool) );
#endif
  
  return mesh;
}
