/** \file
 Mass conservation of inner domain for 2 domains version of spaf.
 */

#include "mass_conservation.h"
#include "memory.h"
#include "linalg.h"
#include "strings.h"
#include "statistics.h"
#include "mesh_search.h"
#include "data_io.h"
#include "shape_functions.h"

//#define DBG

/// Mass conserving operator, for version z
static double
**_normal, ///< vector normal to surface at each nodes
*_weight,  ///< weight to correct normal velocity
**_weighted_normal; ///< vector used to compute integral of normal velocity
static int
_NbOfEntries; //< number of entries of the operator
// this is the value above which the mass conserving filtering is done
static double _net_flux_treshold = 1.e-10;
// length of the vector for checking whether the normal to an element is pointing outside
// note : it has to be at least bigger than the float epsilon machine, otherwise when using ANN library with float types is wrong when checking the tip of the vector
static double _vector_length = 1.e-5;

//=============================================================================
/** Make the mass conserving operator, for version z.
 
 This operator is used to filter the velocity field on the surface of the inner mesh in such a way it is mass conserving.
 We assume that surface elements are flat, this is the case for version z since the inner domain should be a parallelepipedic box, which is required by the convection routine to run much faster. We also assume that the nodes on the boundary of the inner mesh have the highest indexes, as these are on the Dirichlet boundary.
 The operators allow to compute the net flux I over the boundary S, with normal n :
 
 I = \Int_S u \ddot n = \Sum_{i \in nodes on S} u_i \ddot \Sum_{k \in elements that contain i} ( \Int_{element k} shapefunction_i ) n_{element k}
 
 To compute the inner integral, we use a quadrature based on all the local nodes of an element, with equal weights of 1/6 of the element surface area. This correspond to a lumped matrix, the precision of this quadrature is between order 1 and 2, as it blends a quadrature of order 1 (3 vertices) and of order 2 (3 edge mid points). This is enough for this purpose and the operators are diagonals.
 So we have the weighted_normal vector such that : I = \Sum_{i \in nodes on S} u_i \ddot weighted_normal_i
 
 We also use the normal at a node as the normed weighted_normal, and the weight at a node is its norm.
 
 Then the filtering goes as follows, we assume the filtered velocity at node i u_i is such that
 u_i = v_i + lambda c_i n_i
 with v_i the velocity to be filtered, lambda a global correction coefficient and c_i a local correction coefficient, and n_i the normal as defined above.
 Then I_u = 0   =>   lambda = - I_v / \Sum_{i \in nodes on S} c_i n_i \ddot weighted_normal_i
 = - I_v / \Sum_{i \in nodes on S} c_i weight_i
 
 then finally u_i is corrected.
 
 */
//=============================================================================
void
GetMassConservingOperator(
                          const char*   PathToCase, ///< path to case for inner mesh
                          const mesh_t* mesh )      ///< inner mesh
{
  info("\n""Preparing mass conserving operator\n");
  
  // build path to surface connectivity
  char *PathToConnectivityZ = StringConcat(PathToCase, "/mesh/connectivity.z");
  
  // read raw connectivity array
  int number = 0;
  int *array_int = NULL;
  
  ReadArray_int( PathToConnectivityZ, &number, &array_int );
  
  free(PathToConnectivityZ);
  
  // create connectivity array
  int NbOfNodesPerElements = 6,
  NbOfElements = number / NbOfNodesPerElements;
  int **connectivity = NULL;
  
  info("\t"INT_FMT" elements on surface\n",NbOfElements);
  
  AllocMint( NbOfElements, NbOfNodesPerElements, connectivity );
  
  // we need the number of nodes in that surface mesh, we assume that those nodes have been numberer with the highiest indexes in the domain mesh, since this is a boundary, such that we have the following :
  int FirstNodeId = mesh->NbOfFreeVelocityNodes + 1,
  NbOfNodes = mesh->NbOfNodes - FirstNodeId + 1;
  
  info("\t"INT_FMT" nodes on inner mesh surface""\n",NbOfNodes);
  
  // fill it from raw array by offsetting node indexes such that first index is 1, the reason is that GetInverseConnectivity assumes indexes start at 1, moreover it is more convenient for what comes next
  int id = 0;
  for ( int element = 1 ; element <= NbOfElements ; element++ )
    for ( int node = 1 ; node <= NbOfNodesPerElements ; node++ )
      connectivity[element][node] = array_int[id++] - FirstNodeId + 1;
  
  free(array_int);
  
  // get inverse connectivity table
  SparseMatrix_t *InvCon = GetInverseConnectivity(NbOfNodes, NbOfElements,
                             NbOfNodesPerElements, connectivity );
  
  // arrays that contain the elements normals and surface areas
  double **normal = NULL, *area = NULL;
  AllocMdouble(NbOfElements, 3, normal);
  AllocVdouble(NbOfElements, area);
  
  // fill those arrays
  for ( int elem = 1 ; elem <= NbOfElements ; elem++ )
  {
    // indexes of nodes in this element, de-offset index because we access nodes data from the inner mesh
    int nodes[6+1] = {0,0,0,0,0,0,0};
    
    for ( int i = 1 ; i <= NbOfNodesPerElements ; i++)
      nodes[i] = connectivity[elem][i] + FirstNodeId - 1;
    
    // get the normal vector, we assume faces are flat
    // get the coordinates of the 3 vertices of the triangle element
    double *point[3] = {
      mesh->points[nodes[1]],
      mesh->points[nodes[2]],
    mesh->points[nodes[3]]},
    // get the 2 edge vectors from first point
    u[3] = {
      point[1][1]-point[0][1],
      point[1][2]-point[0][2],
    point[1][3]-point[0][3] },
    v[3] = {
      point[2][1]-point[0][1],
      point[2][2]-point[0][2],
    point[2][3]-point[0][3] };
    
    // compute the normal vector
    normal[elem][1] = u[1]*v[2]-u[2]*v[1];
    normal[elem][2] = u[2]*v[0]-u[0]*v[2];
    normal[elem][3] = u[0]*v[1]-u[1]*v[0];
    
    // get the norm
    double norm = nrm23(normal[elem]);
    
    // the norm is twice the surface area of the triangle
    area[elem] = norm / 2.;
    
    // make normal normed
    scal3(1./norm,normal[elem]);
    
    // check whether the normal is pointing outward
    // for this we check the vertex at the tip of the normal is in the mesh, but to avoid that the normal goes through out the mesh, we scale it so that the tip is very close to the surface
    
    // first get the base vertex of the normal, it is the average of all vertex points
    double base[3] = {
      ( point[0][1] + point[1][1] + point[2][1] ) / 3.,
      ( point[0][2] + point[1][2] + point[2][2] ) / 3.,
      ( point[0][3] + point[1][3] + point[2][3] ) / 3. },
    
    // get the tip point
    tip[4] = {0.,
      base[0]+_vector_length*normal[elem][1],
      base[1]+_vector_length*normal[elem][2],
      base[2]+_vector_length*normal[elem][3] };
    
    // now check if the tip is in the mesh
    // start the search in the mesh from the first node of the element
    int closest_node = nodes[1];
    
    // if the tip is inside then reverse normal orientation
    if ( SearchElemContainingPoint(mesh, true, tip, &closest_node, NULL) != 0 )
      scal3(-1.,normal[elem]);
    
    // FOR DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DBG
    debug("elem %4d normal ="DBL_FMT" "DBL_FMT" "DBL_FMT,elem,normal[elem][1],normal[elem][2],normal[elem][3]);
    debug(", area ="DBL_FMT"\n",area[elem]);
#endif
    // FOR DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  }
  
  // now build the operators
  
  // allocate them
  _NbOfEntries = NbOfNodes;
  AllocMdouble( NbOfNodes, 3, _normal);
  AllocMdouble( NbOfNodes, 3, _weighted_normal);
  AllocVdouble( NbOfNodes, _weight );
  
  // for info and checking, the total surface area
  double surface_area = 0.;
  
  // fill them
  for ( int node = 1 ; node <= NbOfNodes ; node++ )
  {
    // get the weighted normal, i.e. the contribution of the node
    double normal_av[4] = {0.,0.,0.,0.},
    weight_av = 0.;
    
    // add up contribution from all elements that contain the current node
    for ( int ElemId = 1 ; ElemId <= NbOfEntriesInRow(InvCon, node-1) ; ElemId++ )
    {
      // get the current element index
      int elem = EntryNode(InvCon, node-1, ElemId-1);
      
      // weight of the lumped mass matrix
      double weight = area[elem] / 6.;
      
      // update weight sum
      weight_av += weight;
      
      // update surface area sum
      surface_area += weight;
      
      // update weighted normal sum
      axpby3( weight, normal[elem], 1., normal_av);
    }
    
    // get weighted normal, and normal as the normed weighted normal, and the weight as this latter norm
    copy3( normal_av, _weighted_normal[node] );
    copy3( _weighted_normal[node], _normal[node] );
    double weighted_normal_norm = nrm23(_normal[node]);
    scal3( 1./weighted_normal_norm, _normal[node] );
    _weight[node] = weighted_normal_norm;
    
    // FOR DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DBG
    int NodeId = node + FirstNodeId - 1; // node index in inner mesh
    debug("node %4d normal          ="DBL_FMT" "DBL_FMT" "DBL_FMT"\n",NodeId,
          _normal[node][1],
          _normal[node][2],
          _normal[node][3]);
    debug("node %4d normal weighted ="DBL_FMT" "DBL_FMT" "DBL_FMT"\n",NodeId,
          _weighted_normal[node][1],
          _weighted_normal[node][2],
          _weighted_normal[node][3]);
    debug("node %4d weigth ="DBL_FMT"\n",NodeId,_weight[node]);
#endif
    // FOR DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  }
  
  info("\n\tsurface area = sum of weights = "DBL_FMT"\n",surface_area);
  
  FreeM(normal);
  free(area);
  free(connectivity);
  FreeSparseMatrix(InvCon);
}
//=============================================================================
/** Compute Z bcs by interpolation from G, this is to go after SetVelBC !
 
 We don't actually fill the velocity structure of BC, but instead prescribe
 values of Vel at the version_z nodes on the corresponding bc. Thus this routine
 MUST be called after SetVelBC, because this one will prescribe Vel for all
 non free nodes, and otherwise will overwrite values.
 */
//=============================================================================
void
MassConserveBC(
const double    position[4],  ///< particle position
const mesh_t*   mesh_o, ///< pointer to outer mesh structure
      double**  Vel_o,  ///< pointer to outer fluid variables structure
const mesh_t*   mesh,   ///< pointer to inner mesh structure
      double**  Vel)    ///< pointer to inner fluid variables structure
{
  // FOR DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if 0
  // set Vel_o = (x,0,0) - > div(Vel_o) = 1
  for ( int node = 1 ; node <= mesh_o->NbOfNodes ; node++ )
  {
    Vel_o[1][node] = pow( mesh_o->points[node][1] - 0. ,1);
    //    Vel_o[1][node] = pow( mesh_o->points[node][1] - 1. ,2);
    Vel_o[2][node] = 0.;
    //    Vel_o[2][node] = pow( mesh_o->points[node][2] - 1. ,2);
    Vel_o[3][node] = 0.;
    //    Vel_o[3][node] = pow( mesh_o->points[node][3] - 1. ,2);
  }
#endif
  // FOR DEBUGGING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  
  // first prescribed node
  int FirstPrescribedNode = mesh->NbOfFreeVelocityNodes + 1;
  
  // lambda = - net_flux / lambda2, correction coef that is node dependent so that = 
  double net_flux = 0.,
  lambda2 = 0.,
  *correction = NULL;
  
  AllocVdouble(_NbOfEntries,correction);
  
  for ( int NodeId = FirstPrescribedNode ; NodeId <= mesh->NbOfNodes ; NodeId++ )
  {    
    // if node is outside, impose 0 bc and skip to next node
    if ( mesh->OutsideNodes[NodeId] == true )
    {
      for ( int dir = 1 ; dir <= 3 ; dir++ )
        Vel[dir][NodeId] = 0.;
      
      continue;
    }
    
    double
    *point = mesh->points[NodeId],
    NodeCoord[4] = {0., position[1]+point[1], position[2]+point[2], position[3]+point[3]}, // node coordinates in mesh_o frame
    node_LocalCoord_o[4] = {0.,0.,0.,0.};  // local G node coordinates
    
    // find outer element containing inner mesh node
    int temp = 0,
    ElemContainer_o = SearchElemContainingPoint(
                                                mesh_o, true, NodeCoord, &temp, &temp );
    
    // check the node is not outside the global mesh
    assert_error( ElemContainer_o != 0, "Node %u is not in global mesh", NodeId );
    
    // get the local coordinates of the Z node in the G element, do checking too
    bool ElemContainerFound = MeshGetLocalPosition(mesh_o, ElemContainer_o,
                              NodeCoord, node_LocalCoord_o);
    
    assert_error(ElemContainerFound, "Node %u is not in the element %u",
                 NodeId, ElemContainer_o );
    
    // interpolate the velocity from the outer element at the Z node
    double v[4] = {0.,0.,0.,0.};
    
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      v[dir] = InterpolateOrder2(mesh_o->ConTab[ ElemContainer_o ],
                                 node_LocalCoord_o, Vel_o[dir]);
    
    // prepare mass conserving filtering
    
    // index of the node on the surface, just offset
    int SurfaceNodeId = NodeId - FirstPrescribedNode + 1;
    
    // get node contribution = u_k \ddot n^{\omega}_k
    double contrib = dot3(v, _weighted_normal[SurfaceNodeId]);
    
    // update numerator, denominator and set node dependent correction coef
    net_flux += contrib;
    correction[SurfaceNodeId] = contrib;
    //correction[SurfaceNodeId] = op->_weight[SurfaceNodeId];
    lambda2 += correction[SurfaceNodeId] * _weight[SurfaceNodeId];
    
    //    printf("node %d contrib = "DBL_FMT"\n",NodeId,contrib);
    
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      Vel[dir][NodeId] = v[dir];
  }
  
  // statistics
  StatAdd("Mass Conservation : net flux before", fabs(net_flux));
  
  // apply filtering ?
  // if net_flux is very small, then we do not need to filter, but also lambda2 could be very small too, hence lambda may be not well controlled
  if ( fabs(net_flux) <= _net_flux_treshold ) return;
  
  // net mass flux before and after filtering
  double net_flux_filtered = 0.,
  
  // get lambda
  lambda = - net_flux / lambda2;
  
  // filter velocity
  for ( int NodeId = FirstPrescribedNode ; NodeId <= mesh->NbOfNodes ; NodeId++ )
  {
    int SurfaceNodeId = NodeId - FirstPrescribedNode + 1;
    
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      Vel[dir][NodeId] += lambda * correction[SurfaceNodeId] * _normal[SurfaceNodeId][dir];
      net_flux_filtered += Vel[dir][NodeId] * _weighted_normal[SurfaceNodeId][dir];
    }
  }
  
  // some info
  debug("\tMass Conservation : net flux : before = "DBL_FMT", after = "DBL_FMT"\n",
        net_flux, net_flux_filtered );
  
  // statistics
  StatAdd("Mass Conservation : net flux after", fabs(net_flux_filtered));
  
  free(correction);
}
