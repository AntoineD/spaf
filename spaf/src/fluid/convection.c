/** \file
 Computation of convection terms with the method of characteristics.

*/ 
#include "includes.h"
#include "logging.h"
#include "statistics.h"
#include "linalg.h"
#include "parse.h"
#include "mesh_search.h"
#include "shape_functions.h"
#include "convection.h"
#include "memory.h"
#include "version_z.h"

/// local parameters
// limit courant number over which bounding box searching is used
static double
_bb_search_courant_limit = 0.,
_bb_search_courant_limit_o = 0.,
// minimum substep allowed for the Lagrange convective sub-stepping
_min_substep_coef = 0.,
// The tolerance for a point to lie on the outlet plane
_outlet_plane_tolerance = 0.;

//=============================================================================
/** Read convection parameters.
 */
//=============================================================================
void
ReadConvectionParameters(
FILE* FileId )  ///< file identifier
{
  // read outlet_plane_tolerance
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("outlet_plane_tolerance"),
    "Missing convection parameter : ""outlet_plane_tolerance");
  _outlet_plane_tolerance = Token2double(3,'[',0.,DBL_MAX,']');

  // read min_substep_coef
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("min_substep_coef"),
    "Missing convection parameter : ""min_substep_coef");
  _min_substep_coef = Token2double(3,'[',0.,DBL_MAX,']');

  // read bb_search_courant_limit
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("bb_search_courant_limit"),
    "Missing convection parameter : ""bb_search_courant_limit");
  _bb_search_courant_limit = Token2double(3,'[',0.,DBL_MAX,']');

  // read bb_search_courant_limit for spafz
#ifdef VERSION_Z
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs("bb_search_courant_limit_o"),
    "Missing convection parameter : ""bb_search_courant_limit_o");
  _bb_search_courant_limit_o = Token2double(3,'[',0.,DBL_MAX,']');
#endif

  // check end of section
  TokenizeFileLine(FileId);
  assert_error( SectionIsEnd(),
    "Convection parameters : missing section = end");

  // print parameters
  info("\nConvection parameters :\n" ); 
  info("\tmin_substep_coef        = " "%g" "\n", _min_substep_coef);
  info("\toutlet_plane_tolerance  = " "%g" "\n", _outlet_plane_tolerance);
  info("\tbb_search_courant_limit = " "%g" "\n", _bb_search_courant_limit);
#ifdef VERSION_Z
  info("\tbb_search_courant_limit_o = " "%g" "\n", _bb_search_courant_limit_o);
#endif
}
//==============================================================================
/// Check if a characteristic has crossed an inlet or outlet.
//==============================================================================
static bool
CheckIOCrossing(
const mesh_t* mesh,     ///< mesh structure
      double  foot[4] ) ///< position of the foot of the characteristic
{
  for ( int i = 0 ; i < mesh->NbOfOpenings ; i++ )
  {
    const io_t* io = &mesh->io[i]; // alias to current io
    
    double FootToIODistance = PlaneToPointDistance( foot, &io->plane );
    
    // if the distance from foot to io plane is negative then characteristic
    // did not cross the boundary, jump to next io
    if ( FootToIODistance < 0. ) continue;
    
    // otherwise we project foot on io plane
    for ( int dir = 0 ; dir < 3 ; dir++ ) 
      foot[dir+1] -= FootToIODistance * io->normal[dir];
    
    // and check whether the projection is within the io circle radius
    double dist = sqrt( pow( foot[1] - io->center[ 0 ], 2 ) +
                        pow( foot[2] - io->center[ 1 ], 2 ) +
                        pow( foot[3] - io->center[ 2 ], 2 ) );
    
    assert_error( dist <= io->radius,
                 "Foot projection outside i/o #%d.\ndist = %g, radius = %g\n",
                 i, dist, io->radius );
    
    return true; // the characteristic crossed the boundary
  }

  return false;
}
//=============================================================================
/** Go backward in time.

 Compute the position of the foot of a characteristic and the local courant number which is a measure of the length of the characteristic in term of the number of elements spanned.
 On return, the foot and the local courant number have been updated.
 */
//=============================================================================
static void
GoBackward(
const int     level,        ///< time level : 1 or 2
const double  courant,      ///< min internodal distance around
const double  head[4],      ///< position of the characteristic head
const double  VelHead[4],   ///< velocity of the characteristic head
const double  FrameVel[4],  ///< moving frame velocity
const double  dt,           ///< time step
const double  sub_timestep, ///< sub-stepping time step
      double* foot,         ///< position of the characteristic foot 
      double* CourantLocal )///< number of element spanned
{
  assert_error(sub_timestep >= _min_substep_coef * dt,
               "Convection sub-stepping requires too many substeps, \
                try to decrease the time step and try again!");
  
  double factor = sub_timestep * level,
  displacement[4] = {0.,0.,0.,0.};
  
  for (int dir = 1 ; dir <= 3 ; dir++ )
  {
    displacement[ dir ] = factor * ( VelHead[ dir ] - FrameVel[dir] );
    foot[dir] = head[ dir ] - displacement[ dir ];
  }
  
  *CourantLocal = nrm23( displacement ) / courant;
}
//==============================================================================
/** Compute extrapolated velocities.
 */
//==============================================================================
static void
ExtrapolateVelocity(
const int     order,
const int     NbOfNodes,
const int     NbOfFreeVelocityNodes,
      double* VelOld2[4],
      double* VelOld1[4],
      double* VelTilde1[4],
      double* VelTilde2[4],
      double* VelExtrapol[4])
{
  double coeff1 = 1., // order 1 by default,
  coeff2 = 0.; // coefficients for the extrapolation scheme
  
  // In case that we start with O velocity field the extrapolation to the
  // third time level can be awfully inacurate thus it is better to use
  // 1-st order extrapolation. So we compute the max norm of the n-2 velocity
  // and decide with respect to its value.
  if ( order == 2 )
  {
    double vel_max = 0.; // maximum of the velocity at level n-2
    
    for ( int node = 1; node <= NbOfNodes; node++ )
    {
      double VelTemp[4];
      
      for ( int dir = 1 ; dir <= 3 ; dir++ )
        VelTemp[dir] = VelOld2[dir][node];
      
      double norm = nrm23( VelTemp );
      
      if ( norm > vel_max ) vel_max = norm ;
    }
    
    // if the velocity is high enough, use second order
    if ( vel_max >= 1e-6 )
    {
      coeff1 = 2.;
      coeff2 = -1.;
    }
  }
  
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    if ( order == 1 )
    {
      // extrapolate free nodes
      for ( int node = 1 ; node <= NbOfNodes ; node++ )
        VelExtrapol[ dir ][ node ] = VelOld1[ dir ][ node ];
      
      // prescribe non-free nodes, the extrapolated velocity may be used
      // for sub-stepping near boundaries
      for ( int node = NbOfFreeVelocityNodes + 1 ; node <= NbOfNodes ; node++ )
        VelTilde1[ dir ][ node ] = VelOld1[ dir ][ node ];
    }
    else
    {
      // extrapolate free nodes
      for ( int node = 1 ; node <= NbOfNodes ; node++ )
        VelExtrapol[ dir ][ node ] =
        coeff1 * VelOld1[ dir ][ node ] +
        coeff2 * VelOld2[ dir ][ node ];
      
      // prescribe non-free nodes, the extrapolated velocity may be used
      // for sub-stepping near boundaries
      for ( int node = NbOfFreeVelocityNodes + 1 ; node <= NbOfNodes ; node++ )
      {
        VelTilde1[ dir ][ node ] = VelOld1[ dir ][ node ];
        VelTilde2[ dir ][ node ] = VelOld2[ dir ][ node ];
      }
    }
  }  
}
//==============================================================================
/** Compute convection velocities.

 Advances convective terms in time using a subcyling marching procedure,
 as described by Maday, Patera and Ronquist, J. Sci. Comp., 5:265, 1990.
 This is for first and/or order splitting using Lagrange integration over
 characteristics.
 */
//==============================================================================
void
Convection(
const int order, ///< order of the discretisation
const mesh_t* mesh, ///< mesh structure
const double dt, ///< time step
const double FrameVel[4], ///< frame velocity
      double** VelOld1, ///< velocity field at the previous time step
      double** VelOld2, ///< velocity field at the previous previous time step
      double** VelTilde1, ///< convected velocity at the previous time step
      double** VelTilde2 ) ///< convected velocity at the previous previous time step
{
  debug( "\nVelocity convection step\n" );
  
  // holds the extrapolated velocity
  double *VelExtrapol[3+1] = {NULL, NULL, NULL, NULL};
  for (int i = 1; i <= 3; i++)
    AllocVdouble(mesh->NbOfNodes, VelExtrapol[i]);
  
  // extrapolate the velocity
  ExtrapolateVelocity(order, mesh->NbOfNodes, mesh->NbOfFreeVelocityNodes,
                      VelOld2, VelOld1, VelTilde1, VelTilde2, VelExtrapol);
  
  double CourantLocalMax = 0.; // maximum local courant number (just for info)
  int SubstepNbMax = 0, // maximum number of substeps (just for info)
  ElemCheckedNb = 0; // total # of elements checked (just for info)

  //--------------------------------------------------------
  // Loop over the free points of the mesh
  //--------------------------------------------------------
#pragma omp parallel for \
default(none) \
shared(mesh,VelExtrapol,VelOld1,VelOld2,VelTilde1,VelTilde2,NbOfNodesConvected) \
shared(_bb_search_courant_limit,_min_substep_coef,FrameVel) \
shared(CourantLocalMax) \
reduction(+:SubstepNbMax,ElemCheckedNb)
  for ( int node = 1 ; node <= (int) mesh->NbOfFreeVelocityNodes ; node++ )
  for ( int level = 1 ; level <= order ; level++ ) // Loop over the level in time
  {
    //=============================================================================
    // Euler backwards integration of characteristic equation
    //=============================================================================
    // position of the foot and head of the characteristic through the current node,
    // for sub-stepping, this is the position of the current and previous 
    // lagrangian node
    double foot[4] = {0.,0.,0.,0.},
    // head is initialized with current node position
    head[4] = { 0.,
      mesh->points[node][1],
      mesh->points[node][2],
      mesh->points[node][3] },
    // velocity at the current head of the characterstic, initialize with the one of the current node
    VelHead[4] = { 0.,
      VelExtrapol[1][node],
      VelExtrapol[2][node],
      VelExtrapol[3][node] },
    CourantLocal = 0., // local Courant number
    // local courant number from the bbox (min internodal distance around each point), initialized with current node
    courant = mesh->courant[ node ],
    sub_time = dt, // sub-time
    sub_timestep = dt, // sub-timestep
    **VelOldPtr = ( level == 1 ? VelOld1 : VelOld2 ), // velocity used for interpolation
    **VelTildePtr = ( level == 1 ? VelTilde1 : VelTilde2 ), // convected velocity
    **VelInterpFromPtr = NULL; // velocity from which the interpolation is done

    GoBackward(level, courant, head, VelHead, FrameVel, dt,
                sub_timestep, foot, &CourantLocal );
    
    // update max local courant number
    if ( CourantLocal > CourantLocalMax ) CourantLocalMax = CourantLocal;
    
    // closest point of the mesh, by default it starts to search from the head of the characteristics
    int ClosestNode = node;
    
    //=============================================================================
    // Here we begin sub-stepping in case that the foot is not found inside the 
    // Eulerian geometry; otherwise we simply make 1 substep
    //=============================================================================
    bool substep = true;
    int SubstepNb = 0;
    
    while ( substep )
    {
      // we use the bbox search only if the local courant number is lower than a treshold, which in general should be 1.
      bool UseNNS = ( CourantLocal > _bb_search_courant_limit );

      // search the element containing the point
      int ElemContainer = SearchElemContainingPoint(
        mesh, UseNNS, foot, &ClosestNode, &ElemCheckedNb );
      
      if ( ElemContainer == 0 )
      {
        // The point is not in the eulerian geometry,
        // we first check whether it crossed an inlet or outlet boundary
        bool crossed = CheckIOCrossing( mesh, foot );
        
        // if the caracteristic crossed an inlet or outlet, it has been projected
        // onto the corresponding boundary, otherwise, we need to do sub-stepping
        // from the original position of the node, in such a way the whole
        // characteristic computed with the same sub-time-step in inside the domain.
        if ( !crossed )
        {
          warning("%s : node %d foot outside""\n", __FUNCTION__, node );
          warning("\tfoot position : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  foot[1],foot[2],foot[3] );
          warning("\thead position : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  head[1],head[2],head[3] );
          warning("\thead velocity : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  VelHead[1],VelHead[2],VelHead[3] );
          warning("\tsub time step : "DBL_FMT"\n", sub_timestep );
          
          // start all over again from the original node position
          for ( int dir = 1 ; dir <= 3 ; dir++ )
          {
            head[dir] = mesh->points[node][dir];
            VelHead[dir] = VelExtrapol[dir][node];
          }
          
          courant = mesh->courant[ node ];
          
          // half sub-timestep and initialize sub-time
          sub_timestep /= 2.;
          sub_time = sub_timestep;

          // do sub-stepping
          GoBackward(level, courant, head, VelHead, FrameVel, dt,
                      sub_timestep, foot, &CourantLocal );
          
          SubstepNb++;
        }
        continue; // continues the while loop (substep)
      }
      
      //=============================================================================
      // the point is in the eulerian geometry, we can interpolate the velocity and check substepping is done
      //=============================================================================
      // The local coordinates of the foot, where we interpolate the velocity.
      double LocalPos[4] = {0.,0.,0.,0.};
      
      // get the local coordinates of the foot
      // also check the point is in the element
      bool ElemCheck = MeshGetLocalPosition(
        mesh, ElemContainer, foot, LocalPos );
      
      // give info if the point is not inside the element
      if ( !ElemCheck )
      {
        double *coord = mesh->points[node];
        
        info("node " INT_FMT " (" DBL_FMT "," DBL_FMT "," DBL_FMT
             ") is not in element " INT_FMT "\n",
             node,coord[1],coord[2],coord[3],ElemContainer );
        info("foot = (" DBL_FMT "," DBL_FMT "," DBL_FMT ")" "\n",
             foot[1], foot[2], foot[3]);
        error("%s : lagrangian point %u not in the element %u",
              __FUNCTION__, node, ElemContainer );
      }
      
      // is the sub-stepping done ?
      if ( fabs( dt - sub_time ) / dt >= _min_substep_coef * 1e-2 )
      {
        //info("sub_time = "DBL_FMT"\n",sub_time);
        
        for ( int dir = 1 ; dir <= 3 ; dir++ )
        {
          // Interpolate the velocity at the current foot and substep and use it for the next sub-stepping
          VelTildePtr[dir][node] = InterpolateOrder2(
            mesh->ConTab[ ElemContainer ], LocalPos, VelExtrapol[dir]);
          
          // update head position with the current foot
          head[dir] = foot[dir];
          
          // update the velocity at the head of the characteristic
          VelHead[dir] = VelTildePtr[dir][node];
        }
        
        // update the local courant number from the closest node
        courant = mesh->courant[ClosestNode];
        
        // update sub time
        sub_time += sub_timestep;

        // do sub-stepping
        GoBackward(level, courant, head, VelHead, FrameVel, dt,
                    sub_timestep, foot, &CourantLocal );
        
        SubstepNb++;
        
        continue; // continue the while loop (substep)
      }
      
      //=============================================================================
      // The sub-stepping is done, interpolate the velocity and store it into VelTilde*
      //=============================================================================
      VelInterpFromPtr = VelOldPtr;

      for ( int dir = 1 ; dir <= 3 ; dir++ )
        VelTildePtr[dir][node] = InterpolateOrder2(
          mesh->ConTab[ ElemContainer ], LocalPos, VelInterpFromPtr[dir]);
      
      substep = false; // done with this node, get out of the while loop (substep)
    }
    
    if ( SubstepNb > SubstepNbMax ) SubstepNbMax = SubstepNb;
  }
  
  StatAdd("Fluid convection : elements checked", ElemCheckedNb);
  StatAdd("Fluid convection : maximum courant number", CourantLocalMax);
  debug( "\t"INT_FMT" elements checked\n", ElemCheckedNb );
  debug( "\tmaximum courant number = "DBL_FMT"\n", CourantLocalMax );
  debug( "\tmaximum number of substeps = "INT_FMT"\n", SubstepNbMax );

  for (int i = 1; i < 4; i++) free(VelExtrapol[i]);
}
//==============================================================================
/** For version z
*/
//==============================================================================
void
ConvectionMicroGrid(
const int     order, ///< order of the discretisation
const mesh_t* mesh, ///< mesh structure
const double  dt, ///< time step
const double  FrameVel[4], ///< frame velocity
const double  position[4], ///< particle position
const mesh_t* mesh_o, ///< pointer to the mesh structure of G
      double** Vel_o, ///< velocity field of G
      double** VelOld1, ///< velocity field at the previous time step
      double** VelOld2, ///< velocity field at the previous previous time step
      double** VelTilde1, ///< convected velocity at the previous time step
      double** VelTilde2 ) ///< convected velocity at the previous previous time step
{
  debug( "\nVelocity convection step\n" );
  
  // holds the extrapolated velocity
  double *VelExtrapol[3+1] = {NULL, NULL, NULL, NULL};
  for (int i = 1; i <= 3; i++)
    AllocVdouble(mesh->NbOfNodes, VelExtrapol[i]);
  
  // extrapolate the velocity
  ExtrapolateVelocity(order, mesh->NbOfNodes, mesh->NbOfFreeVelocityNodes,
                      VelOld2, VelOld1, VelTilde1, VelTilde2, VelExtrapol);
  
  double CourantLocalMax = 0.; // maximum local courant number (just for info)
  int SubstepNbMax = 0, // maximum number of substeps (just for info)
  ElemCheckedNb = 0, // total # of elements checked (just for info)
  NbOfFeetOutsideInnerDomain = 0; // number of feet outside inner domain (just for info)

  //--------------------------------------------------------
  // Loop over the free points of the mesh
  //--------------------------------------------------------
#pragma omp parallel for \
default(none) \
shared(mesh,VelExtrapol,VelOld1,VelOld2,VelTilde1,VelTilde2,NbOfNodesConvected) \
shared(_bb_search_courant_limit,_min_substep_coef,FrameVel) \
shared(_bb_search_courant_limit_o,Vel_o,mesh_o) \
shared(CourantLocalMax,position) \
reduction(+:SubstepNbMax,ElemCheckedNb,NbOfFeetOutsideInnerDomain)
  for ( int node = 1 ; node <= (int) mesh->NbOfNodes ; node++ ) 
  for ( int level = 1 ; level <= order ; level++ ) // Loop over the level in time
  {
    //=============================================================================
    // Euler backwards integration of characteristic equation
    //=============================================================================
    // position of the foot and head of the characteristic through the current node,
    // for sub-stepping, this is the position of the current and previous 
    // lagrangian node
    double foot[4] = {0.,0.,0.,0.},
    // head is initialized with current node position
    head[4] = { 0.,
      mesh->points[node][1],
      mesh->points[node][2],
      mesh->points[node][3] },
    // velocity at the current head of the characterstic, initialize with the one of the current node
    VelHead[4] = { 0.,
      VelExtrapol[1][node],
      VelExtrapol[2][node],
      VelExtrapol[3][node] },
    CourantLocal = 0., // local Courant number
    // local courant number from the bbox (min internodal distance around each point), initialized with current node
    courant = mesh->courant[ node ],
    sub_time = dt, // sub-time
    sub_timestep = dt, // sub-timestep
    **VelOldPtr = ( level == 1 ? VelOld1 : VelOld2 ), // velocity used for interpolation
    **VelTildePtr = ( level == 1 ? VelTilde1 : VelTilde2 ), // convected velocity
    **VelInterpFromPtr = NULL; // velocity from which the interpolation is done

    // if a node is outside the outer domain, set convection velocity to 0 and jump to next node.
    if ( mesh->OutsideNodes[node] )
    {
      for ( int dir = 1 ; dir <= 3 ; dir++ )
        VelTildePtr[dir][node] = 0.;
      
      continue;
    }

    GoBackward(level, courant, head, VelHead, FrameVel, dt,
                sub_timestep, foot, &CourantLocal );
    
    // update max local courant number
    if ( CourantLocal > CourantLocalMax ) CourantLocalMax = CourantLocal;
    
    // the following comment is for version z
    // After this first step, the point can be in one of the following :
    //    - Z
    //    - G but not Z
    //    - outside G and cross G's IO, then projected
    //    - outside G and doesn't cross G's IO, then we do substepping
    //! Note : Z's IO (outlet) are ignored
    
    // tells what is the current mesh
    const mesh_t *meshPtr = mesh;
    
    // closest point of the mesh, by default it starts to search from the head of the characteristics
    int ClosestNode = node;

    //=============================================================================
    // Here we begin sub-stepping in case that the foot is not found inside the 
    // Eulerian geometry; otherwise we simply make 1 substep
    //=============================================================================
    bool substep = true;
    int SubstepNb = 0;

    while ( substep )
    {
      // we use the bbox search only if the local courant number is lower than a treshold, which in general should be 1.
      bool UseNNS = ( CourantLocal > _bb_search_courant_limit );

      // for outer mesh
      if ( meshPtr == mesh_o )
//      {
        UseNNS = ( CourantLocal > _bb_search_courant_limit_o );
//        for ( int dir = 1 ; dir <= 3 ; dir++ )
//          foot[dir] += position[dir];
//      }

      // search the element containing the point
      int ElemContainer = SearchElemContainingPoint(
        meshPtr, UseNNS, foot, &ClosestNode, &ElemCheckedNb );

      if ( !ElemContainer && ( meshPtr == mesh ) )
      {
        // The foot is outside the inner domain, we search in the outer domain
        // change frame of foot from inner to outer
        for ( int dir = 1 ; dir <= 3 ; dir++ )
          foot[dir] += position[dir];

        // use bb search for closest node, because the closest node guess provided by inner mesh (node) is irrelevant to outer one.
        ElemContainer = SearchElemContainingPoint(
          mesh_o, true, foot, &ClosestNode, &ElemCheckedNb );
        
        // we'll use the outer mesh from now on
        meshPtr = mesh_o;
      }

      if ( !ElemContainer )
      {
        // The point is not in the macro grid,
        // has it crossed an inlet or outlet boundary ?
        bool crossed = CheckIOCrossing(mesh_o, foot);
        
        // if the caracteristic crossed an inlet or outlet, it has been projected
        // onto the corresponding boundary, otherwise, we need to do sub-stepping
        // from the original position of the node, in such a way the whole
        // characteristic computed with the same sub-time-step in inside the domain.
        if ( !crossed )
        {
          // THIS IS A HACK TO MAKE IT WORK !!!!!!!!!!!!!!!!!!!!!!!!!!!
          // the foot is outside the outer domain and has not crossed an i/o
          // then it should be in a wall, so we set the convected velocity to 0
          // and jump to next node
          for ( int dir = 1 ; dir <= 3 ; dir++ )
            VelTildePtr[dir][node] = 0.;
          
          substep = false;
          continue;
          // END OF HACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          warning("%s : node %d foot outside""\n", __FUNCTION__, node );
          warning("\tfoot position : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  foot[1],foot[2],foot[3] );
          warning("\thead position : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  head[1],head[2],head[3] );
          warning("\thead velocity : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  VelHead[1],VelHead[2],VelHead[3] );
          warning("\tsub time step : "DBL_FMT"\n", sub_timestep );
          
          // start all over again from the original node position in the micro grid
          for ( int dir = 1 ; dir <= 3 ; dir++ )
          {
            head[dir] = mesh->points[node][dir];
            VelHead[dir] = VelExtrapol[dir][node];
          }
          
          courant = mesh->courant[ node ];
          
          // half sub-timestep and initialize sub-time
          sub_timestep /= 2.;
          sub_time = sub_timestep;

          // framevel here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          // do sub-stepping
          GoBackward(level, courant, head, VelHead, FrameVel, dt,
                      sub_timestep, foot, &CourantLocal );
          
          SubstepNb++;

          // we'll use the inner mesh from now on
          meshPtr = mesh;
        }
        continue; // continues the while loop (substep)
      }
      
      //=============================================================================
      // the point is in the eulerian geometry, we can interpolate the velocity and check substepping is done
      //=============================================================================
      // The local coordinates of the foot, where we interpolate the velocity.
      double LocalPos[4] = {0.,0.,0.,0.};
      
      // get the local coordinates of the foot
      // also check that the point is in the element
      bool ElemCheck = MeshGetLocalPosition(
        meshPtr, ElemContainer, foot, LocalPos);
      
      // give info if the point is not inside the element
      if ( !ElemCheck )
      {
        info("are we in inner mesh? : " "\n", meshPtr==mesh);
        
        double *coord = mesh->points[node];
        
        info("node " INT_FMT " (" DBL_FMT "," DBL_FMT "," DBL_FMT
             ") is not in element " INT_FMT "\n",
             node, coord[1], coord[2], coord[3], ElemContainer );
        info("foot = (" DBL_FMT "," DBL_FMT "," DBL_FMT ")" "\n",
             foot[1], foot[2], foot[3]);
        error("%s : lagrangian point %u not in the element %u",
              __FUNCTION__, node, ElemContainer );
      }
      
      // is the sub-stepping done ?
      if ( fabs( dt - sub_time ) / dt >= _min_substep_coef * 1e-2 )
      {
        //info("sub_time = "DBL_FMT"\n",sub_time);
        
        // introduced for merging spaf and spafz, we interpolate from the extrapolated velocity
        VelInterpFromPtr = VelExtrapol;
        const double *FrameVelPtr = FrameVel;

        // if the foot is outside the inner mesh and inside the outer mesh, then we interpolate from the outer velocity field with a zero frame velocity
        double ZeroFrameVel[4] = {0.,0.,0.,0.};
        
        if ( meshPtr == mesh_o )
        {
          VelInterpFromPtr = Vel_o;
          FrameVelPtr = ZeroFrameVel;
        }

        for ( int dir = 1 ; dir <= 3 ; dir++ )
        {
          // Interpolate the velocity at the current foot and substep and use it for the next sub-stepping
          VelTildePtr[dir][node] = InterpolateOrder2(
            meshPtr->ConTab[ ElemContainer ], LocalPos, VelInterpFromPtr[dir]);
          
          // update head position with the current foot
          head[dir] = foot[dir];
          
          // update the velocity at the head of the characteristic
          VelHead[dir] = VelTildePtr[dir][node];
        }
        
        // update the local courant number from the closest node
        courant = meshPtr->courant[ClosestNode];
        
        // update sub time
        sub_time += sub_timestep;
        
        // framevel here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        // do sub-stepping
        GoBackward(level, courant, head, VelHead, FrameVelPtr, dt,
                    sub_timestep, foot, &CourantLocal );
        
        SubstepNb++;
        continue; // continue the while loop (substep)
      }
      
      //=============================================================================
      // The sub-stepping is done, interpolate the velocity and store it into VelTilde*
      //=============================================================================
      VelInterpFromPtr = VelOldPtr;
      
      // when the foot is outside the inner mesh and inside the outer mesh, we interpolate from the outer velocity field
      if ( meshPtr == mesh_o )
      {
        VelInterpFromPtr = Vel_o;
        NbOfFeetOutsideInnerDomain++;
      }
      
      for ( int dir = 1 ; dir <= 3 ; dir++ )
        VelTildePtr[dir][node] = InterpolateOrder2(
          meshPtr->ConTab[ ElemContainer ], LocalPos, VelInterpFromPtr[dir]);
      
      substep = false; // done with this node, get out of the while loop (substep)
    }
    
    if ( SubstepNb > SubstepNbMax ) SubstepNbMax = SubstepNb;
  }
  
  StatAdd("Fluid convection : elements checked", ElemCheckedNb);
  StatAdd("Fluid convection : maximum courant number", CourantLocalMax);
  StatAdd("Fluid convection : feet outside inner domain", NbOfFeetOutsideInnerDomain);
  debug( "\t"INT_FMT" elements checked\n", ElemCheckedNb );
  debug( "\tmaximum courant number = "DBL_FMT"\n", CourantLocalMax );
  debug( "\tmaximum number of substeps = "INT_FMT"\n", SubstepNbMax );
  debug( "\t" INT_FMT " feet outside inner domain""\n",NbOfFeetOutsideInnerDomain);

  for (int i = 1; i < 4; i++) free(VelExtrapol[i]);
}
//==============================================================================
/** For version z two way, outer and inner are inter-changed.
 */
//==============================================================================
void
ConvectionMacroGrid(
const int order, ///< order of the discretisation
const mesh_t* mesh, ///< mesh structure
const double dt, ///< time step
const double FrameVel[4], ///< frame velocity
const double position[4], ///< particle position
const mesh_t* mesh_o, ///< pointer to the mesh structure of G
      double** Vel_o, ///< velocity field of G
      double** VelOld1, ///< velocity field at the previous time step
      double** VelOld2, ///< velocity field at the previous previous time step
      double** VelTilde1, ///< convected velocity at the previous time step
      double** VelTilde2 ) ///< convected velocity at the previous previous time step
{
  debug( "\nVelocity convection step\n" );
  
  // holds the extrapolated velocity
  double *VelExtrapol[3+1] = {NULL, NULL, NULL, NULL};
  for (int i = 1; i <= 3; i++)
    AllocVdouble(mesh->NbOfNodes, VelExtrapol[i]);
  
  // extrapolate the velocity
  ExtrapolateVelocity(order, mesh->NbOfNodes, mesh->NbOfFreeVelocityNodes,
                      VelOld2, VelOld1, VelTilde1, VelTilde2, VelExtrapol);
  
  double CourantLocalMax = 0.; // maximum local courant number (just for info)
  int SubstepNbMax = 0, // maximum number of substeps (just for info)
  ElemCheckedNb = 0, // total # of elements checked (just for info)
  NbOfFeetInsideInnerDomain = 0; // number of feet outside inner domain (just for info)
  
  //--------------------------------------------------------
  // Loop over the free points of the mesh
  //--------------------------------------------------------
#pragma omp parallel for \
default(none) \
shared(mesh,VelExtrapol,VelOld1,VelOld2,VelTilde1,VelTilde2,NbOfNodesConvected) \
shared(_bb_search_courant_limit,_min_substep_coef,FrameVel) \
shared(_bb_search_courant_limit_o,Vel_o,mesh_o) \
shared(CourantLocalMax,position) \
reduction(+:SubstepNbMax,ElemCheckedNb,NbOfFeetInsideInnerDomain)
  for ( int node = 1 ; node <= (int) mesh->NbOfFreeVelocityNodes ; node++ ) 
  for ( int level = 1 ; level <= order ; level++ ) // Loop over the level in time
  {
    //=============================================================================
    // Euler backwards integration of characteristic equation
    //=============================================================================
    // position of the foot and head of the characteristic through the current node,
    // for sub-stepping, this is the position of the current and previous 
    // lagrangian node
    double foot[4] = {0.,0.,0.,0.},
    // head is initialized with current node position
    head[4] = { 0.,
      mesh->points[node][1],
      mesh->points[node][2],
      mesh->points[node][3] },
    // velocity at the current head of the characterstic, initialize with the one of the current node
    VelHead[4] = { 0.,
      VelExtrapol[1][node],
      VelExtrapol[2][node],
      VelExtrapol[3][node] },
    CourantLocal = 0., // local Courant number
    // local courant number from the bbox (min internodal distance around each point), initialized with current node
    courant = mesh->courant[ node ],
    sub_time = dt, // sub-time
    sub_timestep = dt, // sub-timestep
    **VelOldPtr = ( level == 1 ? VelOld1 : VelOld2 ), // velocity used for interpolation
    **VelTildePtr = ( level == 1 ? VelTilde1 : VelTilde2 ), // convected velocity
    **VelInterpFromPtr = NULL; // velocity from which the interpolation is done
    
    // if a node is outside the outer domain, set convection velocity to 0 and jump to next node.
    if ( mesh->OutsideNodes[node] )
    {
      for ( int dir = 1 ; dir <= 3 ; dir++ )
        VelTildePtr[dir][node] = 0.;
      
      continue;
    }
    
    GoBackward(level, courant, head, VelHead, FrameVel, dt,
                sub_timestep, foot, &CourantLocal );
    
    // update max local courant number
    if ( CourantLocal > CourantLocalMax ) CourantLocalMax = CourantLocal;
    
    // the following comment is for version z
    // After this first step, the point can be in one of the following :
    //    - Z
    //    - G but not Z
    //    - outside G and cross G's IO, then projected
    //    - outside G and doesn't cross G's IO, then we do substepping
    //! Note : Z's IO (outlet) are ignored
    
    // tells what is the micro mesh
    const mesh_t *meshPtr = mesh_o;
    
    // closest point of the mesh, by default it starts to search from the head of the characteristics
    int ClosestNode = node;
    
    //=============================================================================
    // Here we begin sub-stepping in case that the foot is not found inside the 
    // Eulerian geometry; otherwise we simply make 1 substep
    //=============================================================================
    bool substep = true;
    int SubstepNb = 0;
    
    while ( substep )
    {
      // we use the bbox search only if the local courant number is lower than a treshold, which in general should be 1.
      bool UseNNS = ( CourantLocal > _bb_search_courant_limit_o );

      // are we in the micro grid ?
      if ( meshPtr == mesh_o )
      {
        // change frame of foot from macro to micro grid
        for ( int dir = 1 ; dir <= 3 ; dir++ )
          foot[dir] -= position[dir];
      
        // use bb search for closest node, because the closest node guess provided by outer mesh (node) is irrelevant to inner one.
        UseNNS = true;
      }
      
      // search the element containing the point
      int ElemContainer = SearchElemContainingPoint(meshPtr, UseNNS, foot,
                                                    &ClosestNode, &ElemCheckedNb);
      

//      // for outer mesh
//      if ( meshPtr == mesh_o )
//        //      {
//        UseNNS = ( CourantLocal > _bb_search_courant_limit );
//      //        for ( int dir = 1 ; dir <= 3 ; dir++ )
//      //          foot[dir] += position[dir];
//      //      }
//      
//      // search the element containing the point
//      int ElemContainer = SearchElemContainingPoint(
//        meshPtr, UseNNS, foot, &ClosestNode, &ElemCheckedNb );
      
      if ( !ElemContainer && ( meshPtr == mesh_o ) )
      {
        // The foot is outside the micro grid, we search in the macro one
        // change frame of foot from micro to macro
        for ( int dir = 1 ; dir <= 3 ; dir++ )
          foot[dir] += position[dir];

        UseNNS = ( CourantLocal > _bb_search_courant_limit_o );

        // use bb search for closest node, because the closest node guess provided by inner mesh (node) is irrelevant to outer one.
        ElemContainer = SearchElemContainingPoint(
          mesh, UseNNS, foot, &ClosestNode, &ElemCheckedNb );
        
        // we'll use the macro grid from now on
        meshPtr = mesh;
      }
      
      if ( !ElemContainer )
      {
        // The point is not in the macro grid,
        // has it crossed an inlet or outlet boundary ?
        bool crossed = CheckIOCrossing(mesh, foot);
        
        // if the caracteristic crossed an inlet or outlet, it has been projected
        // onto the corresponding boundary, otherwise, we need to do sub-stepping
        // from the original position of the node, in such a way the whole
        // characteristic computed with the same sub-time-step in inside the domain.
        if ( !crossed )
        {
//            // THIS IS A HACK TO MAKE IT WORK !!!!!!!!!!!!!!!!!!!!!!!!!!!
//            // the foot is outside the outer domain and has not crossed an i/o
//            // then it should be in a wall, so we set the convected velocity to 0
//            // and jump to next node
//            for ( int dir = 1 ; dir <= 3 ; dir++ )
//              VelTildePtr[dir][node] = 0.;
//            
//            substep = false;
//            continue;
//            // END OF HACK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
          warning("%s : node %d foot outside""\n", __FUNCTION__, node );
          warning("\tfoot position : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  foot[1],foot[2],foot[3] );
          warning("\thead position : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  head[1],head[2],head[3] );
          warning("\thead velocity : "DBL_FMT" " DBL_FMT" "DBL_FMT"\n",
                  VelHead[1],VelHead[2],VelHead[3] );
          error("\tsub time step : "DBL_FMT"\n", sub_timestep );
          
          // start all over again from the original node position in the micro grid
//            for ( int dir = 1 ; dir <= 3 ; dir++ )
//            {
//              head[dir] = mesh->points[node][dir];
//              VelHead[dir] = VelExtrapol[dir][node];
//            }
//            
//            courant = mesh->courant[ node ];
//            
//            // half sub-timestep and initialize sub-time
//            sub_timestep /= 2.;
//            sub_time = sub_timestep;
//            
//            // framevel here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//            
//            // do sub-stepping
//            Substepping(level, courant, head, VelHead, FrameVel, dt,
//                        sub_timestep, foot, &CourantLocal );
//            
//            SubstepNb++;
//            
//            // we'll use the inner mesh from now on
//            meshPtr = mesh;
        }
        continue; // continues the while loop (substep)
      }
      
      //=============================================================================
      // the point is in the eulerian geometry, we can interpolate the velocity and check substepping is done
      //=============================================================================
      // The local coordinates of the foot, where we interpolate the velocity.
      double LocalPos[4] = {0.,0.,0.,0.};
      
      // get the local coordinates of the foot
      // also check that the point is in the element
      bool ElemCheck = MeshGetLocalPosition(
        meshPtr, ElemContainer, foot, LocalPos);
      
      // give info if the point is not inside the element
      if ( !ElemCheck )
      {
        info("are we in inner mesh? : " "\n", meshPtr==mesh);
        
        double *coord = mesh->points[node];
        
        info("node " INT_FMT " (" DBL_FMT "," DBL_FMT "," DBL_FMT
             ") is not in element " INT_FMT "\n",
             node, coord[1], coord[2], coord[3], ElemContainer );
        info("foot = (" DBL_FMT "," DBL_FMT "," DBL_FMT ")" "\n",
             foot[1], foot[2], foot[3]);
        error("%s : lagrangian point %u not in the element %u",
              __FUNCTION__, node, ElemContainer );
      }
      
      // is the sub-stepping done ?
//        if ( fabs( dt - sub_time ) / dt >= _min_substep_coef * 1e-2 )
      //{
//          //info("sub_time = "DBL_FMT"\n",sub_time);
//          
//          // introduced for merging spaf and spafz, we interpolate from the extrapolated velocity
//          VelInterpFromPtr = VelExtrapol;
//          const double *FrameVelPtr = FrameVel;
//          
//          // if the foot is outside the inner mesh and inside the outer mesh, then we interpolate from the outer velocity field with a zero frame velocity
//          double ZeroFrameVel[4] = {0.,0.,0.,0.};
//          
//          if ( meshPtr == mesh_o )
//          {
//            VelInterpFromPtr = Vel_o;
//            FrameVelPtr = ZeroFrameVel;
//          }
//          
//          for ( int dir = 1 ; dir <= 3 ; dir++ )
//          {
//            // Interpolate the velocity at the current foot and substep and use it for the next sub-stepping
//            VelTildePtr[dir][node] = InterpolateOrder2(
//                                                       meshPtr->ConTab[ ElemContainer ], LocalPos, VelInterpFromPtr[dir]);
//            
//            // update head position with the current foot
//            head[dir] = foot[dir];
//            
//            // update the velocity at the head of the characteristic
//            VelHead[dir] = VelTildePtr[dir][node];
//          }
//          
//          // update the local courant number from the closest node
//          courant = meshPtr->courant[ClosestNode];
//          
//          // update sub time
//          sub_time += sub_timestep;
//          
//          // framevel here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//          
//          // do sub-stepping
//          Substepping(level, courant, head, VelHead, FrameVelPtr, dt,
//                      sub_timestep, foot, &CourantLocal );
//          
//          SubstepNb++;
//          continue; // continue the while loop (substep)
//        }
      
      //=============================================================================
      // The sub-stepping is done, interpolate the velocity and store it into VelTilde*
      //=============================================================================
      VelInterpFromPtr = VelOldPtr;
      
      // when the foot is outside the inner mesh and inside the outer mesh, we interpolate from the outer velocity field
      if ( meshPtr == mesh_o )
      {
        VelInterpFromPtr = Vel_o;
        NbOfFeetInsideInnerDomain++;
      }
      
      for ( int dir = 1 ; dir <= 3 ; dir++ )
        VelTildePtr[dir][node] = InterpolateOrder2(
          meshPtr->ConTab[ ElemContainer ], LocalPos, VelInterpFromPtr[dir]);
      
      substep = false; // done with this node, get out of the while loop (substep)
    }
    
    if ( SubstepNb > SubstepNbMax ) SubstepNbMax = SubstepNb;
  }
  
  StatAdd("Fluid convection : elements checked", ElemCheckedNb);
  StatAdd("Fluid convection : maximum courant number", CourantLocalMax);
  StatAdd("Fluid convection : feet inside inner domain", NbOfFeetInsideInnerDomain);
  debug( "\t"INT_FMT" elements checked\n", ElemCheckedNb );
  debug( "\tmaximum courant number = "DBL_FMT"\n", CourantLocalMax );
  debug( "\tmaximum number of substeps = "INT_FMT"\n", SubstepNbMax );
  debug( "\t" INT_FMT " feet inside inner domain""\n",NbOfFeetInsideInnerDomain);
  
  for (int i = 1; i < 4; i++) free(VelExtrapol[i]);
}
