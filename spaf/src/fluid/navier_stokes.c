/** \file
  Solve fluid flow.
 
 Here we solve the Navier-Stokes equations with a first or second order pressure correction scheme. The computational frame can follow one particle motion in some specified direction.
 */

#include "includes.h"
#include "linalg.h"
#include "memory.h"
#include "convection.h"
#include "particle_mesh_integration.h"
#include "pcg.h"
#include "version_z.h"
#include "logging.h"
#include "statistics.h"
#include "mesh_search.h"
#include "fluid_output_vtk.h"
#include "output.h"
#include "navier_stokes.h"

//=============================================================================
/** Compute the frame velocity.
 It is computed by averaging the particles velocity in the direction that is moving with the particles.
 */
//=============================================================================
static void
FrameVelSet(
const particle_t  particle[],     ///< particles
const bool        FrameVelDir[4], ///< frame fixed directions
      double      FrameVel[4] )   ///< frame velocity
{
  if ( particle == NULL ) return;
  
  // here we use the predicted velocity of the particle (how about the previous step ?)
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    if ( FrameVelDir[dir] == true )
    {
      FrameVel[ dir ] = 0.;
      
      for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
        FrameVel[ dir ] += particle[ iPar ].Vel[ dir ];
      
      FrameVel[ dir ] /= particle[1].ParNb;
    }
}
//=============================================================================
/// Find nodes of inner domain that are outside outer domain.
//=============================================================================
static void
FindNodesOutside(
const particle_t  particle[], ///< particle structure
const mesh_t*     mesh_o,     ///< outer domain mesh
      mesh_t*     mesh )      ///< mesh
{
  int OutsideNodesNb = 0;
  
  for ( int node = 1 ; node <= mesh->NbOfNodes ; node++ )
  {
    const double
    *point = mesh->points[node],
    *position = particle[1].Pos,
    NodeCoord[4] = {0., position[1]+point[1], position[2]+point[2],
                    position[3]+point[3]}; // node coordinates in mesh_o frame
    
    int temp = 0, // temp variable used for SearchElemContainingPoint
    // search an element of the outer domain that contains this node
    ElemContainer = SearchElemContainingPoint(
    mesh_o, true, NodeCoord, &temp, &temp );
    
    // if the node isn t outside, just continue
    if ( ElemContainer != 0 )
    {
      mesh->OutsideNodes[node] = false;
      continue;
    }
    
    // update outside nodes count
    OutsideNodesNb++;
    
    // add node to list
    mesh->OutsideNodes[node] = true;
  }
  
  // statistics
  StatAdd("Nodes outside outer domain", OutsideNodesNb);
  
  debug("\n\t"INT_FMT" nodes outside outer domain\n", OutsideNodesNb);
  
  // if there are nodes outside, then check the particles are fully inside, otherwise stop simulation
  if ( OutsideNodesNb == 0 ) return;
  
  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  for ( int i = 1 ; i <= particle[iPar].mesh->NodesIn[0] ; i++ )
  {
    int node = particle[iPar].mesh->NodesIn[i];
    assert_error( mesh->OutsideNodes[node] == false,
                 "Particle "INT_FMT" touched a boundary at node "INT_FMT"\n",
                 iPar,node);
  }
}
//=============================================================================
/** Solve Navier-Stokes equations.
 
 We solve the following set of equations :
 
 \f{align*}
 \frac{D \mathbf{u}}{D t} = - \nabla p + \frac{1}{Re} \nabla^2 \mathbf{u} + \frac{\rho_r-1}{Fr} \mathbf{e}_g &&
 \text{ and } &&
 \quad \nabla \cdot \mathbf{u} = 0
 \f}
 
 by sub stepping and splitting, the following assumes second order discretization in time.
 
 1 : advection-diffusion sub-step\n
 We solve for \f$\mathbf{u}^*\f$ from
 \f[
 \tau_0 \mathbf{u}^* - \frac{1}{Re} \nabla^2 \mathbf{u}^* = - \tau_1 \tilde{\mathbf{u}}^n -\tau_2  \tilde{\mathbf{u}}^{n-1} - \nabla p^n + \frac{\rho_r-1}{Fr} \mathbf{e}_g
  \f]
 with boundary conditions and where \f$\tilde{\mathbf{u}}^n\f$, \f$\tilde{\mathbf{u}}^{n-1}\f$ are the velocities from time levels \f$n\f$ and \f$n-1\f$, which have been advected alongside an approximation of the characteristics

 2 : projection sub-step\n
 We solve for \f$\mathbf{u}^{**}\f$ and \f$p\f$ from
 \f{align*}
 \tau_0 (\mathbf{u}^{**} - \mathbf{u}^*) &= - \nabla ( p^{n+1} - p^n ), \\
 \nabla \cdot \mathbf{u}^{**} &= 0, \\
 \mathbf{u}^{**} \cdot \mathbf{n} &= 0,
 \f}
 
 */
//=============================================================================
void
FluidNavierStokes(
const int         order,      ///< order of the time discretization
const mesh_t*     mesh,       ///< mesh structure
const bc_t*       bc,         ///< boundary conditions
#ifdef VERSION_Z
      mesh_t*     mesh_o,     ///< mesh structure of outer domain
      double*     Vel_o[4],   ///< fluid velocity of outer domain
#endif          
const fluid_t*    fluid,      ///< fluid parameters
const particle_t  particle[], ///< particles
const double      tau[3],     ///< time discretization coefficients
const double      dt,         ///< time step
      FluidVar_t* FluidVar )  ///< fluid variables
{
  double
  *VelTilde1[4] = {NULL,NULL,NULL,NULL}, // Convected velocity from level n 
  *VelTilde2[4] = {NULL,NULL,NULL,NULL}, // Convected velocity from level n-1
  FrameVel[4] = { 0.,0.,0.,0. }; // frame velocity
  
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    AllocVdouble(mesh->NbOfNodes, VelTilde1[dir]);
    AllocVdouble(mesh->NbOfNodes, VelTilde2[dir]);
  }
  
#ifdef VERSION_Z
  if ( UsingMicroGrid() )
  {
    // find nodes of inner domain that are outside outer domain
    FindNodesOutside(particle, mesh_o, mesh);
    
    // set those nodes velocity to 0.
    for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
      if ( mesh->OutsideNodes[NodeId] == true )
        for ( int dir = 1 ; dir <= 3 ; dir++ )
          FluidVar->Vel[dir][NodeId] = 0.;    
  }
#endif
  
  //----------------------------------------------------------------------------
  // Compute all the velocity rhs contributions
  //----------------------------------------------------------------------------
  // rhs = 0
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    scal( mesh->NbOfNodes+1, 0., FluidVar->VelRHS[dir] );
  
  // compute frame velocity
  FrameVelSet( particle, fluid->FrameVelDir, FrameVel );
  
  //----------------------------------------------------------------------------
  // Pressure gradient contribution
  //----------------------------------------------------------------------------
  // compute and substract pressure gradient to rhs
  ApplyOperator("gradient", &FluidVar->Pre, FluidVar->VelRHS);
  
  //----------------------------------------------------------------------------
  // Convection terms contribution
  //----------------------------------------------------------------------------
#ifdef VERSION_Z
  const double *ParticlePos = (order == 1) ? particle[1].Pos1 : particle[1].Pos2;
//  const double *ParticlePos = particle[1].Pos;
  
  if ( UsingMicroGrid() )
    ConvectionMicroGrid(order, mesh, dt, FrameVel,
             ParticlePos, mesh_o, Vel_o,
             FluidVar->VelOld1, FluidVar->VelOld2, VelTilde1, VelTilde2);
  else
    ConvectionMacroGrid(order, mesh, dt, FrameVel,
                        ParticlePos, mesh_o, Vel_o,
                        FluidVar->VelOld1, FluidVar->VelOld2, VelTilde1, VelTilde2);
#else
  
  Convection(order, mesh, dt, FrameVel,
             FluidVar->VelOld1, FluidVar->VelOld2, VelTilde1, VelTilde2);
#endif 
  
  // compute the convection terms in Acc
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    axpby( mesh->NbOfNodes+1, -tau[1], VelTilde1[dir], 0., FluidVar->Acc[dir] );
  
  if ( order == 2 )
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      axpby( mesh->NbOfNodes+1, -tau[2], VelTilde2[dir], 1., FluidVar->Acc[dir]);
  
  // weight by mass matrix and substract them from rhs
  ApplyOperator("convection", FluidVar->Acc, FluidVar->VelRHS);
  
  // set convected velocity at previous timestep as initial guess for conjugate gradient
  // WARNING : this has to be done before we set boundary conditions, otherwise those latter could get overwritten
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    copy(mesh->NbOfNodes+1, VelTilde1[dir], FluidVar->Vel[dir]);
  
  //----------------------------------------------------------------------------
  // Particle weight contribution
  //----------------------------------------------------------------------------
  if ( UsingMicroGrid() )
    GetParticleMomentumContribution(mesh, particle, fluid, FluidVar->VelRHS);
  
  //----------------------------------------------------------------------------
  // Boundary conditions contribution
  //----------------------------------------------------------------------------
  SetVelBC( mesh->NbOfNodes, bc, FluidVar->Vel );
  
#ifdef VERSION_Z
  // prescribe dirichlet bc by interpolating outer domain velocity, this has to be done after SetVelBC !
  if ( UsingMicroGrid() )
    MassConserveBC( particle[1].Pos, mesh_o, Vel_o, mesh, FluidVar->Vel );
#endif
  
  // weight with stifness matrix and substract to rhs
  ApplyOperator("v_bc", FluidVar->Vel, FluidVar->VelRHS);
  
  //----------------------------------------------------------------------------
  // Solve for the velocity: advection-diffusion step
  //----------------------------------------------------------------------------
  debug( "\nVelocity diffusion step\n" );
  
  SolveOperator("v_stiffness", mesh->OutsideNodes, FluidVar->VelRHS, FluidVar->Vel);

  //----------------------------------------------------------------------------
  // Solve for the pressure star: projection step
  //----------------------------------------------------------------------------
  debug( "\nPressure prediction step\n" );
  
  // initialize pressure rhs
  double *PreStar = NULL, // Predicted presssure
  *PreRHS  = NULL; // Presssure RHS
  
  AllocVdouble( mesh->NbOfPressureNodes,PreStar );
  AllocVdouble( mesh->NbOfPressureNodes,PreRHS );
  
  // compute and add velocity divergence to pressure rhs
  ApplyOperator("divergence", FluidVar->Vel, &PreRHS);
  
  // mutliply pressure rhs by -tau_0
  scal( mesh->NbOfPressureNodes+1, -tau[0], PreRHS );
  
  // set previous step pressure as initial guess for conjugate gradient
  copy( mesh->NbOfFreePressureNodes+1, FluidVar->Pre, PreStar );
  
  SolveOperator("p_stiffness", mesh->OutsideNodes, &PreRHS, &PreStar);

  //----------------------------------------------------------------------------
  // Solve for the velocity: projection step
  //----------------------------------------------------------------------------
  debug( "\nVelocity projection step\n" );
  
  // rhs = 0
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    scal( mesh->NbOfNodes+1, 0., FluidVar->VelRHS[dir] );
  
  // compute and substract pressure star gradient to velocity rhs
  ApplyOperator("gradient", &PreStar, FluidVar->VelRHS);
  
  // VelTilde2 = 0, used here as temporary array, it will store the un-weigthed gradient
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    scal( mesh->NbOfNodes+1, 0., VelTilde2[dir] );
  
  // compute the un-weigthed gradient
  SolveOperator("v_mass", mesh->OutsideNodes, FluidVar->VelRHS, VelTilde2);
  
  // Remove the non solenoidal part of the velocity, i.e. the un-weigthed gradient
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    axpby(mesh->NbOfNodes+1, 1. / tau[0], VelTilde2[dir], 1.,
          FluidVar->Vel[dir]);
  
  //----------------------------------------------------------------------------
  // Solve for the pressure: correction step
  //----------------------------------------------------------------------------
  debug( "\nPressure correction step\n" );
  
  // p = p + p*
  axpby( mesh->NbOfPressureNodes+1, 1., PreStar, 1., FluidVar->Pre );
  
#ifdef PRE_INCREMENTAL_ROTATIONAL
  // then solve for - 1/Re Div(vel) = PreMass^{-1} * PreRHS, solve in p* again
  
  // solve for PreStar
  for ( int i = 1 ; i <= mesh->NbOfPressureNodes ; i++ ) 
    PreStar[i] = PreRHS[i] * FluidOperators->PreMassPrec[i];
  
  // p = p + p* = p - 1/Re Div(vel)
  axpby( mesh->NbOfPressureNodes+1, 1., PreStar, 1., FluidVar->Pre );
#endif
  
  // Compute the acceleration field of the fluid, Acc = tau0 * Vel - Acc
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    axpby(mesh->NbOfNodes+1, tau[0], FluidVar->Vel[dir], -1.,
          FluidVar->Acc[dir]);
  
  // free local arrays
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {  
    free(VelTilde1[dir]);    
    free(VelTilde2[dir]);
  }
  free(PreStar);
  free(PreRHS);  
}
