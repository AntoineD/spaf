/** \file
 Impose the rigid body motion to particles.
*/

#include "includes.h"
#include "memory.h"
#include "output.h"
#include "logging.h"
#include "linalg.h"
#include "particle_mesh_integration.h"
#include "statistics.h"
#include "particle_io.h"
#include "parse.h"
#include "pcg.h"
#include "rigid_body.h"

// file local variable
// by default we enforce rigid motion in particles explicitely
static bool _explicit = true;

// optionnaly, the angular velocity can be prescribed, this variable stores it
static bool _FixedAngVel = false;

//=============================================================================
/** Read parameters for rigid body motion. 
 */
//=============================================================================
void
ReadRigidBodyMotionParameters(
FILE* FileId )  ///< file parameters
{
  // read the way to imose rb motion
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "setting" ),
               "Missing rigid_body parameter: ""setting");
  
  info("\n""Particle rigid body motion is set");
  
  if ( TokenIs("pointwise", 3) )
  {
    _explicit = true;
    info(" pointwise""\n");
  }
  else if ( TokenIs("L2", 3) )
  {
    _explicit = false;
    info(" in L2 sense""\n");
  }
  else
    error("\n""In rigid_body_motion section, setting can pointwise or L2");
  
  // read the fixed angular velocity setting
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "fixed_angular_velocity" ),
               "Missing rigid_body parameter: ""fixed_angular_velocity");
  
  if ( TokenIs("true", 3) )
  {
    info("The angular velocity will be fixed to its initial value.");
    _FixedAngVel = true;
  }
  
  // check end of section
  TokenizeFileLine( FileId );
  assert_error( SectionIsEnd(),
               "Solver parameters : missing section = end");
}
//=============================================================================
/** Set rigid body motion in L2 sense.
 
 Solve ...
 */
//=============================================================================
static void
SetRigidBodyMotionL2(
const mesh_t*       mesh,       ///< mesh structure
const particle_t*   particle,   ///< particles structure
      FluidVar_t*   FluidVar )  ///< fluid velocity
{
  // rhs = 0
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    scal( mesh->NbOfNodes+1, 0., FluidVar->VelRHS[dir] );
  
  // compute rhs
  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
    IntegralOfRigidMotionError(
      mesh, &particle[iPar], FluidVar->Vel, FluidVar->VelRHS );
  
  // initial guess is 0.
  // Acc is used as a temporary array ot store the solution increment
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    scal( mesh->NbOfNodes+1, 0., FluidVar->Acc[dir] );

  SolveOperator("v_mass", mesh->OutsideNodes, FluidVar->VelRHS, FluidVar->Acc);
  
  // update velocity with increment
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    axpby( mesh->NbOfNodes+1, 1., FluidVar->Acc[dir], 1., FluidVar->Vel[dir]);
}
//=============================================================================
/** Set rigid body motion pointwise.
 */
//=============================================================================
static void
SetRigidBodyMotionPointwise(
const particle_t* particle, ///< particles structure
      double**    GeoTab,   ///< geometric table
      double**    Vel )     ///< fluid velocity
{
  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    // just an alias
    const particle_t* p = &particle[iPar];
    
    // impose rigid body motion explicitely and pointwise if required
    for ( int node = 1 ; node <= p->mesh->NodesIn[0] ; node++ )
    {    
      int NodeId = p->mesh->NodesIn[node];
      
      for ( int dir = 1 ; dir <= 3 ; dir++ )
      {
        int next = dir % 3 + 1,
               prev = next % 3 + 1;
        
        Vel[dir][NodeId] = p->Vel[dir] +
          p->AngVel[ next ] * ( GeoTab[NodeId][ prev ] - p->Pos0[ prev ] )
        - p->AngVel[ prev ] * ( GeoTab[NodeId][ next ] - p->Pos0[ next ] );
      }
    }
  }
}
//=============================================================================
/** Set rigid body motion.
 
 This is a wrapper that select the appropriate way of imposition.
 */
//=============================================================================
void
SetRigidBodyMotion(
const mesh_t*     mesh,       ///< mesh structure
const particle_t* particle,   ///< particles structure
      FluidVar_t* FluidVar )  ///< fluid velocity
{
  if ( _explicit )
    // pointwise
    SetRigidBodyMotionPointwise(particle, mesh->points, FluidVar->Vel);
  
  else
    // in L2 sense
    SetRigidBodyMotionL2(mesh, particle, FluidVar);
}
//=============================================================================
/** Solve the rigid body constraint.
 
 Solve the equations : see article for notations.
  
  U^{n+1} = U^p + 1 / \rho_r * ( 1 / V int_{particle} u^{**} - U^p )
 
  U^p is the predicted particle velocity
  u^{**} is the fluid velocity field from Navier-Stokes
*/
//=============================================================================
void
RigidBody(
const mesh_t*     mesh,         ///< mesh structure
const double      tau[3],       ///< time discretization coefficients
      FluidVar_t* FluidVar,     ///< fluid variables
      particle_t  particle[] )  ///< particles
{
  if ( particle == NULL ) return;
  
  debug( "\nParticles rigid body\n" );
  
  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    particle_t* p = &particle[iPar];
    
    double AverageVel[4] = {0.,0.,0.,0.},
           AverageAcc[4] = {0.,0.,0.,0.},
           AverageIw[4]  = {0.,0.,0.,0.},
           AveragedIwdt[4] = {0.,0.,0.,0.};
    
    // get average velocity, acceleration and angular momentum
    GetIntegralOverParticle(
      mesh, p, FluidVar->Acc, FluidVar->Vel,
      AverageVel, AverageAcc, AverageIw, AveragedIwdt );
    
    // get particle final velocity and angular velocity
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      double PredictedVel = ( AverageAcc[dir] -
                             tau[1] * p->Vel1[dir] -
                             tau[2] * p->Vel2[dir] ) / tau[0];
      
      p->Vel[ dir ] = PredictedVel +
        ( AverageVel[ dir ] - PredictedVel ) / p->density;
      
      p->Iw[dir] = ( AveragedIwdt[dir] -
        ( tau[1] * p->Iw1[dir]  + tau[2] * p->Iw2[dir] ) ) / tau[0];
//      p->Iw[dir] = AverageIw[dir];
    }
    
    if ( _FixedAngVel )
      copy3(p->AngVel1, p->AngVel);
    else
      // multiply by inverse inertia tensor to get particle angular velocity
      gemv3(p->InertiaTensorInvPhysFrame, p->Iw, p->AngVel);  
  }

  // Set rigid body motion
  SetRigidBodyMotion(mesh, particle, FluidVar);
  
  
  
  /*
  // The following was used to impose the rigid body motion in L2 sense by iterating. By computing the angular velocity with the inertia tensor, this process has proven to be unstable, thus no longer used.
  {
    int iRBI = 0;
    double RBIError = 0.;
    
    do
    {
      iRBI++;

      SetRigidBodyMotionL2(mesh, particle, FluidVar);
      
      RBIError = 0.;
      
      // update particle angular velocity and
      // compute RB error L2 norm to check for convergence
      for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
      {
        particle_t* p = &particle[iPar];

        error("FIX ME");
        // here is computed the angular velocity, this has not been updated since switching off the usage of the "curl" way of computing the angular velocity.
        
        double error = IntegralOfRigidMotionErrorL2( mesh, p, FluidVar->Vel);
        
        if ( error > RBIError ) RBIError = error;
      }
        
      debug("\tRBI L² error = "DBL_FMT"\n", RBIError);
    }
    while ( ( RBIError > _solver.tolerance ) && ( iRBI < _solver.IterMax ) );
    
    // update stats
    StatAdd("Particle RBI",iRBI);    
  }
   */
}
