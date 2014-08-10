#include "includes.h"
#include "memory.h"
#include "logging.h"
#include "linalg.h"
#include "particle_mesh.h"
#include "collision.h"
#include "particle_io.h"

#include "prediction_correction.h"

//==============================================================================
/** Update the rotation matrix.
 */
//==============================================================================
static void
UpdateOrientation(
const double    RotationVectorIncrement[4],
const double    RotationVectorOld[4],
      double**  RotationMatrix,
      double    RotationVector[4] )
{
#if 1
  // compute the rotation matrix associated with the
  // increment rotation vector
  double **RotationMatrixIncrement = NULL;
  AllocMdouble(3,3,RotationMatrixIncrement);
  RotationVectorToMatrix(RotationVectorIncrement, RotationMatrixIncrement);
  
  // compute the rotation matrix associated with the
  // old rotation vector
  double **RotationMatrixOld = NULL;
  AllocMdouble(3,3,RotationMatrixOld);
  RotationVectorToMatrix(RotationVectorOld, RotationMatrixOld);
  
  // compute new rotation matrix
  gemm3(RotationMatrixIncrement, RotationMatrixOld, RotationMatrix);
  FreeM(RotationMatrixIncrement);
  
  // compute rotation vector to new orientation
  RotationMatrixToVector(RotationMatrix, RotationVector);
  FreeM(RotationMatrixOld);
#else
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    RotationVector[ dir ] = RotationVectorOld[ dir ] + RotationVectorIncrement[ dir ];

  RotationVectorToMatrix(RotationVector, RotationMatrix);
#endif
}
//==============================================================================
/** Prediction of the particle's position and orientation

  If collisions with force are managed, ParColForcePred is called,
  otherwise the positions and orientations are updated first order :
 
  \f$ X_i^{n+1} = X_i^n + U_i^n \Delta t \f$ and
  \f$ \theta_i^{n+1} = theta_i^n + \omega_i^n \Delta t \f$.

  or second order scheme :
  \f$ X_i^{n+1} = X_i^{n-1} + 2 \Delta t U_i^n \f$ and
  \f$ \theta_i^{n+1} = theta_i^{n-1} + 2 \Delta t \omega_i^n \f$
*/
//==============================================================================
void
PredictParticlePosition(
const int      order,        ///< order of time integration
const mesh_t*     mesh,         ///< mesh structure
const double      dt,           ///< time step 
const fluid_t*    fluid,        ///< fluid parameters
      particle_t  particle[] )  ///< particle structure
{
  if ( particle == NULL ) return;

  debug( "\nParticles prediciton step\n" );

  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    particle_t *p = &particle[iPar];
    
    // rotation vector from a previous orientation (n-1 or n-2, according
    // to order ) to the new one
    double RotationVectorIncrement[4] = {0.,0.,0.,0.},
    // previous rotation vector (n-1 or n-2, according to order ) 
          *RotationVectorOld = p->Ang1;

    // Correction of particles position and orientation
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      if ( order == 1 )
      {
        p->Pos[ dir ] = p->Pos1[ dir ] + dt * p->Vel1[ dir ];
        RotationVectorIncrement[dir] = dt * p->AngVel1[ dir ];
      }
      else if ( order == 2 )
      {
        p->Pos[ dir ] = p->Pos2[ dir ] + 2. * dt * p->Vel1[ dir ];
        RotationVectorIncrement[dir] = 2. * dt * p->AngVel1[ dir ];
      }
      else
        error("bad order");

      // update particle-mesh relative position
      if ( fluid->FrameVelDir[dir] == false )
        p->Pos0[ dir ] = p->Pos[ dir ];
    }
    
    if ( order == 2 ) RotationVectorOld = p->Ang2;
    
    // update the rotation matrix
    UpdateOrientation(
      RotationVectorIncrement, RotationVectorOld, p->RotMat, p->Ang );

    // update the rotation matrix
    ChangeMatrixBase(
      p->InertiaTensorInv, p->RotMat, p->InertiaTensorInvPhysFrame );
  }

  // Detect collision and stop computation if required
  DetectCollision( particle );

  // Compute particle mesh, this is needed before solving the flow.
  GetParticleMesh( mesh, particle );
}



//==============================================================================
/** Correction of the particle's position and orientation

  If collisions with force are managed, ParColForceCor is called,
  otherwise the positions and orientations are updated by first order :

  \f$ X_i^{n+1} = X_i^n + U_i^{n+1} \Delta t \f$ and
  \f$ \theta_i^{n+1} = theta_i^n + \omega_i^{n+1} \Delta t \f$.

  or second order scheme :
  \f$ X_i^{n+1} = X_i^n + \Delta t ( U_i^{n+1} + U_i^n ) / 2 \f$ and
  \f$ \theta_i^{n+1} = theta_i^n + \Delta t ( \omega_i^{n+1} + \omega_i^n ) / 2 \f$
*/
//==============================================================================
void
CorrectParticlePosition(
const int      order,        ///< order of time integration
const double      dt,           ///< time step
const double      tau[3],       ///< time discretization coefficients
const fluid_t*    fluid,        ///< fluid structure
      particle_t  particle[] )  ///< particle structure
{
  if ( particle == NULL ) return;

  debug( "\nParticles correction step\n" );

  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    particle_t *p = &particle[iPar];
    
    // vector of rotation from a previous orientation to new one
    double RotationVectorIncrement[4] = {0.,0.,0.,0.},
    // previous votation vector
    *RotationVectorOld = p->Ang1;
    
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      if ( order == 1 )
      {
        p->Pos[ dir ] = p->Pos1[ dir ] + dt * p->Vel[ dir ];
        RotationVectorIncrement[dir] = dt * p->AngVel[ dir ];
      }
      else if ( order == 2 )
      {
        p->Pos[ dir ] = p->Pos1[ dir ] + .5 * dt * ( p->Vel[ dir ]    + p->Vel1[ dir ] );
        RotationVectorIncrement[dir] = .5 * dt * ( p->AngVel[ dir ] + p->AngVel1[ dir ] );
      }
      else
        error("%s bad order", __FUNCTION__ );

      // update particle-mesh relative position
      if ( fluid->FrameVelDir[dir] == false )
        p->Pos0[ dir ] = p->Pos[ dir ];
      
      // get centroid acceleration
      p->Acc[dir] = tau[0] * p->Vel[dir] + tau[1] * p->Vel1[dir] + tau[2] * p->Vel2[dir];
      
      // get angular acceleration
      p->AngAcc[dir] = tau[0] * p->AngVel[dir] + tau[1] * p->AngVel1[dir] + tau[2] * p->AngVel2[dir];
    }

    // update the rotation matrix
    UpdateOrientation(
      RotationVectorIncrement, RotationVectorOld, p->RotMat, p->Ang );
    
    // update the rotation matrix
    ChangeMatrixBase(
      p->InertiaTensorInv, p->RotMat, p->InertiaTensorInvPhysFrame );
  }

  // Detect collision and stop computation if required
  DetectCollision( particle );
}
