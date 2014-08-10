/**
 \file particle_mesh_integration.c 
 
 \brief Compute integrals over particle domains.
 
 The integrals are computed using Gauss cubature, a particle domain has been filled with original elements from the mesh and sub-elements that match the particle shape. The main routine that performs the integration is the marco integrate_this, which prepare everything to account for the contribution of a gauss point in an element to the integral considered, thus only the computing of this contributon has to be provided to this macro.
 */

#include "includes.h"
#include "memory.h"
#include "linalg.h"
#include "logging.h"
#include "particle_mesh_find_elem.h"
#include "mesh_search.h"
#include "shape_functions.h"
#include "particle_mesh_integration.h"


// local file wise variables
/// Scaling factor to possibly overcome round off error summation in integrals.
static double _ScalingCoef = 1.;


//=============================================================================
/** This macro prepare the intregration over a particle domain.
 
 For each element inside a particle, original or from a subdivision, and at all gauss points of this element, we compute the global position, the value of all shape functions and the value of the gauss weight times the element jacobian. 
 */
//=============================================================================

#define integrate_this( code ) \
double \
/* volume of the particle mesh */ \
volume = 0., \
/* weighted jacobian for integration */ \
*wdetJ = NULL, \
/* positions of the gauss points */ \
**GaussPoints = NULL, \
/* shape functions values at the gauss points */ \
**phi = NULL; \
\
gauss_t* gauss = mesh->gauss; \
\
AllocVdouble( gauss->NbOfPoints, wdetJ); \
AllocMdouble( gauss->NbOfPoints, 3, GaussPoints ); \
AllocMdouble( gauss->NbOfPoints, NODES_PER_EL, phi ); \
\
/* loop over the elements in the particle */ \
for ( int iElem = 1 ; iElem <= particle->mesh->ElemFull[0] ; iElem++ ) \
{ \
  /* index of current element */ \
  int ElemId = particle->mesh->ElemFull[ iElem ]; \
   \
  /* Check whether this element is in the original mesh */ \
  if ( ElemId <= mesh->NbOfElements ) \
  { \
    /* Compute weighted jacobians and gauss point positions */ \
    PrepareIntegrationForOriginalElement( mesh, ElemId, wdetJ, GaussPoints ); \
     \
    /* Loop over all gauss points and integrate with quadrature */ \
    for ( int iGp = 1 ; iGp <= gauss->NbOfPoints ; iGp++ ) \
    { \
      /* some pointers */ \
      double \
      *ShapeFunctionsAtGP = gauss->phi[iGp], \
      *GPPosition = GaussPoints[iGp], \
      WeightedJacobianAtGP = wdetJ[iGp]; \
       \
      int *ElemConnectivity = mesh->ConTab[ElemId]; \
       \
      volume += wdetJ[ iGp ]; \
       \
      /* here is done the integrand computation that is provided to the macro */ \
      do{ code } while(0); \
    \
  } \
} \
else /* This element is a new one */ \
{  \
  /* Compute weighted jacobians, gauss point positions, shape functions */ \
  /* at those points and get the original element that contains this one */ \
  int ParentElem = PrepareIntegrationForSubElement( \
    mesh, ElemId, particle->mesh, wdetJ, GaussPoints, phi ); \
   \
  /* skip if element is too small */ \
  if ( ParentElem == 0 ) continue; \
   \
  /* Loop over all gauss points and integrate with quadratute */ \
  for ( int iGp = 1 ; iGp <= gauss->NbOfPoints ; iGp++ ) \
  { \
    double \
    *ShapeFunctionsAtGP = phi[iGp], \
    *GPPosition = GaussPoints[iGp], \
    WeightedJacobianAtGP = wdetJ[iGp]; \
     \
    int *ElemConnectivity = mesh->ConTab[ParentElem]; \
     \
    volume += wdetJ[ iGp ]; \
     \
    do{ code } while(0); \
  } \
} \
} \
 \
free(wdetJ); \
FreeM(GaussPoints); \
FreeM(phi);



//=============================================================================
/** Prepare integration ingredients for an original element.
 */
//=============================================================================
static void
PrepareIntegrationForOriginalElement(
const mesh_t*   mesh,         ///< mesh structure
const int    ElemId,       ///< current element index
      double*     wdetJ,        ///< weighted jacobian at all gauss points
      double**    GaussPoints ) ///< position of all gauss points
{
  // Get global position of the gauss points
  for ( int iGp = 1 ; iGp <= mesh->gauss->NbOfPoints ; iGp++ )
  {
    GetJacobianAtPoint(
      NODES_PER_EL, mesh->points, mesh->ConTab[ElemId],
      mesh->gauss->local_dphi[iGp], mesh->gauss->weight[iGp],
      NULL, &wdetJ[iGp]);
    
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      GaussPoints[iGp][dir] = 0.;

    for ( int LocalNodeId = 1 ; LocalNodeId <= NODES_PER_EL ; LocalNodeId++ )
    {
      int NodeId = mesh->ConTab[ ElemId ][ LocalNodeId ];

      for ( int dir = 1; dir <= 3; dir++ )
        GaussPoints[iGp][ dir ] += mesh->points[ NodeId ][ dir ] *
                                   mesh->gauss->phi[ iGp ][ LocalNodeId ];
    }
  }
}
//=============================================================================
/** Prepare integration ingredients for a new element.
 */
//=============================================================================
static int
PrepareIntegrationForSubElement(
const mesh_t*     mesh,         ///< mesh structure
const int         ElemId,       ///< current element index
const Parmesh_t*  ParMesh,      ///< a particle mesh structure
      double*     wdetJ,        ///< weighted jacobian at all gauss points
      double**    GaussPoints,  ///< position of all gauss points
      double**    PhiAtGP )     ///< shape functions values at all gauss points
{
  // Get jacobian for this new element
  bool integrate = GetParElementJacobian( mesh, ParMesh, ElemId, wdetJ );

	// if the element volume is negligible, return 0 means do not integrate
  if ( integrate == false ) return 0;
  
  int ParentElem = ParMesh->ConTab[ 0 ][ ElemId - mesh->NbOfElements ];

  //Loop over Gauss points
  for ( int iGp = 1 ; iGp <= mesh->gauss->NbOfPoints ; iGp++ )
  {
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      GaussPoints[iGp][dir] = 0.;

    // compute gauss point global position
    for ( int LocalNodeId = 1; LocalNodeId <= NODES_PER_EL; LocalNodeId++ )
    {
      int NodeId = ParMesh->ConTab[ LocalNodeId ][ ElemId - mesh->NbOfElements ];
      
      if ( NodeId <= mesh->NbOfNodes )
        for ( int dir = 1; dir <= 3; dir++ )
          GaussPoints[iGp][ dir ] += mesh->points[ NodeId ][ dir ] *
          mesh->gauss->phi[ iGp ][ LocalNodeId ];
      else
        for ( int dir = 1; dir <= 3; dir++ )
          GaussPoints[iGp][ dir ] += ParMesh->points[ dir ][ NodeId - mesh->NbOfNodes ] *
          mesh->gauss->phi[ iGp ][ LocalNodeId ];
    }
    
    double LocalPos[4] = {0.,0.,0.,0.};
    
    // Find local coordinates of this gp in parent element
    bool Elemfound = MeshGetLocalPosition(mesh, ParentElem, GaussPoints[iGp],
                                          LocalPos);
    
    assert_error( Elemfound, "The gauss point is not in the original element" );
    
    GetShapeFunctionsAtPoint(5, &LocalPos[1], &PhiAtGP[iGp][1]);
  }
  
  return ParentElem;
}
//=============================================================================
/** Add a contribution of an integral to rhs.
 */
//=============================================================================
inline static void
RHSUpdate(
const double  integrand,  ///< the contribution to the integral
const double* phi,        ///< shape functions values at current gauss point
const double  wdetJ,      ///< weighted jacobian at current gauss point
const int*    ElemCon,    ///< local connectivity of current element
      double* rhs )       ///< right hand side
{
  for ( int LocalNodeId = 1 ; LocalNodeId <= NODES_PER_EL ; LocalNodeId++ )
    rhs[ ElemCon[ LocalNodeId ] ] += wdetJ * phi[ LocalNodeId ] * integrand;
}
//=============================================================================
/** Compute the velocity at current gauss point.
 */
//=============================================================================
inline static double
VelAtGaussPoint(
const int* ElemCon,  ///< local connectivity of current element
const double*   phi,      ///< shape functions values at current gauss point
      double*   Vel )     ///< fluid velocity
{
  double VelAtGP = 0.;
  
  for ( int LocalNodeId = 1 ; LocalNodeId <= NODES_PER_EL ; LocalNodeId++ )
    VelAtGP += Vel[ ElemCon[ LocalNodeId ] ] * phi[ LocalNodeId ];
  
  return VelAtGP;
}
//=============================================================================
/** Computes the dir component of the rigid body motion error.
 
 \f$ error = U^{n+1} + \omega \times ( x - X^p ) - u \f$
 */
//=============================================================================
inline static double
RigidMotionError(
const int      dir,        ///< direction
const int*       ElemCon,    ///< local connectivity of current element
const double*     phi,        ///< shape functions values at current gauss point
const double      GpPos[4],   ///< position of current gauss point
      double**    Vel,        ///< fluid velocity
const particle_t* particle )  ///< a particle structure
{
  return VectorProduct(dir,particle->AngVel,GpPos,particle->Pos0)
  + particle->Vel[ dir ] - VelAtGaussPoint(ElemCon, phi, Vel[dir]);
}
//=============================================================================
/** Integrand for rigid motion error.
 */
//=============================================================================
void
IntegralOfRigidMotionError(
const mesh_t*     mesh,     ///< mesh structure
const particle_t* particle, ///< particles
      double**    Vel,      ///< fluid velocity
      double**    rhs )     ///< right hand side
{  
  integrate_this
  (        
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      double integrand = RigidMotionError(
       dir,ElemConnectivity,ShapeFunctionsAtGP,GPPosition,Vel,particle);

      RHSUpdate( integrand, ShapeFunctionsAtGP, WeightedJacobianAtGP,
                ElemConnectivity, rhs[ dir ] );
    }
  )
}
//=============================================================================
/** Compute L2 norm of the rigid body error.
 */
//=============================================================================
double
IntegralOfRigidMotionErrorL2(
const mesh_t*     mesh,     ///< mesh structure
const particle_t* particle, ///< a particle structure
      double**      Vel )   ///< fluid velocity
{
  double L2error = 0.;
  
  integrate_this
  (        
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
     double integrand = RigidMotionError(
       dir,ElemConnectivity,ShapeFunctionsAtGP,GPPosition,Vel,particle);
     
     L2error += integrand * integrand * WeightedJacobianAtGP;
    }
  )
  
  return sqrt(L2error);
}
//=============================================================================
/** Compute the gravity term for bouyant particle.
 */
//=============================================================================
inline static void
MomentumPrediction(
const int*        ElemConnectivity,     ///< local connectivity of current element
const double      WeightedJacobianAtGP, ///< weighted jacobian at current gauss point
const double*     ShapeFunctionsAtGP,   ///< shape functions values at current gauss point
const particle_t* p,      ///< pointer to a particle
const fluid_t*    fluid,  ///< fluid parameters
      double**    rhs )   ///< position of the gauss point
{
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {     
    // gravity
    double integrand = fluid->Gravity[dir] / fluid->Fr;

    // weight contribution
    integrand *= p->density - 1.;
    
    RHSUpdate( integrand, ShapeFunctionsAtGP, WeightedJacobianAtGP,
              ElemConnectivity, rhs[ dir ] );
  }
}
//=============================================================================
/** Compute the gravity term for bouyant particle.
 */
//=============================================================================
static void
GetParticleMomentumContributionWrapped(
const mesh_t*     mesh,     ///< mesh structure
const particle_t* particle, ///< a particle structure
const fluid_t*    fluid,    ///< fluid structure
      double**    rhs )     ///< rhs to be updated with gravity term
{  
  integrate_this
  (
    MomentumPrediction(
      ElemConnectivity, WeightedJacobianAtGP, ShapeFunctionsAtGP,
      particle, fluid, rhs);
  )
}
//=============================================================================
/** Compute the gravity contribution for all bouyant particles.
 */
//=============================================================================
void
GetParticleMomentumContribution(
const mesh_t*     mesh,     ///< mesh structure
const particle_t* particle, ///< particles structure
const fluid_t*    fluid,    ///< fluid structure
      double**    rhs )     ///< rhs to be updated with gravity terms
{
  if ( particle == NULL ) return;
  
  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    if ( particle[iPar].density == 1. ) continue;
    
    GetParticleMomentumContributionWrapped( mesh, &particle[iPar], fluid, rhs );    
  }
}
//=============================================================================
/** Compute the angular velocity.
 
 Compute the contribution of an element to the integral giving the
 angular velocity.
 */
//=============================================================================
inline static void
IntegralOverParticle(
const int*      ElemConnectivity,     ///< local connectivity of current element
const double    WeightedJacobianAtGP, ///< weighted jacobian at current gauss point
const double*   ShapeFunctionsAtGP,   ///< shape functions values at current gauss point
const double    GPPosition[4],        ///< position of the gauss point
const double    ParPos0[4],           ///< position of the particel's centroid
      double**  Acc,                  ///< fluid acceleration
      double**  Vel,                  ///< fluid velocity
      double    AverageVel[4],        ///< average of fluid velocity over particle
      double    AverageAcc[4],        ///< average of fluid acceleration over particle
      double    AverageIw[4],         ///< average of fluid angular momentum over particle
      double    AveragedIwdt[4] )     ///< average of fluid time derivative of angular momentum over particle
{
  double GpAcc[4] = {0., 0., 0., 0.},
         GpVel[4] = {0., 0., 0., 0.};
  
  // get velocity and acceleration at gauss point and gauss point position relative to particle center
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    GpVel[dir] = VelAtGaussPoint(ElemConnectivity, ShapeFunctionsAtGP, Vel[dir]);
    GpAcc[dir] = VelAtGaussPoint(ElemConnectivity, ShapeFunctionsAtGP, Acc[dir]);
  }
  
  // Compute the integrand, i.e.  \f$ ( x - X^p ) \times u \f$ 
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    // minus because of the way VectorProduct is defined, scaling because integrand can be very small, and round off error can accumulate
    AverageIw[ dir ] -= WeightedJacobianAtGP * _ScalingCoef * 
      VectorProduct(dir,GpVel,GPPosition,ParPos0);
    
    AveragedIwdt[ dir ] -= WeightedJacobianAtGP * _ScalingCoef * 
      VectorProduct(dir,GpAcc,GPPosition,ParPos0);
    
    AverageVel[ dir ] += WeightedJacobianAtGP * GpVel[dir];
    
    AverageAcc[ dir ] += WeightedJacobianAtGP * GpAcc[dir];
  }  
}
//=============================================================================
/** Compute integrals over a particle.
 
 We compute the average velocity and acceleration of the fluid over a particle. For the angular momentum, we compute the average time the exact particle volume.
 */
//=============================================================================
void
GetIntegralOverParticle(
const mesh_t*     mesh,          ///< mesh structure
const particle_t* particle,      ///< pointer to a particle
      double**    Acc,           ///< fluid acceleration
      double**    Vel,           ///< fluid velocity
      double      AverageVel[4], ///< average of fluid velocity over particle
      double      AverageAcc[4], ///< average of fluid acceleration over particle
      double      AverageIw[4], ///< average of fluid angular momentum over particle
      double      AveragedIwdt[4] ) ///< average of fluid time derivative of angular momentum over particle

{
  // initialize
  scal3(0., AverageVel);
  scal3(0., AverageAcc);
  scal3(0., AverageIw);
  scal3(0., AveragedIwdt);

  integrate_this
  (
    IntegralOverParticle(
      ElemConnectivity, WeightedJacobianAtGP, ShapeFunctionsAtGP, GPPosition,
      particle->Pos0, Acc, Vel, AverageVel, AverageAcc, AverageIw, AveragedIwdt);
  )
  
  // divide by particle volume
  scal3( 1./volume, AverageVel );
  scal3( 1./volume, AverageAcc );
  scal3( 1./volume, AverageIw );
  scal3( 1./volume, AveragedIwdt );
  
  // multiply back with scaling
  scal3( 1./_ScalingCoef, AverageIw );
  scal3( 1./_ScalingCoef, AveragedIwdt );

  // multiply by ideal particle volume, as the inertia tensor includes it
  scal3( particle->volume, AverageIw );
  scal3( particle->volume, AveragedIwdt );
}
