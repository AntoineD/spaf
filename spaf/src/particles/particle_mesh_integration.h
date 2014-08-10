#ifndef _PARTICLE_DOMAIN_INTEGRATION_H_
#define _PARTICLE_DOMAIN_INTEGRATION_H_

#include "gauss.h"
#include "particle_mesh.h"
#include "particle.h"
#include "fluid.h"
#include "mesh.h"

void
IntegralOfRigidMotionError(
const mesh_t*     mesh,
const particle_t* particle,
      double**      Vel,
      double**      rhs );

double
IntegralOfRigidMotionErrorL2(
const mesh_t*     mesh,
const particle_t* particle,
      double**      Vel );

void
GetParticleMomentumContribution(
const mesh_t*     mesh,    
const particle_t* particle,
const fluid_t*    fluid,   
      double**    rhs );   

void
GetIntegralOverParticle(
const mesh_t*     mesh,    
const particle_t* particle,
      double**    Acc,           
      double**    Vel,           
      double      AverageVel[4], 
      double      AverageAcc[4], 
      double      AverageIw[4],
      double      AveragedIwdt[4] );

void
GetFluidAngularMomentum(
const mesh_t*     mesh,       
const particle_t* particle,   
      double**    Vel,        
      double      AngVel[4] );

#endif
