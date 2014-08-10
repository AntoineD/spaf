#ifndef PARTICLE_H
#define PARTICLE_H

#include "ellipsoid.h"
#include "planes.h"
#include "mesh.h"

/// Node limit on surfNodesInSp
#define NODELMT 5000

/// particle kind definition
typedef enum
{
  KIND_NONE,      ///< no kind defined, default value
  KIND_ELLIPSOID  ///< ellipsoid particle
}
ParticleKind;

/// Particle mesh structure
typedef struct
{
  int
    *ElemFull,   ///< list of elements old and new, that are inside particle
    *ElemPartly, ///< list of elements that are partly inside
    *NodesIn;    ///< list of nodes in the particle, first entry [0] is the number of nodes

//#ifdef SUB_ELEM
  int
    **ConTab; ///< connectivity of new elements created

  double
    **points; ///< coordinates of the new nodes created
//#endif
}
Parmesh_t;

/// Particle structure
typedef struct
{
  int
    ParNb; ///< the total number of particles

  ParticleKind
    kind; ///< particle kind

  ellipsoid_t
    *ellipsoid; ///< pointer to ellipsoid structure
  
  // particle properties
  double
    density,  ///< density relative to fluid density
    volume,   ///< volume
    **InertiaTensorInv,           ///< inverse of inertia tensor
    **InertiaTensorInvPhysFrame;  ///< inverse of inertia tensor in physical frame

  // kinematic part
  double
    Vel[4],       ///< Centroidal velocity of the particles at levels n+1
    Vel1[4],      ///< Centroidal velocity of the particles at levels n
    Vel2[4],      ///< Centroidal velocity of the particles at levels n-1;
    Pos[4],       ///< Centroidal coordinates at time level n+1
    Pos0[4],      ///< Centroidal coordinates relatively to the mesh
    Pos1[4],      ///< Centroidal coordinates at time level n
    Pos2[4],      ///< Centroidal coordinates at time level n-1.
    Ang[4],       ///< Euler angles at level n+1
    Ang1[4],      ///< Euler angles at level n
    Ang2[4],      ///< Euler angles at level n-1
    AngVel[4],    ///< Angular velocity at levels n+1
    AngVel1[4],   ///< Angular velocity at levels n
    AngVel2[4],   ///< Angular velocity at levels n
    Acc[4],       ///< Angular velocity at levels n
    AngAcc[4],    ///< Angular velocity at levels n
    **RotMat,     ///< rotation matrix
    Iw[4],        ///< Inertia tensor times angular velocity at step n
    Iw1[4],       ///< Inertia tensor times angular velocity at step n
    Iw2[4];       ///< Inertia tensor times angular velocity at step n-1
  
  // mesh inside particle
  Parmesh_t
    *mesh; ///< pointer to particle mesh structure
}
particle_t;

particle_t*
ReadParticlesParameters(
const char* FileName );

void
ParShiftVar(
particle_t particle[] );

#endif
