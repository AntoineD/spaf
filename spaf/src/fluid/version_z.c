/** \file
  Stuff for the 2 domains version of spaf.
 */

#include "version_z.h"
#include "logging.h"
#include "shape_functions.h"
#include "mesh_search.h"
#include "linalg.h"

//=============================================================================
/** Set particle initial conditions.
 
 Mesh relative position is set to 0, i.e. the particle is fixed in the inner mesh, at origin.
 Velocity at previous step is set to outer fluid velocity at particle centroid.
 */
//=============================================================================
void
ParticleInitialConditions(
const mesh_t*     mesh_o,
const FluidVar_t* FluidVar_o,
      particle_t* particle )
{
  info("Setting particle initial conditions for version z""\n");
  
  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    scal3(0.,particle[iPar].Pos0);
    
    // find outer element
    int temp = 0,
    ElemContainer_o = SearchElemContainingPoint(
      mesh_o, true, particle[iPar].Pos1, &temp, &temp );

    // check the node is not outside the outer mesh
    assert_error( ElemContainer_o != 0,
      "Particle centroid is not in outer mesh" );
    
    // local coordinate of the point
    double LocalCoord_o[4] = {0., 0., 0., 0.};
    
    // get the local coordinates
    bool ElemContainerFound = MeshGetLocalPosition(
      mesh_o, ElemContainer_o, particle[iPar].Pos1, LocalCoord_o );
    
    assert_error( ElemContainerFound,
      "Particle centroid is not in the element %u", ElemContainer_o );
    
    // interpolate the velocity from the outer element
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      particle[iPar].Vel1[dir] = InterpolateOrder2(
        mesh_o->ConTab[ ElemContainer_o ], LocalCoord_o, FluidVar_o->Vel[dir]);

    // some infos
    info("\t""particle %d initial velocity = %g %g %g""\n", iPar,
         particle[iPar].Vel1[1],
         particle[iPar].Vel1[2],
         particle[iPar].Vel1[3]);    
  }
}
