/** \file
  Manage particle-wall collision.
 */

#include "memory.h"
#include "logging.h"
#include "parse.h"
#include "linalg.h"
#include "particle.h"
#include "planes.h"
#include "mesh.h"
#include "collision.h"

// local file-wise parameters for collision

// flags for detection on or off
static bool _DetectCollision = false;
// limit distance below which collision is considered to be true
static double _threshold = 0.;
// to decide what to do when a collision occurs
static bool _Repulse = false;
// number of step of the repulsive force algorithm
static int _SubStepNb = 0;
/// definition of planes for collision detection
static planes_t* _planes = NULL;


//=============================================================================
/** Read collision planes defintions. 
 */
//=============================================================================
void
GetCollisionPlanes(
const char* FileName )  ///< name of the file that contains planes definitions
{
  _planes = PlanesRead(FileName);
}
//=============================================================================
/** Read parameters for particle collision. 
 */
//=============================================================================
void
ReadCollisionParameters(
FILE* FileId )  ///< file id that contains the collision parameters
{  
  // read the threshold, i.e. the distance below which collision is considered as true
  TokenizeFileLine(FileId);
  assert_error( TokenNameIs( "threshold" ),
    "Missing collision parameter : ""threshold");
  
  _threshold = Token2double( 3, '[', 0., 1.e2,'[' );
  
  // repulsive force is not implemented yet
  _Repulse = false;
      
  // read the action type, i.e. what to do when a collision occurs
//  TokenizeFileLine(FileId);
//  assert_error( TokenNameIs( "action_type" ),
//    "Missing collision parameter : ""action_type");

  // repulsive force not yet implemented, so stop is the only option, this is the default _Repulse
  
//  if ( TokenValueIs( "stop" ) ) _Repulse = false; // stop the computation
      
// not yet implemented
// else if ( TokenValueIs( "force" ) ) _Repulse = true; // stop the computation

//  else error( "bad action type token : %s", Token(3) );
      
// not yet implemented
// read the sub-step number for collision avoiding iterations
// else if ( TokenNameIs( "sub_step_number" ) )
//   _SubStepNb = Token2int( 3, ']', 0, 100,'[' );

  // check section = end
  TokenizeFileLine(FileId);
  assert_error( SectionIsEnd(),
    "Collision parameters : missing section = end");
  
  _DetectCollision = true;
  
  // print some infos
  info( "\nParticles Collision parameters :\n" );
  
  info( "\tThreshold = %g\n", _threshold );
  
  if ( _Repulse == true )
  {
    info( "\tRepulsive force will be applied\n" );
    info( "\tNumber of sub-steps = "INT_FMT"\n", _SubStepNb );    
  }
  else
    info( "\tStop simulation on collision.\n" );
}
/*==============================================================================
Purpose:
  dispatch distance computation to corresponding routine
==============================================================================*/
static double
GetParToPlaneDistance(
const particle_t* particle,
const equation_t* PlaneEquation )
{
  double distance;
  
  // select the particle belonging criterion
  switch ( particle->kind )
  {
  case KIND_ELLIPSOID :
    distance = EllipsoidPlaneDist(
      particle->ellipsoid, particle->Pos, particle->RotMat, PlaneEquation );
    break;
  
  default :
    error( "bad particle kind" );
  }

  return distance;
}
/*==============================================================================
PURPOSE:
  compute the minimum particle-walls distance and check collision
==============================================================================*/
void
DetectCollision(
particle_t  particle[] )  ///< particle structure
{
  if ( _DetectCollision == false ) return;
  
  double distance_min = DBL_MAX;
 
  // checking
  assert_error( _planes != NULL, "No plane specified" );
  assert_error( _planes->n > 0, "No plane specified" );

  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  for ( int iPlane = 1 ; iPlane <= _planes->n ; iPlane++ )
  {
    // get particle plane distance
    double distance = GetParToPlaneDistance( &particle[iPar], &_planes->equation[iPlane] );
  
//    printf("%g\n", distance);
    
    // find minimum distance
    if ( distance < distance_min ) distance_min = distance;
  }

//  printf("%g %g %d\n", distance_min, ParColParam->threshold, distance_min <= ParColParam->threshold);
  
  if ( distance_min <= _threshold ) error( "Collision !""\n");
}


/*==============================================================================
Some common vairables for the next routines
==============================================================================*/
static bool FirstTime = true;

static double **PosCol , **VelCol1, **VelColCor1,
            **PosCol1, **VelCol2, **VelColCor2,
            **VelCol;

static double dts = 0.; // substep size
/*==============================================================================
PURPOSE
  take care of particle collision by applying a repulsive force.
==============================================================================*/
void
ParColForcePred(
const double        dt,
const fluid_t*    fluid,
      particle_t  particle[] )
{
  error("check deeply before any use");
  // check the moving frame stuff

  if ( FirstTime == true )
  {
    AllocMdouble( particle[1].ParNb , 3, PosCol     );
    AllocMdouble( particle[1].ParNb , 3, PosCol1    );
    AllocMdouble( particle[1].ParNb , 3, VelCol     );
    AllocMdouble( particle[1].ParNb , 3, VelCol1    );
    AllocMdouble( particle[1].ParNb , 3, VelCol2    );
    AllocMdouble( particle[1].ParNb , 3, VelColCor1 );
    AllocMdouble( particle[1].ParNb , 3, VelColCor2 );
    
    FirstTime = false;
  }

  // collision sub iterations timestep
  dts = dt / _SubStepNb;

  // PosCol1 = Pos1
  // VelCol1 = Vel1
  int dir, iPar;

  for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    copy3( particle[ iPar ].Pos1, PosCol1[ iPar ] ); // check Pos1 and frame change
    copy3( particle[ iPar ].Vel1, VelCol1[ iPar ] );
    RotationVectorToMatrix(particle[ iPar ].Ang, particle[ iPar ].RotMat );
  }

  //=============================================================================
// prediction of particles collisions
//=============================================================================
  for ( int i = 1 ; i <= _SubStepNb ; i++ )
  {
////////////////////////////////////////////////////////////////////////////////
//    ApplyForce( fluid, dts, PosCol1, WallPlanes, VelColCor1 );
////////////////////////////////////////////////////////////////////////////////
//    for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
//    for ( dir = 1 ; dir <= 3 ; dir++ )
//    {
      if ( i == 1 )
        // Predict the sub-step velocity - first order
        // VelCol = VelCol1 + VelColCor1
        for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
        for ( dir = 1 ; dir <= 3 ; dir++ )
          VelCol[ iPar ][ dir ] = VelCol1[ iPar ][ dir ] + VelColCor1[ iPar ][ dir ];
      else
        // Predict the sub-step velocity - second order
        // VelCol = VelCol2 + 2 * VelColCor1
        for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
        for ( dir = 1 ; dir <= 3 ; dir++ )
          VelCol[ iPar ][ dir ] = VelCol2[ iPar ][ dir ] + 2.0 * VelColCor1[ iPar ][ dir ];
      
      // Predict the sub-step position - second order
      // PosCol = PosCol1 + ( VelCol + VelCol1 ) * dts / 2, except for moving frame
      for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
      for ( dir = 1 ; dir <= 3 ; dir++ )
        if ( fluid->FrameVelDir[dir] == false )
          PosCol[ iPar ][ dir ] = PosCol1[ iPar ][ dir ]
            + 0.5 * dts * ( VelCol[ iPar ][ dir ] + VelCol1[ iPar ][ dir ] );
//    }

    // Find collision correction at sub-step position
////////////////////////////////////////////////////////////////////////////////
//      ApplyForce( fluid, dts, PosCol, WallPlanes, VelColCor2 );
////////////////////////////////////////////////////////////////////////////////
    for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
    for ( dir = 1 ; dir <= 3 ; dir++ )
    {
      // Correct the sub-step velocity
      // VelCol = VelCol1 + ( VelColCor1 + VelColCor2 ) / 2
      VelCol[ iPar ][ dir ] = VelCol1[ iPar ][ dir ]
             + 0.5 * ( VelColCor1[ iPar ][ dir ] + VelColCor2[ iPar ][ dir ] );
      // Correct the sub-step position
      // PosCol = PosCol1 + ( VelCol + VelCol1 ) * dts / 2, except for moving frame
      if ( fluid->FrameVelDir[dir] == false )
        PosCol[ iPar ][ dir ] = PosCol1[ iPar ][ dir ]
        + 0.5 * dts * ( VelCol[ iPar ][ dir ] + VelCol1[ iPar ][ dir ] );
  
      // Backup sub-step velocities and position for next step
      VelCol2[ iPar ][ dir ] = VelCol1[ iPar ][ dir ];
      VelCol1[ iPar ][ dir ] = VelCol[ iPar ][ dir ];
      PosCol1[ iPar ][ dir ] = PosCol[ iPar ][ dir ];
    }
  } // end for sub-stepping
  
  // Final predicted position and velocity are in PosCol and VelCol
  // Pos = PosCol
  // Vel = VelCol
  for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    copy3( PosCol[ iPar ], particle[ iPar ].Pos );
    copy3( VelCol[ iPar ], particle[ iPar ].Vel );
  }
}

void
ParColForceCor(
const fluid_t*    fluid,
      particle_t  particle[] )
{
  error("check deeply before any use");
  // check frame change

  int dir, iPar;

  // initialize position
  for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    // update the particles rotation matrix
    RotationVectorToMatrix(particle[ iPar ].Ang, particle[ iPar ].RotMat );

    copy3( particle[ iPar ].Pos1, PosCol[ iPar ] ); // check Pos1
  }

  // Correct position of the particle in sub-steps, except for moving frame
  for ( int i = 1 ; i <= _SubStepNb ; i++ )
  {
    for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
    for ( dir = 1 ; dir <= 3 ; dir++ )
      if ( fluid->FrameVelDir[dir] == false )
        PosCol1[ iPar ][ dir ] = PosCol[ iPar ][ dir ]
          + .5 * dts * ( particle[ iPar ].Vel1[ dir ] + particle[ iPar ].Vel [ dir ] );

////////////////////////////////////////////////////////////////////////////////
//      ApplyForce( fluid, dts, PosCol , WallPlanes, VelColCor1 );
//      ApplyForce( fluid, dts, PosCol1, WallPlanes, VelColCor2 );
////////////////////////////////////////////////////////////////////////////////

    for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
    for ( dir = 1 ; dir <= 3 ; dir++ )
      if ( fluid->FrameVelDir[dir] == false )
        PosCol[ iPar ][ dir ] = PosCol1[ iPar ][ dir ]
          + .5 * dts
          * .5 * ( VelColCor1[ iPar ][ dir ] + VelColCor2[ iPar ][ dir ] );
  }
  
  // Final predicted position and velocity are in PosCol and VelCol
  for ( iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
    copy3( PosCol[ iPar ], particle[iPar].Pos );
}
