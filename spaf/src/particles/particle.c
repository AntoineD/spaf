/** \file
 Create particles and shift their kinetic data.
 */

#include "includes.h"
#include "memory.h"
#include "linalg.h"
#include "strings.h"
#include "output.h"
#include "logging.h"
#include "parse.h"
#include "collision.h"
#include "file_system.h"
#include "shape_functions.h"
#include "particle.h"

// maximum number of particles
#define PAR_NUMBER_MAX 1

//=============================================================================
/// Allocate and initialize an array of particle structure.
//=============================================================================
static particle_t*
ParCreate(
const int n ) ///< number of particle
{
  particle_t* particle = ( particle_t* ) calloc( n + 1 , sizeof( particle_t ) );

  assert_error( particle != NULL, "Memory allocation error" );

  info( "\nParCreate : allocating %u particle structure \n", n );

  // initialization of every members
  for ( int i = 1 ; i <= n ; i++ )
  {
    // just an alias
    particle_t* p = &particle[i];
    
    // particle number
    p->ParNb = n;

    // particle kind
    p->kind = KIND_NONE;
    p->ellipsoid = NULL;

    // particle properties
    p->density = 0.;
    p->volume  = 0.;
    AllocMdouble( 3,3, p->InertiaTensorInv );
    AllocMdouble( 3,3, p->InertiaTensorInvPhysFrame );

    // kinematic part
    scal3(0., p->Vel  );
    scal3(0., p->Vel1 );
    scal3(0., p->Vel2 );
    scal3(0., p->AngVel );
    scal3(0., p->AngVel1);
    scal3(0., p->AngVel2);
    scal3(0., p->Pos0 );
    scal3(0., p->Pos  );
    scal3(0., p->Pos1 );
    scal3(0., p->Pos2 );
    scal3(0., p->Ang  );
    scal3(0., p->Ang1 );
    scal3(0., p->Ang2 );
    scal3(0., p->Acc );
    scal3(0., p->AngAcc );
    AllocMdouble( 3,3, p->RotMat );
    scal3(0., p->Iw );
    scal3(0., p->Iw1 );
    scal3(0., p->Iw2 );
    
    // mesh
    p->mesh = NULL;
  }
  
  return particle;
}
//=============================================================================
/// Read particles parameters.
//=============================================================================
static void
ReadOneParticleParameters(
      FILE*       FileId,       ///< file identifier
const char*       ParSet,       
const int         ParNb,
      particle_t  particle[] )
{
  // first, get the particle id and from now on, work with Par as a shortcut
  int ParId = String2int( ParSet, '[', 1, ParNb,']' );
  particle_t *Par = &particle[ParId];
  
  while ( TokenizeFileLine( FileId  ) )
  {
    // read density
    if (      TokenNameIs( "density" ) )
      Par->density = Token2double( 3, ']', -DBL_MAX, DBL_MAX,'[' );
    
    // read initial position
    else if ( TokenNameIs( "position" ) )
      for ( int i = 1 ; i <= 3 ; i++ )
        Par->Pos1[i] = Token2double( i+2 , ']', -DBL_MAX, DBL_MAX,'[' );
    
    // read initial velocity
    else if ( TokenNameIs( "velocity" ) )
      for ( int i = 1 ; i <= 3 ; i++ )
        Par->Vel1[i] = Token2double( i+2, ']', -DBL_MAX, DBL_MAX,'[' );
    
    // read initial angular position
    else if ( TokenNameIs( "angular_position" ) )
    {
      for ( int i = 1 ; i <= 3 ; i++ )
        Par->Ang1[i] = Token2double( i+2 , ']', -DBL_MAX, DBL_MAX,'[' );
      
      // check if the angle - vector of rotation form is provided
      if ( TokenCount() == 6 )
      {
        // read the angle
        double angle = Token2double( 6 , '[', -DBL_MAX, DBL_MAX, ']' );
        
        if ( angle == 0. )
        {
          scal3(0.,Par->Ang1);
          
          continue;
        }
        
        // convert angle from degree to radian
        angle *= acos(-1.) / 180.;
        
        // norm of this vector
        double norm = 0.;
        
        for ( int i = 1 ; i <= 3 ; i++ )
          norm += pow( Par->Ang1[i], 2 );
        
        norm = sqrt( norm );
        
        // compute the angular position as the angle time the normed rotation vector
        for ( int i = 1 ; i <= 3 ; i++ )
          Par->Ang1[i] *= angle / norm;
      }
    }
    
    // read initial angular velocity
    else if ( TokenNameIs( "angular_velocity" ) )
      for ( int i = 1 ; i <= 3 ; i++ )
        Par->AngVel1[i] = Token2double( i+2 , ']', -DBL_MAX, DBL_MAX,'[' );
    
    // read parameter for an ellipsoid
    else if ( SectionIs( "ellipsoid" ) )
    {
      Par->kind = KIND_ELLIPSOID;
      Par->ellipsoid = EllipsoidRead( FileId );
      Par->volume = EllipsoidVolume(Par->ellipsoid);
      EllipsoidInertiaTensorInv(Par->ellipsoid,Par->density,Par->InertiaTensorInv);
    }
    
    // end of section ?
    else if ( SectionIsEnd() )
      break;
    
    else
      error( "bad name token : %s", Token(1) );
  }
} 
//=============================================================================
///Read particles parameters.
//=============================================================================
particle_t*
ReadParticlesParameters(
const char* FileName ) ///< name of the file that contains parameters
{
  if ( FileName == NULL ) return NULL;
  
  PrintTitle("Reading particles parameters");
  info("from %s\n", FileName );

  FILE* FileId = fopen(FileName, "r");
  assert_error(FileId, "Failed to open %s\n", FileName);
  
  particle_t* particle = NULL;
  int ParNb = 0;
  
  while ( TokenizeFileLine( FileId ) )
  {        
    // first, get the number of particles and then allocate and initialize the
    // particles structure
    if ( TokenNameIs( "number" ) )
    {
      ParNb = Token2int( 3 , ']', 0, PAR_NUMBER_MAX,']' );
      particle = ParCreate( ParNb );
    }
    
    // Collision parameters
    else if ( SectionIs( "collision" ) )
      ReadCollisionParameters( FileId );
    
    // read particles parameters, the second value token contains the desciption
    // of which particle those parameters apply to, i.e. 
    // section = particle 1, applies to particle 1 only
    else if ( SectionIs( "particle" ) )
      ReadOneParticleParameters( FileId, Token(4), ParNb, particle );
    
    else
      error( "bad name token : %s", Token(1) );
  }
  
  
  // print particle structure infos
  for ( int iPar = 1 ; iPar <= ParNb ; iPar++ )
  {
    info( "\nParticle %u properties and initial conditions :\n", iPar );
    
    // set up a shortcut
    const particle_t *Par = &particle[iPar];
    
    // print kind
    switch ( Par->kind )
    {
      case KIND_ELLIPSOID :
        EllipsoidPrint( Par->ellipsoid );
        break;
        
      default :
        error( "%s : bad particle kind", __FUNCTION__ );
        break;
    }
    
    info( "\trelative density = %g\n", Par->density );
    info( "\tvolume           = %g\n", Par->volume );
    info( "\tInertia tensor inverted = |\t%g\t%g\t%g\t|\n",
         Par->InertiaTensorInv[1][1],
         Par->InertiaTensorInv[1][2],
         Par->InertiaTensorInv[1][3] );
    info( "\t                          |\t%g\t%g\t%g\t|\n",
         Par->InertiaTensorInv[2][1],
         Par->InertiaTensorInv[2][2],
         Par->InertiaTensorInv[2][3] );
    info( "\t                          |\t%g\t%g\t%g\t|\n",
         Par->InertiaTensorInv[3][1],
         Par->InertiaTensorInv[3][2],
         Par->InertiaTensorInv[3][3] );
    
    info( "\tposition         = %g\t%g\t%g\n",
         Par->Pos1[1],Par->Pos1[2],Par->Pos1[3] );
    info( "\tvelocity         = %g\t%g\t%g\n",
         Par->Vel1[1],Par->Vel1[2],Par->Vel1[3] );
    info( "\torientation      = %g\t%g\t%g\n",
         Par->Ang1[1],Par->Ang1[2],Par->Ang1[3] );
    info( "\tangular velocity = %g\t%g\t%g\n",
         Par->AngVel1[1],Par->AngVel1[2],Par->AngVel1[3] );
  }

  
  // post reading stuff, initialization
  for ( int iPar = 1 ; iPar <= ParNb ; iPar++ )
  {
    // just an alias
    particle_t *p = &particle[iPar];    
    
    // initialize Inertia tensor times angular velocity
    // we need the rotation matrix
    RotationVectorToMatrix( p->Ang1, p->RotMat );
    
    // then the inertia tensor inverted in physical frame
    ChangeMatrixBase(
      p->InertiaTensorInv, p->RotMat, p->InertiaTensorInvPhysFrame );
    
    // get the inverse
    double **I = NULL;
    AllocMdouble(3,3,I);
    inv3by3(p->InertiaTensorInvPhysFrame, I);

    // then multiply
    gemv3(I, p->AngVel1, p->Iw1);
    
    FreeM(I);
    
    copy3( p->Pos1, p->Pos0 );
    copy3( p->Pos1, p->Pos );
    copy3( p->Vel1, p->Vel );
    copy3( p->AngVel1, p->AngVel );
    copy3( p->Ang1, p->Ang );
    
    // allocate particle mesh
    p->mesh = (Parmesh_t*) calloc( 1, sizeof(Parmesh_t) );
    
    // just a alias
    Parmesh_t* mesh = particle[iPar].mesh;

    // in general, first item is the number of elements in the array
    AllocVint(1-1,mesh->ElemFull);
    AllocVint(1-1,mesh->ElemPartly);
    
    mesh->NodesIn = NULL;
    
    // those arrays are not usual 2d arrays since will need to reallocate
    // with incresing nodes number
    mesh->points = (double**) calloc( 4, sizeof(double*) );
    assert_error( mesh->points != NULL, "Memory allocation error");
    
    for ( int dir = 0 ; dir <= 3 ; dir++ )
    {
      mesh->points[dir] = (double*) calloc(1,sizeof(double));
      assert_error( mesh->points[dir] != NULL, "Memory allocation error");
    }

    mesh->ConTab = (int**) calloc( 11, sizeof(int*) );
    assert_error( mesh->ConTab != NULL, "Memory allocation error");
    
    for ( int i = 0 ; i <= 10 ; i++ )
    {
      mesh->ConTab[ i ] = (int*) calloc( 1, sizeof(int) );
      assert_error( mesh->ConTab[ i ] != NULL, "Memory allocation error");
    }
  }

  return particle;
} 
//=============================================================================
/// Prepare next time step by copying current variables into those at previous time step
//=============================================================================
void
ParShiftVar(
particle_t particle[] )  ///< particle structure
{
  if ( particle == NULL ) return;

  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    copy3( particle[ iPar ].Pos1  , particle[ iPar ].Pos2    );
    copy3( particle[ iPar ].Pos   , particle[ iPar ].Pos1    );
    copy3( particle[ iPar ].Vel1  , particle[ iPar ].Vel2    );
    copy3( particle[ iPar ].Vel   , particle[ iPar ].Vel1    );
    copy3( particle[ iPar ].Ang1  , particle[ iPar ].Ang2    );
    copy3( particle[ iPar ].Ang   , particle[ iPar ].Ang1    );
    copy3( particle[ iPar ].AngVel1, particle[ iPar ].AngVel2 );
    copy3( particle[ iPar ].AngVel, particle[ iPar ].AngVel1 );
    copy3( particle[ iPar ].Iw1   , particle[ iPar ].Iw2     );
    copy3( particle[ iPar ].Iw    , particle[ iPar ].Iw1     );
  }
}
