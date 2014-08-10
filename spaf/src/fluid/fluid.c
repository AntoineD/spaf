/** \file
 Manages fluid variables.
 
 Here fluid physical parameters are read from a .flu file, fluid fields variables are created and allocated as well as shifted in time.
*/

#include "parse.h"
#include "output.h"
#include "logging.h"
#include "memory.h"
#include "linalg.h"
#include "file_system.h"
#include "fluid.h"

//==============================================================================
/** Read fluid physical parameters.
 
 Create, read from a .flu file and return a fluid structure.
 */
//==============================================================================
fluid_t*
ReadFluidParameters(
const char* FileName )  ///< name of the file that contains parameters
{
  PrintTitle("Reading fluid parameters");
  info("from %s\n", FileName );
  
  FILE* FileId = fopen(FileName, "r");
  assert_error(FileId, "Failed to open %s""\n", FileName);
  
  // allocation
  fluid_t* fluid = (fluid_t*) calloc( 1, sizeof(fluid_t) );
  assert_error(fluid, "Memory allocation error" );
  
  // parse file
  while ( TokenizeFileLine( FileId  ) )
  {        
    // Reynolds number
    if ( TokenNameIs( "Reynolds" ) )
      fluid->Re = Token2double( 3 , ']', 0., 1.e4,'[' );
    
    // Froude number
    else if ( TokenNameIs( "Froude" ) )
      fluid->Fr = Token2double( 3 , ']', 0., 1.e4,'[' );
    
    // gravity direction
    else if ( TokenNameIs( "gravity_direction" ) )
    {
      //fluid->GravityDir = Token2int( 3 , '[', -3, 3,']' );
      
      for ( int dir = 1 ; dir <= 3 ; dir++ )
        fluid->Gravity[dir] = Token2double( 2+dir , '[', -DBL_MAX, DBL_MAX,']' );
      
      // make it unitary
      double norm = nrm23(fluid->Gravity);
      assert_error( norm != 0, "Gravity vector is null" );
      
      for ( int dir = 1 ; dir <= 3 ; dir++ )
        fluid->Gravity[dir] /= norm;
    }
    
    // moving frame direction
    else if ( TokenNameIs( "moving_frame_direction" ) )
    {
      int dir = Token2int( 3 , '[', 1, 3,']' );
      fluid->FrameVelDir[dir] = true;
    }
    
    else
      error("bad name token : %s", Token(1) );
  }
  
  // print some info
  info("\nFluid parameters :\n"); 
  info("\tReynolds number   = %g\n", fluid->Re ); 
  info("\tFroude number     = %g\n", fluid->Fr ); 
  info("\tGravity direction =" ); 
  
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    info(" %g", fluid->Gravity[dir] ); 
  
  info("\n"); 
  
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    if ( fluid->FrameVelDir[dir] )
      info("\tMoving frame direction = %d\n", dir );
  
  return fluid;
}
//==============================================================================
/** Allocate a fluid variables structure instance.
 Do not allocate or fill members though.
 */
//==============================================================================
FluidVar_t*
FluidVarAlloc( void )
{
  // allocate structure
  FluidVar_t* FluidVar = (FluidVar_t*) calloc( 1 , sizeof(FluidVar_t) );
  
  assert_error( FluidVar, "Memory allocation error" );

  // initialize members
  FluidVar->Pre = NULL;

  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    FluidVar->VelRHS[dir]  = NULL;
    FluidVar->Vel[dir]     = NULL;
    FluidVar->VelOld1[dir] = NULL;
    FluidVar->VelOld2[dir] = NULL;
    FluidVar->Acc[dir]     = NULL;    
  }
  
  return FluidVar;
}
//==============================================================================
/// Create and initialize the fluid variables structure.
//==============================================================================
FluidVar_t*
FluidVarCreate(
const int NbOfPreNodes,   ///< number of pressure nodes
const int NbOfVelNodes )  ///< number of velocity nodes
{
  // allocate structure
  FluidVar_t* FluidVar = FluidVarAlloc();
  
  // allocate members
  AllocVdouble( NbOfPreNodes, FluidVar->Pre );

  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    AllocVdouble(NbOfVelNodes, FluidVar->VelRHS[dir]  );
    AllocVdouble(NbOfVelNodes, FluidVar->Vel[dir]     );
    AllocVdouble(NbOfVelNodes, FluidVar->VelOld1[dir] );
    AllocVdouble(NbOfVelNodes, FluidVar->VelOld2[dir] );
    AllocVdouble(NbOfVelNodes, FluidVar->Acc[dir]     );
  }
  
  FluidVar->NbOfPreNodes = NbOfPreNodes;
  FluidVar->NbOfVelNodes = NbOfVelNodes;
  
  return FluidVar;
}
//==============================================================================
/// Prepare next time step by copying current variables into those at previous time step.
//==============================================================================
void
FluidShiftVar(
FluidVar_t* FluidVar )  ///< fluid variables
{
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  {
    copy(FluidVar->NbOfVelNodes+1, FluidVar->VelOld1[dir], FluidVar->VelOld2[dir] );
    copy(FluidVar->NbOfVelNodes+1, FluidVar->Vel[dir], FluidVar->VelOld1[dir]);
  }
}
