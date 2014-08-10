/** \file
  Read or write info and data of a simulation.
 
 This file is used to read or write informations or data of a simulation, it can:
 - read or write checkpoints
 - write visualisation data
 - write particle kinetic data
 - write statistics
 - write residuals for fluid only simulation
 - print titles
 - print time step header 
*/

#include "includes.h"
#include "parse.h"
#include "data_io.h"
#include "strings.h"
#include "logging.h"
#include "particle_io.h"
#include "statistics.h"
#include "linalg.h"
#include "memory.h"
#include "file_system.h"
#include "output.h"
#include "checkpoint.h"
#include "fluid_io.h"

// file local variables, names say it all, default is to write nothing
static int
  _period_of_checkpoint_save = 0,
  _period_of_particle_save = 0,
  _period_of_stat_save = 0,
  _period_of_fluid_save = 0,
  _period_of_residuals_save = 0;

// flag to specify data format
// this flag is not read from parameters file
static bool
  _write_in_binary = true,
  _write_as_double = false;

// stores the output directory
static char* _path_to_output_dir = NULL;

//=============================================================================
/** Set the path to the output directory.
*/
//=============================================================================
void
SetOutputDirectory(
char* PathToOutput ) ///< path to output directory
{
  _path_to_output_dir = PathToOutput;

  // check this directory exists
  assert_error(CheckFileExistence(_path_to_output_dir),
               "Output directory %s does not exist", _path_to_output_dir);
}
//=============================================================================
/** Go to the output directory.
 
 This function set the current working directory to the output directory, it
 also return the previous working directory. If something goes wrong, it
 returns NULL.
*/
//=============================================================================
char*
GoToOutputDirectory( void )
{
  return GoToDirectory(_path_to_output_dir);
}
//=============================================================================
/** Read parameters for writing ouput.
 
*/
//=============================================================================
void
ReadOutputParameters(
FILE* FileId )  ///< file identifier
{
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "save_period_checkpoint" ),
               "Missing output parameter : ""save_period_checkpoint");
  _period_of_checkpoint_save = Token2int( 3 , '[', 0, INT_MAX,']' );
//------------------------------------------------------------------------------
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "save_period_particles" ),
    "Missing output parameter : ""save_period_particles");
  _period_of_particle_save = Token2int( 3 , '[', 0, INT_MAX,']' );
//------------------------------------------------------------------------------
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "save_period_fluid" ),
    "Missing output parameter : ""save_period_fluid");
  _period_of_fluid_save = Token2int( 3 , '[', 0, INT_MAX,']' );
//------------------------------------------------------------------------------
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "save_period_stat" ),
    "Missing output parameter : ""save_period_stat");
  _period_of_stat_save = Token2int( 3 , '[', 0, INT_MAX,']' );
//------------------------------------------------------------------------------
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "write_as_double" ),
    "Missing output parameter : ""write_as_double");
  if ( StringCompare(Token(3),"double") )
    _write_as_double = true;
  else if ( StringCompare(Token(3),"float") )
    _write_as_double = false;
  else
    error("bad output parameter : write_as_double : %s", Token(3));
//------------------------------------------------------------------------------
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "write_in_binary" ),
    "Missing output parameter : ""write_in_binary");
  if ( StringCompare(Token(3),"binary") )
    _write_in_binary = true;
  else if ( StringCompare(Token(3),"ascii") )
    _write_in_binary = false;
  else
    error("bad output parameter : write_in_binary : %s", Token(3));  
//------------------------------------------------------------------------------
  TokenizeFileLine( FileId );
  assert_error( SectionIsEnd(),
    "Output parameters : missing section = end");

//------------------------------------------------------------------------------
  // print infos
  info("\nOutput parameters :\n");
  
  info("\tPeriod to save checkpoint                = %lu\n",
         _period_of_checkpoint_save );
  
  info("\tPeriod to save particles kinetic history = %lu\n",
         _period_of_particle_save );
  
  info("\tPeriod to save fluid kinetic history     = %lu\n",
         _period_of_fluid_save );
  
  info("\tPeriod to save statistics history        = %lu\n",
         _period_of_stat_save );
  
  if ( _write_in_binary )
    info("\tData files will be written in binary and ");
  else
    info("\tData files will be written in ascii and ");
    
  if ( _write_as_double )
    info("as double.\n");
  else
    info("as float.\n");
}
//=============================================================================
/** Print a title.
 */
//=============================================================================
void
PrintTitle(
const char title[] )  ///< string to be written
{
  info( "\n\n---------------------------------------------------------------------------------\n" );
  info( "\t%s", title );
  info( "\n---------------------------------------------------------------------------------\n" );
}
//=============================================================================
/** Print a timestep title.
 */
//=============================================================================
void
PrintTimeStep(
const double  time,       ///< time
const int     timestep )  ///< time step
{
  char string[200] = "";
  sprintf( string,"TIMESTEP "INT_FMT" , time = "DBL_FMT, timestep, time );
  debug( "\n\n-----------------------------------------------------\n" );
  debug( "\t%s", string );
  debug( "\n-----------------------------------------------------\n" );
}  
//==============================================================================
/** Manage writting to data to disc.

  Depending on the timestep and file local parameters, data for visualization, checkpoint, particle tracking, statistics, fluid only convergence history, are written.
*/
//==============================================================================
void
WriteOutput(
const FluidVar_t* FluidVar,   ///< fluid variables
const particle_t* particle,   ///< particles structure
const double      time,       ///< time
const int         TimeStep )  ///< time step
{
  // visualization data
  if ( _period_of_fluid_save != 0 &&
       TimeStep % _period_of_fluid_save == 0 )
    WriteVizData(_write_in_binary, _write_as_double, time, TimeStep, FluidVar);

  // checkpoint
  if ( TimeStep != 0 &&
       _period_of_checkpoint_save != 0 &&
       TimeStep % _period_of_checkpoint_save == 0 )
    WriteCheckPoint(_write_in_binary, TimeStep, FluidVar, particle);
  
  
  // particles data
  if ( _period_of_particle_save != 0 &&
       TimeStep % _period_of_particle_save == 0 )
    WriteParticleResults(_write_as_double, time, TimeStep, particle );
  
  // write stats
  if ( _period_of_stat_save != 0 &&
       TimeStep % _period_of_stat_save == 0 )
    StatFlush("stat.log");
  
  // for fluid only if there is no particle
  // Compute time max error
  if ( particle == NULL &&
       _period_of_residuals_save != 0 &&
       TimeStep % _period_of_residuals_save == 0 )
  {
    static bool FirstTime = true;
    FILE* FileId = NULL;
    
    if ( FirstTime )
    {
      FileId = fopen("residuals.csv", "w");
      assert_error(FileId, "Failed to open residuals.csv""\n");
      fprintf(FileId, "timestep v1 v2 v3""\n");
      FirstTime = false;
    }
    else
    {
      FileId = fopen("residuals.csv", "a");
      assert_error(FileId, "Failed to open residuals.csv""\n");
    }
      
    debug("\n" "Max error on velocity update" "\n");
    fprintf(FileId,INT_FMT, TimeStep);
    
    // measure fluid time error
    for ( int dir = 1 ; dir <= 3 ; dir++ )
    {
      double error_sup = 0.;
      
      for ( int NodeId = 1 ; NodeId <= FluidVar->NbOfVelNodes ; NodeId++ )
      {
        double error = fabs(FluidVar->Vel[dir][NodeId] -
                            FluidVar->VelOld1[dir][NodeId] );
        error_sup = ( error_sup > error ) ? error_sup : error;
      }
      
      debug("\t" "dir " INT_FMT " = " DBL_FMT,dir,error_sup);
      fprintf(FileId," "DBL_FMT, error_sup);
    }
    
    debug("\n");
    fprintf(FileId,"\n");
    fclose(FileId);
  }
}
