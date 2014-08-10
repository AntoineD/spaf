/** \file
  Manages checkpointing.
 
 A checkpoint contains all necessary data for continuing a simulation that has been stopped.
 All checkpoints are written in a directory called 'checkpoint' in the simulation output directory. A simulation will keep at most the last 2 checkpoints.
 A checkpoint is a directory named after the current timestep. It contains the full fluid field at the previous time step, the fluid field at the previous previous time step (the velocity field only, the pressure is not required for the pressure correction scheme).
 It also contains a dump of the particles variables at previous time steps required for continuing a simulation, see particle_io.c for more details.
 When a simulation starts, if the 'checkpoint' directory is found then the checkpoint with the highest index is read.
 */

#include "checkpoint.h"
#include "logging.h"
#include "file_system.h"
#include "parse.h"
#include "strings.h"
#include "particle_io.h"
#include "fluid_io.h"

//=============================================================================
/** Write a checkpoint.
 
 Writes fluid and particles data in a checkpoint named after the current timestep.
 The data can be written in ascii or binary.
 */
//=============================================================================
void
WriteCheckPoint(
const bool        WriteInBinary,  ///< write in binary or ascii
const int         TimeStep,       ///< timestep
const FluidVar_t* FluidVar,       ///< fluid variables
const particle_t* particles )     ///< particles variables
{  
  // store the previous and current checkpoint directory path
  static char previous[STRING_LEN] = "",
  current[STRING_LEN] = "";
  
  // go to checkpoint directory
  assert_MakeAndGoToDirectory( "checkpoint" );
  
  // newer checkpoint directory
  char newer[STRING_LEN] = "";
  sprintf(newer,INT_FMT,TimeStep);
  
  // rename previous checkpoint directory if it exist
  if ( strlen(previous) != 0 ) rename(previous, newer);
  
  // create directory name based on iteration number
  char DirName[STRING_LEN] = "";
  
  sprintf(DirName, INT_FMT, TimeStep);
  
  // create and/or go to directory DirName
  assert_MakeAndGoToDirectory( DirName );
  
  // fields at previous time step
  WriteField(WriteInBinary, true, "1", FluidVar->NbOfPreNodes, FluidVar->Pre,
             FluidVar->NbOfVelNodes, FluidVar->Vel);
  
  // fields at previous previous time step, no pressure
  WriteField(WriteInBinary, true, "2", 0, NULL, FluidVar->NbOfVelNodes,
             FluidVar->VelOld1);
  
  // write particle data
  WriteParticleCheckpoint(true, particles);
  
  assert_GoToDirectoryUp();
  
  // set current checkpoint as previous and newer one as current
  strcpy(previous, current);
  strcpy(current, newer);
  
  // go back to initial directory
  assert_GoToDirectoryUp();
}
//=============================================================================
/** Read a checkpoint.
 
 Check if a checkpoint is present, if not, do nothing but return false.
 Otherwise, look for the latest checkpoint, then read it into the fluid and particles variables, and return true. 
 */
//=============================================================================
bool
ReadCheckPoint(
int*        TimeStep,   ///< timestep
FluidVar_t* FluidVar,   ///< fluid variables
particle_t* particles ) ///< particles variables
{ 
  // move to output directory and save current working directory
  char* WorkingDirectory = GetWorkingDirectory();
  char* CheckpointDirectory = GoToDirectory("checkpoint");
  
  // is there a checkpoint directory ?
  if ( ! CheckpointDirectory )
  {
    info("\n""\n""No checkpoint found.""\n""\n");
    
    // move back to original directory
    assert_GoBackToDirectory(WorkingDirectory);
    free(WorkingDirectory);
    free(CheckpointDirectory);
    
    return false;
  }
  
  free(CheckpointDirectory);
  
  // list entries in checkpoint directory
  char *FileList = ListFilesInDirectory(".");
  
  TokenSetCheck(false);
  int NbOfEntries = TokenizeString(FileList);
  
  char last[STRING_LEN] = "";
  
  // loop over entries and find last one, i.e. higher number
  // filter numeric entries as well
  for ( int entry = 1 ; entry <= NbOfEntries ; entry++ )
  {
    char *name = Token(entry);
    
    if ( ! StringIsNumber(name) ) continue;
    
    if ( atoi(name) > atoi(last) ) strcpy(last,name);
  }
  
  free(FileList);
  
  // check we have a valid entry
  assert_error(strspn(last, "0123456789") != '\0',
               "No entry in checkpoint directory, remove it or check entries");
  
  // go to that last checkpoint
  assert_GoBackToDirectory( last );
  
  info("Found checkpoint at timestep %s""\n", last);
  
  // get timestep from directory name
  *TimeStep = atoi(last);
  
  // read previous timesteps fluid fields
  ReadField("1", FluidVar->NbOfPreNodes, FluidVar->Pre, FluidVar->NbOfVelNodes,
            FluidVar->VelOld1 );
  
  ReadField("2", 0, NULL, FluidVar->NbOfVelNodes, FluidVar->VelOld2 );
  
  // read previous timesteps particle data
  ReadParticleCheckpoint(particles);
  
  // move back to original directory
  assert_GoBackToDirectory(WorkingDirectory);
  free(WorkingDirectory);
  
  return true;
}
