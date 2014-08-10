/** \file
 Read and write fluid fields to and from disc, and visualization data as well.
 
 A fluid field is a directory which contains a velocity field and/or a pressure field.
 The name of the fluid field is the name of the directory.
 The pressure field is a data file named 'pressure'.
 The velocity field is composed of 3 data files named 'velocity.1', 'velocity.2', 'velocity.3'.
 See data_io.c for a description of a data file.
 
 For visualization data, write the fluid field named after the current timestep, and the current time is also written in a data file.
 */

#include "fluid_io.h"
#include "file_system.h"
#include "data_io.h"
#include "logging.h"

//=============================================================================
/** Write a fluid field data to disc.
 
 If the velocity pointer is not NULL, write the velocity field.
 Same for the pressure field if the pressure.
 Those data could be written in binary or ascii, and as float or double.
 */
//=============================================================================
void
WriteField(
const bool      WriteInBinary,  ///< write in binary or ascii
const bool      WriteAsDouble,  ///< write as double or float
const char*     FieldName,      ///< name of the fluid field
const int       NbOfPreNodes,   ///< number of pressure nodes
const double*   Pre,            ///< pressure field
const int       NbOfVelNodes,   ///< number of velocity nodes
      double**  Vel )           ///< velocity field
{
  // save current directory then make and go to the fluid field directory
  char* WorkingDirectory = GetWorkingDirectory();
  assert_MakeAndGoToDirectory( FieldName );
  
  if ( Vel )
  {
    WriteArray_double("velocity.1", NbOfVelNodes, &Vel[1][1], WriteInBinary,
                      WriteAsDouble);
    
    WriteArray_double("velocity.2", NbOfVelNodes, &Vel[2][1], WriteInBinary,
                      WriteAsDouble);
    
    WriteArray_double("velocity.3", NbOfVelNodes, &Vel[3][1], WriteInBinary,
                      WriteAsDouble);
  }
  
  if ( Pre )
    WriteArray_double("pressure", NbOfPreNodes, &Pre[1], WriteInBinary,
                      WriteAsDouble);
  
  // go back to saved directory
  assert_GoBackToDirectory(WorkingDirectory);
  free(WorkingDirectory);
}
//=============================================================================
/** Read fluid data from disc.
 
 Read the velocity and/or pressure fields from a fluid field.
 They have to be already allocated.
 */
//=============================================================================
void
ReadField(
const char*     FieldName,    ///< name of the fluid field 
const int       NbOfPreNodes, ///< number of pressure nodes
      double*   Pre,          ///< pressure field
const int       NbOfVelNodes, ///< number of velocity nodes
      double**  Vel )         ///< velocity field
{
  // save current directory then make and go to the fluid field directory
  char *WorkingDirectory = GoToDirectory(FieldName);  
  assert_Directory(WorkingDirectory, "field");
  
  if ( Vel )
  {
    ReadAndCopyArray_double("velocity.1", NbOfVelNodes, &Vel[1][1]);
    ReadAndCopyArray_double("velocity.2", NbOfVelNodes, &Vel[2][1]);      
    ReadAndCopyArray_double("velocity.3", NbOfVelNodes, &Vel[3][1]);
  }
  
  if ( Pre )
    ReadAndCopyArray_double("pressure", NbOfPreNodes, &Pre[1]);    
  
  // go back to saved directory
  assert_GoBackToDirectory(WorkingDirectory);
  free(WorkingDirectory);
}
//=============================================================================
/** Write data for visualization.
 
 Write the fluid field named after the timestep, write the current time as well.
 */
//=============================================================================
void
WriteVizData(
const bool        WriteInBinary,  ///< write in binary or ascii
const bool        WriteAsDouble,  ///< write as double or float
const double      time,           ///< time
const int         TimeStep,       ///< timestep
const FluidVar_t* FluidVar)       ///< fluid fields
{
  assert_MakeAndGoToDirectory( "viz" );
  
  // create a directory named after the timestep
  char DirName[STRING_LEN] = "";
  sprintf(DirName, INT_FMT, TimeStep);
  assert_MakeAndGoToDirectory( DirName );
  
  // write fluid field in current directory
  WriteField(WriteAsDouble, WriteInBinary, ".", FluidVar->NbOfPreNodes,
             FluidVar->Pre, FluidVar->NbOfVelNodes, FluidVar->Vel);
  
  // write time too in ascii and float
  WriteArray_double( "time", 1, &time, false, true);
  
  assert_GoToDirectoryUp();
  assert_GoToDirectoryUp();
}
