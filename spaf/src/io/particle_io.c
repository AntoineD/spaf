/** \file
 Read and write particle data from and to file.
 
 Here checkpoint can be read and written, results data can be stored in a csv file as well as elements inside particles in legacy vtk file for visualization.
 See the functions descriptions for details.
 */

#include "includes.h"
#include "strings.h"
#include "logging.h"
#include "output.h"
#include "linalg.h"
#include "parse.h"
#include "file_system.h"
#include "particle_io.h"
#include "memory.h"
#include "data_io.h"

/// delimiter used between each fields for the results file
#define DELIMITER " "
/// short format for writing data for the results file
#define DATA_FMT "% .3e"
/// long double precision format for writing data for the results file
#define DATA_FMT_FULL "% .18e"

//=============================================================================
/** Write particles results data to file.
 
 The results are written in a file named particle.x.csv where 'x' is the index of the particle, this file is stored in the output directory. It contains a header with the name of the fields written, and a line of data for each iteration written. The fields are separeted by a blank space. If a simulation is restarted from an iteration that is lower than the last iteration stored in the result file, the results will be written for the first iteration that is higher than the last stored. In other words, if 10 iterations are stored and the simulation is restarted from iteration 5, then no result will be written until iteration 11. By default, data are written with 3 significant digits, this can be overidden to be 18 significant digits.
 */
//=============================================================================
void
WriteParticleResults(
const bool        WriteAsDouble,  ///< write as double or in short format ?
const double      time,           ///< time
const int         TimeStep,       ///< time step
const particle_t  particle[] )    ///< particles structure
{
  if ( ! particle ) return;
  
  // to check the first time this routine is called
  static bool FirstTime = true;
  // last time step recorded when restarting computation
  static int LastTimeStepRecorded = 0;
  
  FILE* FileId = NULL;

  debug( "\nSaving particles data\n");
  
  assert_error( particle[1].ParNb <= 1,
    "Modify previous time step checking for more than one particle");

  for ( int iPar = 1 ; iPar <= particle[1].ParNb ; iPar++ )
  {
    // name of the file
    char FileName[30] = "";
    
    // create file name
    sprintf(FileName,"particle." INT_FMT ".csv",iPar);

    if ( TimeStep == 0 )
    {
      // first time step, create file from scratch
      FileId = fopen( FileName, "w" );
      assert_error(FileId, "Failed to open %s""\n", FileName);

      // Write header with data names
      fprintf( FileId,
               "TimeStep" DELIMITER
               "Time" DELIMITER
               "PositionX" DELIMITER
               "PositionY" DELIMITER
               "PositionZ" DELIMITER
               "AngleX" DELIMITER
               "AngleY" DELIMITER
               "AngleZ" DELIMITER
               "VelocityX" DELIMITER
               "VelocityY" DELIMITER
               "VelocityZ" DELIMITER
               "AngularVelocityX" DELIMITER
               "AngularVelocityY" DELIMITER
               "AngularVelocityZ" "\n" );
    }
    else
    {
      // in case we are restarting from a checkpoint and some values have been written
      // after the checkpoint was created, we don't want to write twice the same 
      // values, thus, the first time this file is accessed, we have to find
      // at which time step the last values were written.
      
      // check the file already exists
      if ( ( FirstTime == true ) && ( ( FileId = fopen( FileName, "r" ) ) != NULL ) )
      {
        // removed token checking just for the following parsing
        TokenSetCheck( false );

        // read the last time step written on the file
        while( TokenizeFileLine( FileId ) )
          LastTimeStepRecorded = Token2int( 1, '[', 0, INT_MAX,'[' );

        // restore token checking
        TokenSetCheck( true );
        
        FirstTime = false;
      }
    
      // check the current time step is greater then the last written, if any
      if ( TimeStep > LastTimeStepRecorded )
        // append to file
        FileId = fopen( FileName, "a" );
      else
        return;
    }
      
    fprintf( FileId, "%d" DELIMITER, TimeStep );
    
    fprintf( FileId, DATA_FMT DELIMITER, time );
    
    // building string format
    // the following is for one datum
    char* format1 = NULL;
    
    if ( WriteAsDouble == true )
      format1 = StringConcat(DATA_FMT_FULL, DELIMITER);
    else
      format1 = StringConcat(DATA_FMT, DELIMITER);

    // this is for 3 data
    char* format3 = StringDuplicate(format1);
    StringAppend( &format3, format1 );
    StringAppend( &format3, format1 );
    free(format1);
    
    // just an alias to current particle
    const particle_t *p = &particle[iPar];
    
    fprintf( FileId, format3, p->Pos[1], p->Pos[2], p->Pos[3] );
    fprintf( FileId, format3, p->Ang[1], p->Ang[2], p->Ang[3] );
    fprintf( FileId, format3, p->Vel[1], p->Vel[2], p->Vel[3] );
    fprintf( FileId, format3, p->AngVel[1], p->AngVel[2], p->AngVel[3] );
    fprintf( FileId, "\n" );
    free(format3);
    
    fclose( FileId );
  }
}
//=============================================================================
/** Write the mesh inside a particle into vtk legacy file.
 
 The file is named particle.x.vtk, where 'x' is the index of the particle. All elements inside a particle are stored, original as well as new elements. 
 */
//=============================================================================
void
WriteParMeshVTK(
const int         iPar,     ///< particle index
const mesh_t*     mesh,     ///< mesh structure
const Parmesh_t*  ParMesh ) ///< particle mesh structure
{
  // name of the file
  char FileName[30] = "";
  
  // create file name
  sprintf(FileName,"particle." INT_FMT ".vtk",iPar);
  
  FILE* FileId = fopen( FileName, "w" );
  assert_error(FileId, "Failed to open %s""\n", FileName);

  int NodeNb = (int) ParMesh->points[ 0 ][ 0 ] + mesh->NbOfNodes,
  ElemNb = ParMesh->ElemFull[ 0 ];
  
  // Write header
  fprintf( FileId, "# vtk DataFile Version 2.0" "\n" );
  fprintf( FileId, "particle " INT_FMT " domain mesh" "\n", iPar );
  fprintf( FileId, "ASCII" "\n" );
  fprintf( FileId, "DATASET UNSTRUCTURED_GRID" "\n" );
  
  // Write geometric table
  fprintf( FileId, "POINTS " INT_FMT " float" "\n", NodeNb );
  
  // write all nodes because new elements have some of the old nodes
  for ( int i = 1 ; i <= mesh->NbOfNodes ; i++ )
    fprintf( FileId, "%g %g %g" "\n",
            mesh->points[ i ][ 1 ],
            mesh->points[ i ][ 2 ],
            mesh->points[ i ][ 3 ] );
  
  for ( int i = 1 ; i <= (int) ParMesh->points[ 0 ][ 0 ] ; i++ )
    fprintf( FileId, "%g %g %g" "\n",
            ParMesh->points[ 1 ][ i ],
            ParMesh->points[ 2 ][ i ],
            ParMesh->points[ 3 ][ i ] );
  
  // Write connectivity table for new full elements only
  fprintf( FileId, "CELLS " INT_FMT " " INT_FMT "\n", ElemNb, 11*ElemNb );
  
  for ( int i = 1 ; i <= ElemNb ; i++ )
  {
    int theElem = ParMesh->ElemFull[ i ];
    if ( theElem <= mesh->NbOfElements )
      fprintf( FileId, "10 %d %d %d %d %d %d %d %d %d %d" "\n",
              mesh->ConTab[ theElem ][ 1 ]-1,
              mesh->ConTab[ theElem ][ 2 ]-1,
              mesh->ConTab[ theElem ][ 3 ]-1,
              mesh->ConTab[ theElem ][ 4 ]-1,
              mesh->ConTab[ theElem ][ 5 ]-1,
              mesh->ConTab[ theElem ][ 6 ]-1,
              mesh->ConTab[ theElem ][ 7 ]-1,
              mesh->ConTab[ theElem ][ 8 ]-1,
              mesh->ConTab[ theElem ][ 9 ]-1,
              mesh->ConTab[ theElem ][ 10 ]-1 );
    else
    {
      int iElem = theElem - mesh->NbOfElements;
      fprintf( FileId, "10 %d %d %d %d %d %d %d %d %d %d" "\n",
              ParMesh->ConTab[ 1 ][ iElem ]-1,
              ParMesh->ConTab[ 2 ][ iElem ]-1,
              ParMesh->ConTab[ 3 ][ iElem ]-1,
              ParMesh->ConTab[ 4 ][ iElem ]-1,
              ParMesh->ConTab[ 5 ][ iElem ]-1,
              ParMesh->ConTab[ 6 ][ iElem ]-1,
              ParMesh->ConTab[ 7 ][ iElem ]-1,
              ParMesh->ConTab[ 8 ][ iElem ]-1,
              ParMesh->ConTab[ 9 ][ iElem ]-1,
              ParMesh->ConTab[ 10 ][ iElem ]-1 );
    }
  }
  
  // Write cell type
  fprintf( FileId, "CELL_TYPES " INT_FMT "\n", ElemNb );
  for ( int i = 1 ; i <= ElemNb ; i++ )
    fprintf( FileId, "24" "\n" );
  
  fclose( FileId );
}
//=============================================================================
/** Write particle checkpoint.
 
 For each particles, a directory named particle.x is created, where 'x' is the index of a particle. A data file for each of these vectors is written: Pos1, Pos2, Vel1, Vel2, Ang1, Ang2, AngVel1, AngVel2, Iw1, Iw2. 
 */
//=============================================================================
void
WriteParticleCheckpoint(
const bool        WriteInBinary,  ///< write data in binary ?
const particle_t* particles )     ///< particles structure
{
  // just return if no particle
  if ( ! particles ) return;
  
  for ( int i = 1 ; i <= particles[1].ParNb ; i++ )
  {
    const particle_t* p = &particles[i];
    
    // create directory name as particle.'particle index'
    char suffix[10];
    sprintf(suffix,INT_FMT,i);   
    char* DirName = StringDuplicate("particle.");
    StringAppend( &DirName, suffix );
    
    // create and go to directory DirName
    assert_MakeAndGoToDirectory( DirName );
    
    free(DirName);
    
    WriteArray_double("Pos1", 3, &p->Pos[1], WriteInBinary, true);
    WriteArray_double("Pos2", 3, &p->Pos1[1], WriteInBinary, true);
    WriteArray_double("Vel1", 3, &p->Vel[1], WriteInBinary, true);
    WriteArray_double("Vel2", 3, &p->Vel1[1], WriteInBinary, true);
    WriteArray_double("Ang1", 3, &p->Ang[1], WriteInBinary, true);
    WriteArray_double("Ang2", 3, &p->Ang1[1], WriteInBinary, true);
    WriteArray_double("AngVel1", 3, &p->AngVel[1], WriteInBinary, true);
    WriteArray_double("AngVel2", 3, &p->AngVel1[1], WriteInBinary, true);
    WriteArray_double("Iw1" , 3, &p->Iw[1], WriteInBinary, true);
    WriteArray_double("Iw2" , 3, &p->Iw1[1], WriteInBinary, true);
    
    // go back to initial directory
    assert_GoToDirectoryUp();
  }
}
//=============================================================================
/** Read particle data for checkpoint.
 
 For each particle, read the data stored by WriteParticleCheckpoint. Then for consistency, copy the variables at the previous time step into the current time step ones.
 */
//=============================================================================
void
ReadParticleCheckpoint(
particle_t* particles ) ///< particles structure
{
  // just return if no particle
  if ( ! particles ) return;

  for ( int i = 1 ; i <= particles[1].ParNb ; i++ )
  {
    particle_t* p = &particles[i];
    
    // create directory name as particle.'particle index'
    char sufix[10];
    sprintf(sufix,INT_FMT,i);   
    char* DirName = StringDuplicate("particle.");    
    StringAppend( &DirName, sufix );
    
    // go to directory DirName
    assert_GoBackToDirectory( DirName );
    
    free(DirName);
    
    // read particle variables
    ReadAndCopyArray_double( "Pos1"   , 3, &p->Pos1[1] );
    ReadAndCopyArray_double( "Pos2"   , 3, &p->Pos2[1] );
    ReadAndCopyArray_double( "Vel1"   , 3, &p->Vel1[1] );
    ReadAndCopyArray_double( "Vel2"   , 3, &p->Vel2[1] );
    ReadAndCopyArray_double( "Ang1"   , 3, &p->Ang1[1] );
    ReadAndCopyArray_double( "Ang2"   , 3, &p->Ang2[1] );
    ReadAndCopyArray_double( "AngVel1", 3, &p->AngVel1[1] );
    ReadAndCopyArray_double( "AngVel2", 3, &p->AngVel2[1] );
    ReadAndCopyArray_double( "Iw1"    , 3, &p->Iw1[1] );
    ReadAndCopyArray_double( "Iw2"    , 3, &p->Iw2[1] );
    
    // if the frame is moving, the position in the direction of the movement
    // is not updated and keep its initial value. As we read only the previous
    // fields, this value is not prescribed for the current position, hence we have
    // to copy it back from previous positions.
    // We also fill current Vel and AngVel with previous values, this is safer, and
    // is anyway the actual state in a non-restarted simulation.
    // Note : Vel is to be used to compute the frame velocity and needs to be
    // adequately prescibed.
    copy3( p->Pos1, p->Pos );
    copy3( p->Ang1, p->Ang );
    copy3( p->Vel1, p->Vel );
    copy3( p->AngVel1, p->AngVel );    
    copy3( p->Iw1, p->Iw );    
    
    // go back to timestep directory
    assert_GoToDirectoryUp();
  }
}
