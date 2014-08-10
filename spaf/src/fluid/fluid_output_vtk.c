/** \file
  Write fluid flow to a vtk legacy file.
 */

#include "logging.h"
#include "file_system.h"

#include "fluid_output_vtk.h"


/// Write velocity and pressure to a legacy vtk file
void
FluidWriteVTK(
const char*   FileName,       ///< file name
const bool    FormatIsBinary, ///< file format : ascii or FormatIsBinary
const bool    WritePre,       ///< save or not pressure field
const mesh_t* mesh,           ///< mesh structure
const double**  Vel,            ///< velocity field
const double*   Pre )           ///< pressure field
{  
  // check file name
  assert( ( FileName != NULL ) && ( strlen(FileName) > 0 ) );
  // check pointers
  assert( mesh != NULL );
  assert( Vel != NULL );
  assert( Pre != NULL );
  
  assert_error( FormatIsBinary == false,
           " Binary format does not work, probably needs big endian");

//=============================================================================

  // full name of the file
  char FullFileName[30] = "";
  
  // create full file name
  sprintf(FullFileName,"%s" ".vtk",FileName);
  
  // open file
  FILE* FileId = fopen( FullFileName, "w" );
  assert_error(FileId, "Failed to open %s", FullFileName);

//=============================================================================

  // Write header
  fprintf( FileId, "# vtk DataFile Version 3.0" "\n" );
  fprintf( FileId, "spaf fluid" "\n" );
  if ( FormatIsBinary ) fprintf( FileId, "BINARY" "\n" );
  else          fprintf( FileId, "ASCII" "\n" );
  fprintf( FileId, "DATASET UNSTRUCTURED_GRID" "\n" );
    
//=============================================================================

  // Write geometric table
  fprintf( FileId, "POINTS " INT_FMT " float" "\n", mesh->NbOfNodes);
  
  if ( FormatIsBinary )
    for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
    for ( int dir = 1 ; dir <= 3 ; dir++ )
      fwrite( (float*) &(mesh->points[NodeId][dir]), sizeof(float), 1, FileId );
  else
    for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
      fprintf( FileId, "%g %g %g" "\n",
               (float) mesh->points[ NodeId ][ 1 ],
               (float) mesh->points[ NodeId ][ 2 ],
               (float) mesh->points[ NodeId ][ 3 ] );
  
//=============================================================================

  // Write connectivity table with 0 offset numbered nodes (hence the "-1")
  fprintf( FileId, "\nCELLS " INT_FMT " " INT_FMT "\n",
           mesh->NbOfElements, 11*mesh->NbOfElements );
  
  if ( FormatIsBinary )
    for ( int ElemId = 1 ; ElemId <= mesh->NbOfElements ; ElemId++ )
    {
      // array containing the element connectivity in proper format for binary legacy vtk
      int ElemCon[11] = {10,0,0,0,0,0,0,0,0,0,0};
      for ( int node = 1 ; node <= 10 ; node++ )
        ElemCon[node] = (int) mesh->ConTab[ElemId][node]-1;
        
      // write element connectivity
      fwrite( ElemCon , sizeof(int) , 11 , FileId );
    }
  else
    for ( int ElemId = 1 ; ElemId <= mesh->NbOfElements ; ElemId++ )
    {
      fprintf( FileId, "10 " );

      for ( int LocalNode = 1 ; LocalNode <= 10 ; LocalNode++ )
        fprintf( FileId, "%d ", (int) mesh->ConTab[ ElemId ][ LocalNode ]-1 );

      fprintf( FileId, "\n" );
    }
  
  // Write cell types
  fprintf( FileId, "\nCELL_TYPES " INT_FMT "\n", mesh->NbOfElements );
  
  if ( FormatIsBinary )
  {
    int CellType = 24;
    
    for ( int ElemId = 1 ; ElemId <= mesh->NbOfElements ; ElemId++ )
      fwrite( &CellType , sizeof(int) , 1 , FileId );
  }
  else
    for ( int ElemId = 1 ; ElemId <= mesh->NbOfElements ; ElemId++ )
      fprintf( FileId, "24" "\n" );
  
//=============================================================================

  // Write data
  fprintf( FileId, "\nPOINT_DATA %d" "\n",mesh->NbOfNodes);

  // Write velocity field
  fprintf( FileId, "VECTORS velocity float" "\n");
  
  if ( FormatIsBinary )
    error("Binary not implemented");
  else
    for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
      fprintf( FileId, "%g %g %g" "\n",
               (float) Vel[ 1 ][ NodeId ],
               (float) Vel[ 2 ][ NodeId ],
               (float) Vel[ 3 ][ NodeId ] );
  
//=============================================================================
  
  // Write pressure field if required
  if ( WritePre == true )
  {
    // vector of pressure value over all nodes
    double *PreInterp = (double*) calloc(mesh->NbOfNodes+1,sizeof(double));
    
    for ( int iElem = 1; iElem <= mesh->NbOfElements; iElem++ )
    for ( int iintocal1 = 1; iintocal1 <= PRE_NODES_PER_EL; iintocal1++ )
    {
      // For pressure nodes, just get pressure from press vector
      int iNodeGlobal1 = mesh->ConTab[ iElem ][ iintocal1 ];
      PreInterp[ iNodeGlobal1 ] = Pre[ mesh->VelToPreNodeMap[ iNodeGlobal1 ] ];
      
      // For edges nodes, interpolate pressure from corner nodes
      for ( int iintocal2 = iintocal1 + 1; iintocal2 <= PRE_NODES_PER_EL; iintocal2++ )
      {
        int iNodeGlobal2    = mesh->ConTab[ iElem ][ iintocal2 ];
        int iNodeGlobalEdge = mesh->ConTab[ iElem ][ mesh->ConTabLocal[ iintocal1 ][ iintocal2 ] ];
        
        PreInterp[ iNodeGlobalEdge ] =
          .5 * ( PreInterp[ iNodeGlobal1 ] +
                 Pre[ mesh->VelToPreNodeMap[ iNodeGlobal2 ] ] );
      }
    }

    fprintf( FileId, "\nSCALARS pressure float 1" "\n");
    fprintf( FileId, "LOOKUP_TABLE default" "\n");

    // write pressure field  
    if ( FormatIsBinary )
      for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
        fwrite( (float*) &(PreInterp[NodeId]) , sizeof(float), 1 , FileId );
    else
      for ( int NodeId = 1 ; NodeId <= mesh->NbOfNodes ; NodeId++ )
        fprintf( FileId, "%g" "\n", (float) PreInterp[ NodeId ] );
  
    free(PreInterp);
  }

  fclose( FileId );
}
