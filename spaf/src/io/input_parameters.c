/** \file
  Read numerical parameters from a .num file.
 */

#include "logging.h"
#include "output.h"
#include "parse.h"
#include "planes.h"
#include "fluid.h"
#include "output.h"
#include "legacy.h"
#include "convection.h"
#include "mesh_search.h"
#include "pcg.h"
#include "rigid_body.h"
#include "file_system.h"
#include "input_parameters.h"
#include "operators.h"

//=============================================================================
/** Read time parameters.
*/
//=============================================================================
static void
ReadTimeParameters(
FILE*   FileId,       ///< file id
double* dt,           ///< timestep interval
int*    TimeStepNb )  ///< number of time steps
{
  // read time step
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "dt" ), "Missing time parameter : ""dt");
  *dt = Token2double( 3, ']', 0., 10., ']' );
  
  // read number of time steps
  TokenizeFileLine( FileId );
  assert_error( TokenNameIs( "time_step_number" ),
               "Missing time parameter : ""time_step_number");
  *TimeStepNb = Token2int( 3, '[', 1, INT_MAX, ']' );
  
  // check end of section
  TokenizeFileLine( FileId );
  assert_error( SectionIsEnd(), "Time parameters : missing section = end");
  
  // log infos
  info("\nTime parameters :\n");
  info("\tNumber of time steps = " INT_FMT "\n", *TimeStepNb );
  info("\tTime step = %g\n", *dt );
}
//=============================================================================
/** Read numerical parameters.
 */
//=============================================================================
void
ReadNumericalParameters(
const char*   PathToFile,     ///< path to parameter file
      double* dt,             ///< time step interval
      int*    NbOfTimeStep )  ///< number of time step
{
  PrintTitle("Reading numerical parameters");
  info("From %s\n",PathToFile);
  
  FILE* FileId = fopen( PathToFile, "r" );

  assert_error( FileId, "Cannot open numerical parameters""\n");
  
  // read time param
  TokenizeFileLine( FileId );
  assert_error( SectionIs("time"), "time section not found");
  ReadTimeParameters( FileId, dt, NbOfTimeStep );
  
  // read output param
  TokenizeFileLine( FileId );
  assert_error( SectionIs("output"), "output section not found");
  ReadOutputParameters( FileId );
  
  // read conjugate gradient param
  TokenizeFileLine( FileId );
  assert_error( SectionIs("conjugate_gradient"),
               "conjugate_gradient section not found");
  ReadConjugateGradientParameters( FileId );

  // read conjugate gradient param
  TokenizeFileLine( FileId );
  assert_error( SectionIs("preconditioner"),
               "preconditioner section not found");
  ReadPreconditionerParameters( FileId );
  
  // read rigid body param
  TokenizeFileLine( FileId );
  assert_error( SectionIs("rigid_body"), "rigid_body section not found");
  ReadRigidBodyMotionParameters( FileId );
  
  // read bounding box param
  TokenizeFileLine( FileId );
  assert_error( SectionIs("nearest_neighbor_search"),
               "nearest_neighbor_search section not found");
  ReadNNSParameters( FileId );
  
  // read convection param
  TokenizeFileLine( FileId );
  assert_error( SectionIs("convection"), "convection section not found");
  ReadConvectionParameters( FileId );
  
  // read mesh_search param
  TokenizeFileLine( FileId );
  assert_error( SectionIs("mesh_search"), "mesh_search section not found");
  ReadMeshSearchParameters( FileId );
  
  fclose( FileId );
}
