/** \file
 Manage boundary conditions.
 
 Here we read the bc from the disc and create a bc structure, which can then be used for imposing the bc.
 
 For the velocity, Dirichlet or homogeneous Neumann type of bc can be prescribed, but for all the 3 components at the same nodes. In other words, at any one node, the velocity bc must be of the same type for all the 3 components. The nodes at which homogeneous Neumann bc is defined are considered as non-prescribed nodes, thus no special treatment is required in practical terms. This bc type is only used for outflow bc.
 
 For the pressure, only homogeneous Dirichlet bc can be prescribed. So far in the code the only situation that require a pressure bc is the outflow bc, which is treated with zero pressure. This means in practical terms that no special treatment has to be done for the pressure, and that the pressure members of the bc structure are not used in the current code. Nevertheless the current code can read and store pressure Dirichlet bc, but no function has been written to impose them. If there is need for such a function, just adapt the one used for the velocity.
 */

#include "memory.h"
#include "logging.h"
#include "strings.h"
#include "data_io.h"
#include "file_system.h"
#include "output.h"
#include "bc.h"

/// Boundary conditions structure
struct bc_struct
{
  // velocity
  int NbOfVelNodes;       ///< number of prescribed nodes
  double *velocity[3+1];  ///< values of the prescribed nodes, one array per component
  
  // pressure
  int NbOfPreNodes; ///< number of prescribed nodes
  double *pressure; ///< values of the prescribed nodes
};

//==============================================================================
/** Read, initialize and prepare bc.
 
  This function reads the bc from the disc, then fill in the bc structures for the
  velocity and pressure.
 Returns the bc structure.
*/
//==============================================================================
bc_t*
ReadBC(
const char* PathToCase ) ///< path to directory containing bc files
{
  PrintTitle("Reading boundary conditions");

  // get the path to the mesh directory inside the case one
  char *PathToBC = StringConcat(PathToCase, "/boundary_conditions");
  
  info("from %s\n", PathToBC );
  
  // move to bc directory
  char *WorkingDirectory = GoToDirectory(PathToBC);  

  assert_error(WorkingDirectory, "Cannot open bc""\n");

  free(PathToBC);
  
  bc_t* bc = (bc_t*) calloc(1, sizeof(bc_t));

  int number = 0;
  
  ReadArray_double("velocity.1", &number, &bc->velocity[1]);
  bc->NbOfVelNodes = number;
  
  ReadArray_double("velocity.2", &number, &bc->velocity[2]);
  assert_error( number == bc->NbOfVelNodes,
    "Number of prescribed nodes do not match");
  
  ReadArray_double("velocity.3", &number, &bc->velocity[3]);
  assert_error( number == bc->NbOfVelNodes,
               "Number of prescribed nodes do not match");
  
  ReadArray_double("pressure", &number, &bc->pressure);
  bc->NbOfPreNodes = number;

  // print some infos
  info(INT_FMT" prescribed velocity nodes""\n", bc->NbOfVelNodes);
  info(INT_FMT" prescribed pressure nodes""\n", bc->NbOfPreNodes);
  
  // go back to previous workind directory
  assert_GoBackToDirectory(WorkingDirectory);
  free(WorkingDirectory);
  
  return bc;
}
//=============================================================================
/** Prescibe bc for the velocity.
 */
//=============================================================================
void
SetVelBC(
const int       NbOfNodes,  ///< total number of velocity nodes
const bc_t*     bc,         ///< boundary conditions
      double**  Vel )       ///< velocity field
{
  // first prescribed node
  int FirstNode = NbOfNodes - bc->NbOfVelNodes + 1;
  
  for ( int dir = 1 ; dir <= 3 ; dir++ )
  for ( int node = 0 ; node < bc->NbOfVelNodes ; node++ )
    Vel[ dir ][ FirstNode + node ] = bc->velocity[ dir ][node];
}
