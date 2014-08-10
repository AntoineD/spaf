/** \file
 Create a planes structure that hold planes geometric definition, those can be used for collision checking.
 */

#include "parse.h"
#include "logging.h"
#include "file_system.h"
#include "planes.h"

// maximum number of planes
#define PLANES_MAX_NUMBER 10

//=============================================================================
/** Read the definitions of planes from disc.
 
 Return NULL if there is no such file.
 */
//=============================================================================
planes_t*
PlanesRead(
const char* FileName )  ///< name of the file with definitions
{
  FILE* FileId = fopen(FileName, "r");
  
  if ( ! FileId ) return NULL;
  
  planes_t *planes = (planes_t*) calloc( 1 , sizeof(planes_t) );

	// read the number of planes
  TokenizeFileLine( FileId );

  assert_error( TokenNameIs( "number" ), "Misssing number of planes" );

  planes->n = Token2int( 3 , '[', 0, PLANES_MAX_NUMBER,']' );
  planes->equation = (equation_t*) calloc( planes->n + 1 , sizeof(equation_t) );

  // read every plane coefficients
  for ( int iPlane = 1 ; iPlane <= planes->n ; iPlane++ )
	{
		TokenizeFileLine( FileId );
		
		assert_error( TokenNameIs( "plane" ), "Missing plane" );

    for ( int iCoef = 0 ; iCoef < 4 ; iCoef++ )
      planes->equation[iPlane].coef[iCoef] =
        Token2double( iCoef+3, ']', -DBL_MAX, DBL_MAX,'[' );
	}
  
  fclose(FileId);
  
  // print some infos
  info("\nPlanes equation :\n"); 
  
  for ( int i = 1 ; i <= planes->n ; i++ )
  	info("\tplane %d: 0 = %g \t+ %g x \t+ %g y \t+ %g z\n", i,
         planes->equation[i].coef[0],
         planes->equation[i].coef[1],
         planes->equation[i].coef[2],
         planes->equation[i].coef[3] );
  
  return planes;
}
//=============================================================================
/** Compute and return the distance from a point to a plane.
 */
//=============================================================================
double
PlaneToPointDistance(
const double      point[4],   ///< position of a point
const equation_t* equation )  ///< plane equation
{
  return
  equation->coef[0] +
  equation->coef[1] * point[ 1 ] +
  equation->coef[2] * point[ 2 ] +
  equation->coef[3] * point[ 3 ];
}
