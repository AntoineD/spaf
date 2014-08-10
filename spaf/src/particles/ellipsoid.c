/** \file
  Manages ellipsoid particles.
*/

#include "includes.h"
#include "logging.h"
#include "linalg.h"
#include "parse.h"

#include "ellipsoid.h"

// minimum value of the semi axes
#define SEMI_AXES_MIN -1.e3
// maximum value of the semi axes
#define SEMI_AXES_MAX 1.e3

/*--------------------------------------------------------------------
Purpose:
  allocate and initialize an ellipsoid structure.

Input:
  none
  
Output:
  ellipsoid structure allocated and initialized.

By: Antoine
Update: 27/09/06
--------------------------------------------------------------------*/
ellipsoid_t*
EllipsoidCreate( void )
{
  short int i;
  ellipsoid_t* ellipsoid = ( ellipsoid_t* ) malloc( sizeof( ellipsoid_t ) );

  assert_error( ellipsoid != NULL, "Memory allocation error" );

  for ( i = 1 ; i <= 3 ; i++ )
    ellipsoid->SemiAxes[i] = 0.0;

  return ellipsoid;
}
//=============================================================================
// Print a ellipsoid infos in a file or in std output
// Antoine
//=============================================================================
void
EllipsoidPrint(
const ellipsoid_t* ellipsoid )
{
  info( "\tEllipsoid semi-axes = %g\t%g\t%g\n",
    ellipsoid->SemiAxes[1],
    ellipsoid->SemiAxes[2],
    ellipsoid->SemiAxes[3] );
}
/*==============================================================================
PURPOSE
  Parse fluid properties
 
INPUT
  FILE* FileId : file to read from   

OUTPUT :
  fluid_t* fluid : fluid structure filled
  
Antoine, 10/11/06
==============================================================================*/
ellipsoid_t*
EllipsoidRead(
FILE* FileId )
{
  ellipsoid_t* ellipsoid = EllipsoidCreate(); // first, allocate and initialize structure

  while ( TokenizeFileLine( FileId ) )
  {
    if ( TokenNameIs( "semi_axes" ) )
      for ( int i = 1 ; i <= 3 ; i++ )
        ellipsoid->SemiAxes[i] = Token2double( i+2, ']', SEMI_AXES_MIN, SEMI_AXES_MAX,'[' );
    
    else if ( SectionIsEnd() )
      break;
      
    else
      error( "EllipsoidRead : bad name token : %s\n", Token(1) );
  }
  
//  EllipsoidPrint( ellipsoid );

  return ellipsoid;
} 
/*--------------------------------------------------------------------
Purpose:
  Is a point in an ellipsoid ?

Input:
  point vector
  ellipsoid structure
  
Output:
  FALSE if the point is not in the ellipsoid.

By: Antoine
Update: 27/09/06
--------------------------------------------------------------------*/
bool
PointIsInEllipsoid(
const ellipsoid_t*  ellipsoid,
const double          position[ 4 ] )
{
  double dist_sq = 0.0;
  
  // compute (x/a)^2 + (y/b)^2 + (z/c)^2
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    dist_sq += pow( position[ dir ] / ellipsoid->SemiAxes[ dir ], 2 );

  // (x/a)^2 + (y/b)^2 + (z/c)^2 <= 1.0 ?
  return dist_sq <= 1.;
}
/*--------------------------------------------------------------------
PURPOSE:
  To solve for the intersection point between ellipsoid and edge.
  Let P = (x,y,z) be the intersection point, Pin and Pout be the edge
  end points coordinates.
  Equation of ellipsoid in its coordinate system (X,Y,Z):
  (X/e1)^2 + (Y/e2)^2 + (Z/e3)^2 = 1
  ( In order to solve the intersection point, we better use this simple form
  of the ellipsoid equation, thus Pin and Pout need to be expressed in terms
  of (X,Y,Z): Pin -> Min  and  Pout -> Mout. ) already done in this version.
  Equation of edge: vect(MinM) = k * vect(MinMout)
  then xi = ai + bi*t, i = 1,2,3 ; with a = Pin and b = Pout - Pin.
  Take t = 0 for Pin, t = 1 for Pout ==> get ai and bi
  Substitute xi in particle equation and solve for t in [0,1] :
    P * t^2  +  Q * t  +  R  =  0
  From t ==> coordinates of new node (intersection point)

INPUT:
  ellipsoid: the particle infos
  Pin, Pout: points coordinates

OUTPUT:
  t value
  Pin and Pout are modified

By: Veeramani
Update: 02/10/06 Antoine
--------------------------------------------------------------------*/
double
GetEdgeIntersectParamEllipsoid(
const ellipsoid_t*  ellipsoid,
      double          Pin[4],
      double          Pout[4] )
{
  double t,
        P = 0.0, Q = 0.0, R = 0.0, D = 0.0;
  int dir;

  for ( dir = 1 ; dir <= 3 ; dir++ )
    Pout[ dir ] -= Pin[ dir ]; // Pout holds b = Pout - Pin

  for ( dir = 1 ; dir <= 3 ; dir++ )
  {
    Pin[ dir ]  /= ellipsoid->SemiAxes[dir]; // Pin  holds ai / ei
    Pout[ dir ] /= ellipsoid->SemiAxes[dir]; // Pout holds bi / ei
  }
  
  // set polynomial coefficients
  for ( dir = 1 ; dir <= 3 ; dir++ )
  {
    P += pow( Pout[ dir ] , 2 );    // P = sum ( bi / ei )^2
    Q += Pin[ dir ] * Pout[ dir ];  // Q = sum ai * bi / ei^2
    R += pow( Pin[ dir ] , 2 );     // R = sum ( ai / ei )^2
  }

  Q *= 2;   // Q = 2 * sum ai * bi / ei^2
  R -= 1.0; // R = sum ( ai / ei )^2 - 1

  // Now to solve P*t^2 + Q*t + R = 0
  D = Q * Q - 4 * P * R;

  if ( D == 0.0 )
  { // One repeated double root exist
    t = - Q / ( 2 * P );
    if ( 0.0 <= t && t <= 1.0 ) // Proper value of t found
      return t;
  }
  else if ( D > 0.0 )
  { // Two different double roots exist
    D = sqrt( D );
    t = ( D - Q ) / ( 2 * P ); // get the first root
    
    if ( 0.0 <= t && t <= 1.0 ) // Proper value of t found
      return t;
    else
    { // get the second root
      t = - ( D + Q ) / ( 2 * P );
      if ( 0.0 <= t && t <= 1.0 ) // Proper value of t found
        return t;
    }
  } // end else D > 0.0
  else // No double roots
    error("GetEdgeIntersectParamEllipsoid : D < 0 => No double solution exist. \n" );

  // If no proper value of t returned by now, return error
  return -1.0;
}
//=============================================================================
//
//
//
//
//=============================================================================
double
EllipsoidVolume(
const ellipsoid_t* ellipsoid )
{
  return 4. / 3. * PI * ellipsoid->SemiAxes[ 1 ]
                      * ellipsoid->SemiAxes[ 2 ]
                      * ellipsoid->SemiAxes[ 3 ];
}
//=============================================================================
/** Compute the inverse of the inertia tensor.
*/
//=============================================================================
void
EllipsoidInertiaTensorInv(
const ellipsoid_t*  ellipsoid,
      double        density,
      double**      InertiaTensorInv )
{
  // compute inertia tensor
  InertiaTensorInv[1][1] = pow(ellipsoid->SemiAxes[2],2) +
                        pow(ellipsoid->SemiAxes[3],2);
  InertiaTensorInv[2][2] = pow(ellipsoid->SemiAxes[1],2) +
                        pow(ellipsoid->SemiAxes[3],2);
  InertiaTensorInv[3][3] = pow(ellipsoid->SemiAxes[2],2) +
                        pow(ellipsoid->SemiAxes[1],2);

  double volume = EllipsoidVolume(ellipsoid);

  InertiaTensorInv[1][1] *= density * volume / 5.;
  InertiaTensorInv[2][2] *= density * volume / 5.;
  InertiaTensorInv[3][3] *= density * volume / 5.;
  
  /// invert it
  InertiaTensorInv[1][1] = 1. / InertiaTensorInv[1][1];
  InertiaTensorInv[2][2] = 1. / InertiaTensorInv[2][2];
  InertiaTensorInv[3][3] = 1. / InertiaTensorInv[3][3];
}
/*=============================================================================
PURPOSE
compute ellipsoid - plane algebraic distance, also return the closest point on the ellipsoid.

If P = (x,y,z) or (X,Y,Z) is a point,
  F_e( P ) = ( X / e_x )^2 + ( Y / e_y )^2 + ( Z / e_z )^2 - 1 = 0
is ellipsoid equation in its own frame,
  F_p( P ) = a * x + b * y + c * z + d = 0
is plane equation in the physical frame

Then, one condition for the minimum distance is
  n_p = coef * n_e, where
  n_p = (a,b,c), normal to plane vector rotated to be in the ellipsoid frame
  n_e = \nabla F_e = 2 * ( X / e_x^2 + Y / e_y^2 + Z / e_z^2 )
this condition leads to the coordinates of the point satisfying the condition on F_e
  P_e = ( a * e_x^2 , b * e_y^2 , c * e_z^2 ) / coef

The coef can be solved for by using F_e( P_e ) = 0, thus
  coef = +/- sqrt( ( a * e_x )^2 + ( b * e_y )^2 + ( c * e_z )^2 )
because of the symmetry of the ellipsoid, we have 2 solutions : the closest
and the farest point from the plane, on the ellipsoid.

Finally, the algebraic distance is the shortest of F_p( P_e ) / norm( n_p ).
============================================================================*/
double
EllipsoidPlaneDist(
const ellipsoid_t*  ellipsoid,
const double        ParPos[4],
      double**      RotMat,
const equation_t*   PlaneEquation )
{
  double vector[4]; // contains the rotated plane normal vector
  
  // rotation of the normal to plane to get in the particle frame
  gemtv3( RotMat, PlaneEquation->coef, vector );

  // solution point coordinates, there are 2 different possible solutions
  double point[2][4], dist[2], signed_coef;
  
  // coef in front of P_e : 1 / coef
  signed_coef = sqrt( pow( vector[1] * ellipsoid->SemiAxes[1], 2 ) +
                      pow( vector[2] * ellipsoid->SemiAxes[2], 2 ) +
                      pow( vector[3] * ellipsoid->SemiAxes[3], 2 ) );

  signed_coef = 1. / signed_coef;
  
  // compute distances for the 2 solutions
  for ( int i = 0 ; i <= 1 ; i++ )
  {
    // point coordinates on ellipsoid in ellipsoid frame
    for ( int j = 1 ; j <= 3 ; j++ )
      point[i][j] = signed_coef * vector[j] * pow( ellipsoid->SemiAxes[j], 2 );

    // express this point in the physical frame
    ChangeFrameReverse( RotMat, ParPos, point[i] );

    // compute plane equation applied to the point in the physical frame
    dist[i] = PlaneEquation->coef[1] * point[i][1] +
              PlaneEquation->coef[2] * point[i][2] +
              PlaneEquation->coef[3] * point[i][3] +
              PlaneEquation->coef[0];
    
    // divide by plane normal vector norm to get ditance
    dist[i] /= sqrt( pow( PlaneEquation->coef[1], 2 ) +
                     pow( PlaneEquation->coef[2], 2 ) +
                     pow( PlaneEquation->coef[3], 2 ) );
  
    signed_coef *= - 1.; // change sign to get second point
                         
//    debug("dist = %g\n", dist[i]);
  }

  // return minimum algebraic distance
  // thus, if it is negative, we know the particle is outside
  if ( dist[0] < dist[1] ) return dist[0];
  else                     return dist[1];  
}
