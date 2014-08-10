/** \file
 Linear algebra computations.
 
 This file can be used for computing linear algebra directly or as a wrapper to blas library. All the functions that act on non 3D arrays are 0 offset, but 3D algebra is 1 offset.
 */

#include "includes.h"
#include "logging.h"
#include "linalg.h"
#include "memory.h"

// blas library
#ifdef USE_ACCELERATE
# include "Accelerate/Accelerate.h"
# define USE_BLAS
#elif USE_MKL
# include "mkl_cblas.h"
# define USE_BLAS
#elif USE_ACML
# include "acml.h"
# define USE_BLAS
# define cblas_dscal dscal
# define cblas_daxpy daxpy
# define cblas_ddot  ddot
# define cblas_dcopy dcopy
#endif

//=============================================================================
// Permute 2 vectors.
//=============================================================================
inline void
permute(
const int     size, ///< size of vectors
const int*    permutation, ///< size of vectors
const double* x,    ///< vector
      double* y )   ///< vector
{
#pragma omp parallel for \
default(none) \
shared(x,y,permutation)
  for ( int i = 0 ; i < size ; i++ )
    y[i] = x[ permutation[i] ];
}
//=============================================================================
// Permute inversely 2 vectors.
//=============================================================================
inline void
permuteinv(
const int     size, ///< size of vectors
const int*    permutation, ///< size of vectors
const double* x,    ///< vector
      double* y )   ///< vector
{
#pragma omp parallel for \
default(none) \
shared(x,y,permutation)
  for ( int i = 0 ; i < size ; i++ )
    y[ permutation[i] ] = x[i];
}
//=============================================================================
// y = a*x + b*y
//=============================================================================
inline void
axpby(
const int     size, ///< size of vectors
const double  a,    ///< coefficient
const double* x,    ///< vector
const double  b,    ///< coefficient
      double* y )   ///< vector
{
#ifdef USE_BLAS
  cblas_dscal(size, b, y, 1);
  cblas_daxpy(size, a, x, 1, y, 1);
#else
#pragma omp parallel for \
default(none) \
shared(x,y)
  for ( int i = 0 ; i < size ; i++ )
    y[ i ] = a * x[ i ] + b * y[ i ];
#endif
}
//=============================================================================
// y = a*x + b*y in 3D
//=============================================================================
inline void
axpby3(
const double  a,    ///< coefficient
const double* x,    ///< vector
const double  b,    ///< coefficient
      double* y )   ///< vector
{
  y[1] = a * x[1] + b * y[1];
  y[2] = a * x[2] + b * y[2];
  y[3] = a * x[3] + b * y[3];
}
//=============================================================================
// x = a*x
//=============================================================================
inline void
scal(
const int     size, ///< size of vectors
const double  a,    ///< coefficient
      double* x )   ///< vector
{
#ifdef USE_BLAS
  cblas_dscal( size, a, x, 1 );
#else
#pragma omp parallel for \
default(none) \
shared(x)
  for ( int i = 0 ; i < size ; i++ )
    x[i] *= a;
#endif
}
//=============================================================================
// x = a*x in 3D
//=============================================================================
inline void
scal3(
const double  a,  ///< coefficient
      double* x ) ///< vector
{
  x[1] *= a;
  x[2] *= a;
  x[3] *= a;
}
//=============================================================================
// Scalar product
//=============================================================================
inline double
dot(
const int     size, ///< size of vectors
const double* x,    ///< vector
const double* y )   ///< vector
{
#ifdef USE_BLAS
  return cblas_ddot(size, x, 1, y, 1);
#else
  double ddot = 0.;
  
#pragma omp parallel for \
default(none) \
shared(x,y) \
reduction(+:ddot)
  for ( int i = 0 ; i < size ; i++ )
    ddot += x[i] * y[i];
  
  return ddot;
#endif
}
//=============================================================================
// Scalar product in 3D
//=============================================================================
inline double
dot3(
const double* x,  ///< vector
const double* y ) ///< vector
{
  return x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
}
//=============================================================================
// y = x
//=============================================================================
inline void
copy(
const int     size, ///< size of vectors
const double* x,    ///< vector
      double* y )   ///< vector
{
#ifdef USE_BLAS
  cblas_dcopy( size, x, 1, y, 1);
#else
#pragma omp parallel for \
default(none) \
shared(x,y)
  for ( int i = 0 ; i < size ; i++ )
    y[i] = x[i];
#endif
}
//=============================================================================
// y = x in 3D
//=============================================================================
inline void
copy3(
const double* x,  ///< vector
      double* y ) ///< vector
{
  y[1] = x[1];
  y[2] = x[2];
  y[3] = x[3];
}
//=============================================================================
// L2 Norm of a vector
//=============================================================================
inline double
nrm2(
const int     size, ///< size of vectors
const double* x )   ///< vector
{
#ifdef USE_BLAS
  return cblas_dnrm2( size, x, 1 );
#else
  return sqrt( dot(size, x, x) );
#endif
}
//=============================================================================
// L2 Norm of a vector in 3D
//=============================================================================
inline double
nrm23(
const double* x ) ///< vector
{
  return sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
}
//=============================================================================
// z[i] = x[i] * y[i] for all i, 0-offset.
//=============================================================================
inline void
TensorProduct(
const int     size, ///< size of vectors
const double* x,    ///< vector
const double* y,    ///< vector
      double* z )   ///< vector
{
#pragma omp parallel for \
default(none) \
shared(x,y,z)
  for ( int i = 0 ; i < size ; i++ )
    z[i] = x[i] * y[i];
}
//=============================================================================
//=============================================================================
// MATRICES - VECTORS COMPUTATIONS, 1 offset.
//=============================================================================
//=============================================================================
void
setZerosM3(
double** mat )
{
  for ( int i = 1 ; i <= 3 ; i++ )
  for ( int j = 1 ; j <= 3 ; j++ )
    mat[ i ][ j ] = 0.;
}
//=============================================================================
/// Matrix vector product
//=============================================================================
void
gemv3(
      double** A,
const double b[4],
      double c[4] )
{
  for ( int i = 1 ; i <= 3 ; i++ )
  {
    c[ i ] = 0.;
    
    for ( int j = 1 ; j <= 3 ; j++ )
      c[ i ] += A[ i ][ j ] * b[ j ];
  }
}
//=============================================================================
/// Matrix transpose vector product
//=============================================================================
void gemtv3(
        double** A,
  const double b[4],
        double c[4] )
{  
  for ( int i = 1 ; i <= 3 ; i++ )
  {
    c[ i ] = 0.;
    
    for ( int j = 1 ; j <= 3 ; j++ )
      c[ i ] += A[ j ][ i ] * b[ j ];
  }
}
//=============================================================================
/// Matrix matrix product.
//=============================================================================
void gemm3(
double** A,
double** B,
double** C )
{ 
  for ( int i = 1 ; i <= 3 ; i++ )
  for ( int j = 1 ; j <= 3 ; j++ )
  {
    C[ i ][ j ] = 0.;
    for ( int k = 1 ; k <= 3 ; k++ )
      C[ i ][ j ] += A[ i ][ k ] * B[ k ][ j ];
  }
}
//=============================================================================
/// Matrix matrix transposed product.
//=============================================================================
void gemmt3(
double** A,
double** B,
double** C )
{ 
  for ( int i = 1 ; i <= 3 ; i++ )
  for ( int j = 1 ; j <= 3 ; j++ )
  {
    C[ i ][ j ] = 0.;
    for ( int k = 1 ; k <= 3 ; k++ )
      C[ i ][ j ] += A[ i ][ k ] * B[ j ][ k ];
  }
}
//=============================================================================
/// Matrix transpose matrix product
//=============================================================================
void gemtm3(
double** A,
double** B,
double** C )
{ 
  for ( int i = 1 ; i <= 3 ; i++ )
  for ( int j = 1 ; j <= 3 ; j++ )
  {
    C[ i ][ j ] = 0.;
    for ( int k = 1 ; k <= 3 ; k++ )
      C[ i ][ j ] += A[ k ][ i ] * B[ k ][ j ];
  }
}
//=============================================================================
//=============================================================================
// ROTATIONS AND FRAME CHANGE
//=============================================================================
//=============================================================================
//=============================================================================
/// Compute the rotation matrix from a rotation vector.
/**
 Compute the rotation matrix from the rotation vector, i.e. a vector whose norm
 is the angle of rotation and direction is the axis of rotation.
 We use the Rodrigues formula, 

       | Ux²+Ca*(1-Ux²)        Ux*Uy*(1-Ca)-Uz*Sa     Ux*Uz*(1-Ca)+Uy*Sa |
 [R] = | Ux*Uy*(1-Ca)+Uz*Sa    Uy²+Ca*(1-Uy²)         Uy*Uz*(1-Ca)-Ux*Sa |
       | Ux*Uz*(1-Ca)-Uy*Sa    Uy*Uz*(1-Ca)+Ux*Sa     Uz²+Ca*(1-Uz²)     |
 
*/
//=============================================================================
void
RotationVectorToMatrix(
const double* vector,   ///< rotation vector
      double** matrix ) ///< rotation matrix
{                                                                     
  double
  angle = nrm23( vector ); // angle of rotation = norm of rotation vector

  // initialize to zero
  setZerosM3(matrix);
  
  // if angle is zero, we can't get the axis, we directly set the matrix to Id
  if ( angle == 0. )
  {
    for ( int i = 1 ; i <= 3 ; i++ ) matrix[i][i] = 1.;
    return;
  }
  
  // compute the unit vector of the rotation axis
  double axis[4];
  for ( int i = 1 ; i <= 3 ; i++ ) axis[i] = vector[i] / angle;
  
  // rotation matrix
  // R = cos( angle ) Id + ( 1 - cos( angle ) ) R_a + sin( angle ) R_s
  
  double
  Cos = cos( angle ),
  Sin = sin( angle );
  
  // cos( angle ) Id terms
  for ( int i = 1 ; i <= 3 ; i++ ) matrix[i][i] = Cos;
  
  // sin( angle ) R_s terms (off diagonal matrix)
  matrix[1][2] = - Sin * axis[3];
  matrix[1][3] = + Sin * axis[2];
  matrix[2][3] = - Sin * axis[1];
  matrix[2][1] = - matrix[1][2];
  matrix[3][1] = - matrix[1][3];
  matrix[3][2] = - matrix[2][3];
  
  // ( 1 - cos( angle ) ) R_a terms
  for ( int i = 1 ; i <= 3 ; i++ )
  for ( int j = 1 ; j <= 3 ; j++ )
    matrix[i][j] += ( 1. - Cos ) * axis[i] * axis[j];
}
//=============================================================================
/// Compute the rotation vector from a rotation matrix.
/**
 */
//=============================================================================
void
RotationMatrixToVector(
double**  matrix,     ///< rotation matrix
double    vector[4] ) ///< rotation vector
{
  // trace = 1 + 2 * cos ( angle )
  double angle = acos( ( matrix[1][1] + matrix[2][2] + matrix[3][3] - 1.) / 2. );
  
  // if the angle is zero, we can't compute the axis, then we set it directly.
  if ( angle == 0. )
  {
    scal3(0.,vector);
    return;
  }
  
  // two cases according to singularity
  if ( sin( angle ) != 0. )
  {
    // get unit vector from :
    // M - MT = 2 * sin ( angle ) * |  0  -n3  n2 |
    //                              |  n3  0  -n1 |
    //                              | -n2  n1  0  |
    
    double coef = angle / ( 2. * sin( angle ) );
    
    vector[1] = ( matrix[3][2] - matrix[2][3] ) * coef;
    vector[2] = ( matrix[1][3] - matrix[3][1] ) * coef;
    vector[3] = ( matrix[2][1] - matrix[1][2] ) * coef;    
  }
  else
  {
    // first singularity is angle = 0, already done
    // second one is angle = Pi

    // go back to Rodrigues formula : 
    // R = -Id + 2 * R_a where R_a = U Ut
    // compute positive components only
    vector[1] = angle * sqrt ( ( pow(matrix[1][1],2) + 1. ) / 2. );
    vector[2] = angle * sqrt ( ( pow(matrix[2][2],2) + 1. ) / 2. );
    vector[3] = angle * sqrt ( ( pow(matrix[3][3],2) + 1. ) / 2. );

    // treshold
    double epsilon = DBL_EPSILON;
    
    // deal with the signs
    
    if ( fabs(vector[1]) <  epsilon &&
         fabs(vector[2]) >= epsilon &&
         fabs(vector[3]) >= epsilon )
    {
		  if ( matrix[2][3] > 0 ) vector[2] = - vector[2];
		}
    else if ( fabs(vector[2]) <  epsilon &&
              fabs(vector[3]) >= epsilon )
    {
		  if ( matrix[1][3] <= 0 ) vector[3] = - vector[3];
		}
    else if ( fabs(vector[3]) < epsilon )
    {
		  if ( matrix[1][2] <= 0 ) vector[1] = - vector[1];
		}
  }  
}

//=============================================================================
/// Computes the position of a point from current to target frame.
/**
 The target frame origin is defined in the current frame by the vector origin,
 the rotation of the target frame in the current frame is defined in the
 rotation matrix.
 The point to be expressed in the target frame is translated to the origin of
 it, then rotated
 The array point is overwritten.
*/
//=============================================================================
void
ChangeFrame(
      double** matrix,  ///< rotation matrix
const double origin[4],     ///< origin of target frame
      double point[4] )     ///< position of a point
{
  // translation to target frame ogirin
  double temp[4] = {0.,
  point[1] - origin[1],
  point[2] - origin[2],
  point[3] - origin[3]};
  
  // rotation
  gemtv3( matrix, temp, point);
}
//=============================================================================
/// Computes the position of a point from target to current frame.
//=============================================================================
void
ChangeFrameReverse(
      double** matrix,  ///< rotation matrix
const double origin[4],     ///< origin of target frame
      double point[4] )     ///< position of a point
{
  double temp[4] = {0.,0.,0.,0.};
  
  // inverse rotation
  gemv3( matrix, point, temp);
  
  // translation to current frame origin
  for ( int dir = 1 ; dir <= 3 ; dir++ )
    point[dir] = temp[dir] + origin[dir];
}
//=============================================================================
/** Compute the determinant of a 3 by 3 matrix.
 */
//=============================================================================
double
det3by3(
double **myMat )
{
  return
  myMat[ 1 ][ 1 ] * myMat[ 2 ][ 2 ] * myMat[ 3 ][ 3 ] +
  myMat[ 2 ][ 1 ] * myMat[ 3 ][ 2 ] * myMat[ 1 ][ 3 ] +
  myMat[ 3 ][ 1 ] * myMat[ 1 ][ 2 ] * myMat[ 2 ][ 3 ] -
  (	myMat[ 1 ][ 3 ] * myMat[ 2 ][ 2 ] * myMat[ 3 ][ 1 ] +
   myMat[ 2 ][ 3 ] * myMat[ 3 ][ 2 ] * myMat[ 1 ][ 1 ] +
   myMat[ 3 ][ 3 ] * myMat[ 1 ][ 2 ] * myMat[ 2 ][ 1 ] );
}
//=============================================================================
/** Inverte a 3 by 3 matrix. 
 */
//=============================================================================
double
inv3by3(
double** J,       ///< a matrix
double** J_inv )  ///< its inverse
{
  double det = det3by3( J );
  
  assert_error( det != 0., "Determinant is 0");
  
  J_inv[1][1] = (J[2][2]*J[3][3]-J[2][3]*J[3][2]) / det;
  J_inv[1][2] = (J[1][3]*J[3][2]-J[1][2]*J[3][3]) / det;
  J_inv[1][3] = (J[1][2]*J[2][3]-J[1][3]*J[2][2]) / det;
  J_inv[2][1] = (J[2][3]*J[3][1]-J[2][1]*J[3][3]) / det;
  J_inv[2][2] = (J[1][1]*J[3][3]-J[1][3]*J[3][1]) / det;
  J_inv[2][3] = (J[1][3]*J[2][1]-J[1][1]*J[2][3]) / det;
  J_inv[3][1] = (J[2][1]*J[3][2]-J[2][2]*J[3][1]) / det;
  J_inv[3][2] = (J[1][2]*J[3][1]-J[1][1]*J[3][2]) / det;
  J_inv[3][3] = (J[1][1]*J[2][2]-J[1][2]*J[2][1]) / det;
  
  return det;
}
//=============================================================================
/** Compute a  \f$ \times ( b - c ) \f$ .
 */
//=============================================================================
double
VectorProduct(
const int     dir,    ///< direction
const double  a[4],   ///< vector
const double  b[4],   ///< vector
const double  c[4] )  ///< vector
{
  int next = dir % 3 + 1,
      prev = next % 3 + 1;
  
  return a[ next ] * ( b[ prev ] - c[ prev ] ) -
         a[ prev ] * ( b[ next ] - c[ next ] ); 
}
//=============================================================================
/** Computes the dir component of  \f$ a \times ( b - c ) \f$ .
 */
//=============================================================================
double
VectorProduct2(
const int     dir,    ///< direction
const double  a[4],   ///< vector
const double  b[4] )  ///< vector
{
  int next = dir % 3 + 1,
      prev = next % 3 + 1;
  
  return a[ next ] * b[ prev ] - a[ prev ] * b[ next ];
}
//==============================================================================
/// Update the inverse of inertia tensor in the physical frame.
//==============================================================================
void
ChangeMatrixBase(
double** InertiaTensorInv,
double** RotMat,
double** InertiaTensorInvPhysFrame )             
{
  // update inertia tensor in physical frame
  double **temp = NULL;
  AllocMdouble(3,3,temp);
  
  gemmt3(InertiaTensorInv,RotMat, temp);
  gemm3(RotMat,temp,InertiaTensorInvPhysFrame);
  
  FreeM(temp);
}
