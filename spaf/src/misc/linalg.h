#ifndef LINALG_H
#define LINALG_H

void
permute(
const int     size,
const int*    permutation,
const double* x,
      double* y );

void
permuteinv(
const int     size,
const int*    permutation,
const double* x,
      double* y );

void
axpby(
const int     size,
const double  a,
const double* x,
const double  b,
      double* y );

void
axpby3(
const double  a,
const double* x,
const double  b,
double* y );

void
scal(
const int     size,
const double  coef,
      double* array_out );
  
void
scal3(
const double  coef,
      double* array_out );

double
dot(
const int     size,
const double* x,
const double* y );

double
dot3(
const double* x,
const double* y );

void
copy(
const int     size,
const double* x,
      double* y );

void
copy3(
const double* x,
double* y );

double
nrm2(
const int     size,
const double* x );

double
nrm23(
const double* x );

void
TensorProduct(
const int     size,
const double* x,
const double* y,
      double* z );

//=============================================================================
// MATRICES
//=============================================================================

void setZerosM3( double** mat );

//=============================================================================
// MATRICES - VECTORS COMPUTATIONS
//=============================================================================

void
gemv3(
      double** A,
const double b[4],
      double c[4] );

void gemtv3(
        double** A,
  const double b[4],
        double c[4] );

void
gemm3(
double** A,
double** B,
double** C );

void
gemmt3(
         double** A,
         double** B,
         double** C );

void
gemtm3(
         double** A,
         double** B,
         double** C );

//=============================================================================
// ROTATIONS AND FRAME CHANGE
//=============================================================================
void
RotationVectorToMatrix(
const double* vector,       
      double** matrix );

void
RotationMatrixToVector(
double** matrix,
double vector[4] );

void
ChangeFrame(
      double** RotMat,
const double ParPos[4],
      double point[4] );

void
ChangeFrameReverse(
      double** RotMat,
const double ParPos[4],
      double point[4] );

double
det3by3(
double** myMat );

double
inv3by3(
double** J,
double** J_inv );
  
double
VectorProduct(
const int dir,   
const double  a[4], 
const double  b[4], 
const double  c[4] );

double
VectorProduct2(
const int dir,   
const double  a[4], 
const double  b[4] );

void
ChangeMatrixBase(
double** InertiaTensorInv,
double** RotMat,
double** InertiaTensorInvPhysFrame );

#endif
