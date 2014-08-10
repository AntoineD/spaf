#ifndef SMVP_LEGACY_H
#define SMVP_LEGACY_H

void
GenAdd(
const int     n,      
const int*    offset, 
const int*    columns,
const int*    rows,   
const double* entries,
const double* in,     
      double* out );

void
GenSuccAdd(
const int     n,      
const int*    offset, 
const int*    columns,
const double* entries,
const double* in,     
      double* out );

void
GenSuccAdd_T(
const int     nc,
const int     n,      
const int*    offset, 
const int*    columns,
const double* entries,
const double* in,     
      double* out );

void
GenSub(
const int     n,      
const int*    offset, 
const int*    columns,
const int*    rows,   
const double* entries,
const double* in,     
      double* out );

void
SymDiagFirstSuccAdd(
const int     n,
const int*    offset,
const int*    columns,
const double* entries,
const double* diagonal, 
const double* in,
      double* out );

#endif
