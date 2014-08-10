#ifndef DATA_IO_H
#define DATA_IO_H

#include "includes.h"

bool
LittleEndian(void);

void
ReadAndCopyArray_double(
const char*   FileName,
const int     NbOfEntries,
      double* data );

void
ReadArray_int(
const char* FileName,
      int*  NbOfEntries,
      int** data );

void
ReadArray_double(
const char*     FileName,
      int*      NbOfEntries,
      double**  data );

void
WriteArray_int(
const char* FileName,     
const int   NbOfEntries,  
const int*  data,         
const bool  BinaryFormat);

void
WriteArray_double(
const char*   FileName,       
const int     NbOfEntries,    
const double* data,           
const bool    BinaryFormat,   
const bool    FullPrecision );

#endif
