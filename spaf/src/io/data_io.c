/** \file
 Basic data input and output.
 
 Read/write an array from/to a file. The format of a file defined by a header followed by the data.
 The header could be for instance:\n
 number = 123\n
 format = binary\n
 byte_order = little\n
 type = float\n
 which means that the file contains an array of 123 float written in binary with little endian byte ordering.
 The format can be 'ascii' or 'binary'.
 The byte order can be 'little' or 'big' and must only be present when the format is binary.
 The type can be any C data type, but only functions for reading or writing 'int' or 'double' are provided here.
 Finally, after the header come the data of the array, with one line per entry in ascii format.
 
 Platform byte order can be checked too.
*/
#include "parse.h"
#include "logging.h"
#include "linalg.h"
#include "strings.h"
#include "file_system.h"
#include "data_io.h"

//=============================================================================
/// Check byteorder, return true if little endian, false otherwise.
//=============================================================================
bool
LittleEndian(void)
{
  int i = 1;
  char *p = (char *)&i;

  if ( p[0] == 1 ) return true;
  else return false;
}
//=============================================================================
/** Read an array from a file.
 
 This macro can generate a function for any data type.
 */
//=============================================================================
#define ReadArray(type,type_min,type_max) \
void \
ReadArray_##type( \
const char*   FileName, \
      int*    NbOfEntries, \
      type**  data ) \
{ \
  *NbOfEntries = 0; \
   \
  /* return if no filename */ \
  if ( FileName == NULL ) return; \
  \
  FILE* FileId = fopen(FileName, "r"); \
   \
  /* return if no file */ \
  if ( FileId == NULL ) return; \
   \
  /* read number of elements */ \
  TokenizeFileLine(FileId); \
  assert_error( TokenIs("number",1), "No number provided" ); \
  \
  int _NbOfEntries = Token2int(3, '[', 1, INT_MAX, ']'); \
  \
  /* allocate memory */ \
  assert_error( (*data) == NULL, \
          "Data array has to be NULL, because it'll be allocated"); \
  \
  (*data) = (type*) malloc(_NbOfEntries*sizeof(type)); \
  \
  /* read format */ \
  TokenizeFileLine(FileId); \
  assert_error( TokenIs("format",1), "No format provided" ); \
  bool BinaryFormat = true; \
  if ( TokenIs("ascii",3) ) BinaryFormat = false; \
  \
  /* read type of data */ \
  TokenizeFileLine(FileId); \
  assert_error( TokenIs("type",1), "No type provided" ); \
  assert_error( TokenIs(#type,3), "Type is not "#type ); \
  \
  /* for binary format, we need to read the type of data and byte order */ \
  if ( BinaryFormat ) \
  { \
    /* read byte order */ \
    TokenizeFileLine(FileId); \
    assert_error( TokenIs("byte_order",1), "No byte_order provided" ); \
    assert_error( TokenIs("little",3), "Byte order conversion not implemented" ); \
    \
    /* allocate memory and read data */ \
    fread((*data), sizeof(type), _NbOfEntries, FileId); \
  } \
  else \
  { \
    size_t counter = 0; \
     \
    while ( counter < (size_t) _NbOfEntries ) \
    { \
      /* parse line and get number of elements on the line */ \
      /* and check we have the correct number of elements */ \
      assert_error( TokenizeFileLine(FileId) == 1, "More than 1 datum per line"); \
       \
      /* read entry on the line */ \
      (*data)[counter] = Token2##type(1, '[', (type_min), (type_max), ']'); \
      counter++; \
    } \
  } \
   \
  /* Everything went fine, fill return variables */ \
  *NbOfEntries = _NbOfEntries; \
   \
  fclose(FileId); \
} \

ReadArray(int,0,INT_MAX)
ReadArray(double,-DBL_MAX,DBL_MAX)


//=============================================================================
/** Read and copy an array of double from a file to another array of double.
 
 This function is usefull when a certain number of entries is expected.
 */
//=============================================================================
void
ReadAndCopyArray_double(
const char*   FileName,     ///< file name to read from
const int     NbOfEntries,  ///< number of entries expected
      double* data )        ///< array of data
{
  // check arguments
  assert_error( data != NULL, "data is NULL");
  
  int number = 0;
  double *array = NULL;

  ReadArray_double( FileName, &number, &array );

  assert_error( number == NbOfEntries, "Bad number of entries");
  
  copy( number, array, data );
  
  free(array);
}
//=============================================================================
/** Write an array of int to a file.
 */
//=============================================================================
void
WriteArray_int(
const char* FileName,     ///< name of the file to write to
const int   NbOfEntries,  ///< number of entries in data
const int*  data,         ///< array of data
const bool  BinaryFormat) ///< format to write to
{
  FILE* FileId = fopen(FileName, "w");
  assert_error(FileId, "Failed to open %s""\n", FileName);

  // write header
  fprintf(FileId,"number = "INT_FMT"\n",NbOfEntries);
  char format[10] = "ascii";
  if ( BinaryFormat ) strcpy(format, "binary");
  fprintf(FileId,"format = %s""\n",format);
  fprintf(FileId,"type = int""\n"); 
  
  if ( BinaryFormat )
  {
    if ( LittleEndian() == true )
      fprintf(FileId,"byte_order = little""\n");
    else
      fprintf(FileId,"byte_order = big""\n");

    fwrite(data, sizeof(int), NbOfEntries, FileId);
  }
  else
  {
    for ( int i = 0 ; i < NbOfEntries ; i++ )
      fprintf(FileId, INT_FMT "\n", data[i]);
  }
  
  fclose(FileId);
}
//=============================================================================
/** Write an array of double to a file.
 */
//=============================================================================
void
WriteArray_double(
const char*   FileName,     ///< name of the file to write to
const int     NbOfEntries,  ///< number of entries in data
const double* data,         ///< array of data
const bool    BinaryFormat, ///< format to write
const bool    WriteAsDouble ) ///< write as double or float ?
{
  FILE* FileId = fopen(FileName, "w");
  assert_error(FileId, "Failed to open %s""\n", FileName);
  
  fprintf(FileId,"number = "INT_FMT"\n",NbOfEntries);
  
  char format[10] = "ascii";
  if ( BinaryFormat ) strcpy(format, "binary");
  fprintf(FileId,"format = %s""\n",format);
  
  char type[10] = "float";
  if ( WriteAsDouble ) strcpy(type, "double");
  fprintf(FileId,"type = %s""\n",type);

  if ( BinaryFormat )
  {    
    if ( LittleEndian() == true )
      fprintf(FileId,"byte_order = little""\n");
    else
      fprintf(FileId,"byte_order = big""\n");
    
    if ( WriteAsDouble )
      fwrite(data, sizeof(double), NbOfEntries, FileId);
    else
    {    
      for ( int i = 0 ; i < NbOfEntries ; i++ )
      {
        float datum = (float) data[i];
        fwrite(&datum, sizeof(float), 1, FileId);
      }      
    }
  }
  else
  {    
    if ( WriteAsDouble )
      for ( int i = 0 ; i < NbOfEntries ; i++ )
        fprintf(FileId, DBL_FMT_FULL "\n", data[i]);
    else
      for ( int i = 0 ; i < NbOfEntries ; i++ )
        fprintf(FileId, DBL_FMT "\n", data[i]);
  }
  
  fclose(FileId);
}
