// Copyright 2017 Fred Hutchinson Cancer Research Center

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utilsnetcdf.h"


bool UtilsNetcdf::inquireVariable(nc_type *xtype, int ncid, const char *variable)
{

  int nret;
  int varid;
  nret = nc_inq_varid(ncid, variable, &varid);
  FALSE_ON_NETCDF_ERROR(nret, "INFO could not find variable %s\n", variable);


  nret = nc_inq_var(ncid, varid, 0, xtype, 0, 0, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);

  return true;
}



bool UtilsNetcdf::inquireMatrixVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable, nc_type expectedXtype)
{
  int nret;
  nret = nc_inq_varid(ncid, variable, varid);
  FALSE_ON_NETCDF_ERROR(nret, "INFO could not find variable %s\n", variable);


  nc_type xtype;
  const int expectedDims = 2;
  int nDims;

  nret = nc_inq_var(ncid, *varid, 0, &xtype, &nDims, 0, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);

  if (expectedXtype != xtype) {
    fprintf(stderr, "ERROR unexpected type %i for variable %s\n", xtype, variable);
    return false;
  }


  if (expectedDims != nDims) {
    fprintf(stderr, "ERROR unexpected dimensions %i for variable %s\n", nDims, variable);
    return false;
  }


  int dimids[expectedDims];
  nret = nc_inq_var(ncid, *varid, 0, 0, 0, dimids, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);

  size_t sizes[expectedDims];
  for (int i = 0; i < expectedDims; i++) {
    nret = nc_inq_dim(ncid, dimids[i], 0, sizes+i);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about dimension %i\n", dimids[i]);
  }


  *nRows = sizes[0];
  *nCols = sizes[1];

  return true;
}


bool UtilsNetcdf::inquireMatrix8bitVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable)
{
  return inquireMatrixVariable(varid, nRows, nCols, ncid, variable, NC_BYTE);
}


bool UtilsNetcdf::inquireMatrixTextVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable)
{
  return inquireMatrixVariable(varid, nRows, nCols, ncid, variable, NC_CHAR);
}


bool UtilsNetcdf::inquireMatrixFloatVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable)
{
  return inquireMatrixVariable(varid, nRows, nCols, ncid, variable, NC_FLOAT);
}

bool UtilsNetcdf::inquireMatrixDoubleVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable)
{
  return inquireMatrixVariable(varid, nRows, nCols, ncid, variable, NC_DOUBLE);
}



bool UtilsNetcdf::inquireVectorVariable(int *varid, size_t *nRows, int ncid, const char *variable, nc_type expectedXtype)
{
  int nret;
  nret = nc_inq_varid(ncid, variable, varid);
  FALSE_ON_NETCDF_ERROR(nret, "INFO could not find variable %s\n", variable);


  nc_type xtype;
  const int expectedDims = 1;
  int nDims;

  nret = nc_inq_var(ncid, *varid, 0, &xtype, &nDims, 0, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);

  if (expectedXtype != xtype) {
    fprintf(stderr, "ERROR unexpected type %i for variable %s\n", xtype, variable);
    return false;
  }


  if (expectedDims != nDims) {
    fprintf(stderr, "ERROR unexpected dimensions %i for variable %s\n", nDims, variable);
    return false;
  }


  int dimids[expectedDims];
  nret = nc_inq_var(ncid, *varid, 0, 0, 0, dimids, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);

  size_t sizes[expectedDims];
  for (int i = 0; i < expectedDims; i++) {
    nret = nc_inq_dim(ncid, dimids[i], 0, sizes+i);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about dimension %i\n", dimids[i]);
  }


  *nRows = sizes[0];

  return true;
}



bool UtilsNetcdf::inquireVectorIntVariable(int *varid, size_t *nRows, int ncid, const char *variable)
{
  return inquireVectorVariable(varid, nRows, ncid, variable, NC_INT);
}


bool UtilsNetcdf::inquireVector8bitVariable(int *varid, size_t *nRows, int ncid, const char *variable)
{
  return inquireVectorVariable(varid, nRows, ncid, variable, NC_BYTE);
}


bool UtilsNetcdf::inquireVectorDoubleVariable(int *varid, size_t *nRows, int ncid, const char *variable)
{
  return inquireVectorVariable(varid, nRows, ncid, variable, NC_DOUBLE);
}




bool UtilsNetcdf::inquireMatrixStorage(bool *SamplebySNPs, int varid, int ncid)
{
  int nret;
  const int expectedDims = 2;

  int dimids[expectedDims];
  nret = nc_inq_var(ncid, varid, 0, 0, 0, dimids, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable id %i\n", varid);


  char names[expectedDims][NC_MAX_NAME+1]; //row major order


  for (int i = 0; i < expectedDims; i++) {
    nret = nc_inq_dim(ncid, dimids[i], names[i], 0);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about dimension %i\n", dimids[i]);
  }

  if ( (0 == strcmp(names[0], "Samples")) && (0 == strcmp(names[1], "SNPs")) ) {
    *SamplebySNPs = true;
    return true;
  }

  if ( (0 == strcmp(names[1], "Samples")) && (0 == strcmp(names[0], "SNPs")) ) {
    *SamplebySNPs = false;
    return true;
  }

  fprintf(stderr, "one or both dimension names are unknown %s, %s", names[0], names[1]);

  return false;
}






//try to handle unsigned char and int types transparently ...
bool UtilsNetcdf::load(sspt_Array< int > *vec, int ncid, const char *variable)
{
  int nret;

  int varid;
  nret = nc_inq_varid(ncid, variable, &varid);
  FALSE_ON_NETCDF_ERROR(nret, "INFO could not find variable %s\n", variable);


  nc_type xtype;
  const int expectedDims = 1;
  int nDims;

  nret = nc_inq_var(ncid, varid, 0, &xtype, &nDims, 0, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);

  if (!(NC_BYTE == xtype || NC_INT == xtype)) {
    fprintf(stderr, "ERROR unexpected type %i for variable %s\n", xtype, variable);
    return false;
  }


  if (1 != nDims) {
    fprintf(stderr, "ERROR unexpected dimensions %i for variable %s\n", nDims, variable);
    return false;
  }


  int dimids[expectedDims];
  nret = nc_inq_var(ncid, varid, 0, 0, 0, dimids, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);


  size_t sizes[expectedDims];
  size_t total = 1;
  for (int i = 0; i < expectedDims; i++) {
    nret = nc_inq_dim(ncid, dimids[i], 0, sizes+i);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about dimension %i\n", dimids[i]);
    total *= sizes[i];
  }

  if (NC_BYTE == xtype) {
    unsigned char *array = new unsigned char[total];
    nret = nc_get_var_uchar(ncid, varid, array);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read variable %s\n", variable);
    (*vec) = sspt_Array<int>( sizes[0] );
    for (size_t i = 0; i < sizes[0]; i++) {
      (*vec)[i] = array[i];
    }
    delete[] array;
  }
  else if (NC_INT == xtype) {
    int *array = new int[total];
    nret = nc_get_var_int(ncid, varid, array);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read variable %s\n", variable);
    (*vec) = sspt_Array<int>( sizes[0] );
    for (size_t i = 0; i < sizes[0]; i++) {
      (*vec)[i] = array[i];
    }
    delete[] array;

  }

  return true;
}




bool UtilsNetcdf::load(sspt_Array<const char *> *vec, char **buffer, int ncid, const char *variable)
{
  int nret;

  int varid;
  nret = nc_inq_varid(ncid, variable, &varid);
  FALSE_ON_NETCDF_ERROR(nret, "INFO could not find variable %s\n", variable);


  nc_type xtype;
  const int expectedDims = 2;
  int nDims;

  nret = nc_inq_var(ncid, varid, 0, &xtype, &nDims, 0, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);

  if (NC_STRING == xtype) {
    return loadString(vec, buffer, ncid, variable);
  }


  if (NC_CHAR != xtype) {
    fprintf(stderr, "ERROR unexpected type %i for variable %s\n", xtype, variable);
    return false;
  }


  if (2 != nDims) {
    fprintf(stderr, "ERROR unexpected dimensions %i for variable %s\n", nDims, variable);
    return false;
  }


  int dimids[expectedDims];
  nret = nc_inq_var(ncid, varid, 0, 0, 0, dimids, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);


  size_t sizes[expectedDims];
  size_t total = 1;
  for (int i = 0; i < expectedDims; i++) {
    nret = nc_inq_dim(ncid, dimids[i], 0, sizes+i);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about dimension %i\n", dimids[i]);
    total *= sizes[i];
  }

  char *array = new char[total];
  nret = nc_get_var_text(ncid, varid, array);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read variable %s\n", variable);
  (*vec) = sspt_Array<const char *>( sizes[0] );

  //get strings in vector ...
  for (size_t i = 0; i < sizes[0]; i++) {
    (*vec)[i] = array+i*sizes[1];
  }
  *buffer = array;

  return true;
}



//load NetCDF4 strings
bool UtilsNetcdf::loadString(sspt_Array<const char *> *vec, char **buffer, int ncid, const char *variable)

{
  int nret;

  int varid;
  nret = nc_inq_varid(ncid, variable, &varid);
  CLOSE_ON_NETCDF_ERROR(nret, ncid, "ERROR could not find variable %s\n", variable);


  nc_type xtype;
  int nDims;

  nret = nc_inq_var(ncid, varid, 0, &xtype, &nDims, 0, 0);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about variable %s\n", variable);

  if (NC_STRING != xtype) {
    fprintf(stderr, "ERROR unexpected type %i for variable %s\n", xtype, variable);
    return false;
  }


  if (1 != nDims) {
    fprintf(stderr, "ERROR unexpected dimensions %i for variable %s\n", nDims, variable);
    return false;
  }


  int dimids[2];
  nret = nc_inq_var(ncid, varid, 0, 0, 0, dimids, 0);
  FALSE_ON_NETCDF_ERROR(nret,  "ERROR could not read information about variable %s\n", variable);


  size_t sizes[2];
  size_t total = 1;
  for (int i = 0; i < nDims; i++) {
    nret = nc_inq_dim(ncid, dimids[i], 0, sizes+i);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read information about dimension %i\n", dimids[i]);
    total *= sizes[i];
  }

  char **arrayOfStrings = new char*[total];
  nret = nc_get_var_string(ncid, varid, arrayOfStrings);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not read variable %s\n", variable);

  //figure out total size of memory for strings
  size_t n = 0;
  for (size_t i = 0; i < sizes[0]; i++) {
    n += strlen( arrayOfStrings[i] ) + 1;
  }

  //get strings in vector, and store in contigous allocated buffer
  vec->resize( sizes[0] );
  char *array = new char[n];
  char *current = array;
  for (size_t i = 0; i < sizes[0]; i++) {
    (*vec)[i] = strcpy(current, arrayOfStrings[i]);
    current += strlen(arrayOfStrings[i]) + 1;
  }

  //need to free netCDF/HDF string memory
  nret = nc_free_string( sizes[0], arrayOfStrings );
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not free HDF NC_STRING memory\n");

  delete[] arrayOfStrings;

  *buffer = array;

  return true;
}




// Definitions
template<>
int get_var<char>(int ncid, int varid, char *tp)
{
  return nc_get_var_text(ncid, varid, tp);
}

template<>
int get_var<unsigned char>(int ncid, int varid, unsigned char *tp)
{
  return nc_get_var_uchar(ncid, varid, tp);
}

template<>
int get_var<signed char>(int ncid, int varid, signed char *tp)
{
  return nc_get_var_schar(ncid, varid, tp);
}

template<>
int get_var<int>(int ncid, int varid, int *tp)
{
  return nc_get_var_int(ncid, varid, tp);
}

template<>
int get_var<double>(int ncid, int varid, double *tp)
{
  return nc_get_var_double(ncid, varid, tp);
}

template<>
int get_var<float>(int ncid, int varid, float *tp)
{
  return nc_get_var_float(ncid, varid, tp);
}


template<>
int put_var<char>(int ncid, int varid, const char *tp)
{
  return nc_put_var_text(ncid, varid, tp);
}

template<>
int put_var<unsigned char>(int ncid, int varid, const unsigned char *tp)
{
  return nc_put_var_uchar(ncid, varid, tp);
}

template<>
int put_var<signed char>(int ncid, int varid, const signed char *tp)
{
  return nc_put_var_schar(ncid, varid, tp);
}

template<>
int put_var<int>(int ncid, int varid, const int *tp)
{
  return nc_put_var_int(ncid, varid, tp);
}

template<>
int put_var<double>(int ncid, int varid, const double *tp)
{
  return nc_put_var_double(ncid, varid, tp);
}


template<>
int put_var<float>(int ncid, int varid, const float *tp)
{
  return nc_put_var_float(ncid, varid, tp);
}


template<>
int get_vara<char>(int ncid, int varid, const size_t start[], const size_t count[],  char *tp)
{
  return nc_get_vara_text(ncid, varid, start, count, tp);
}

template<>
int get_vara<unsigned char>(int ncid, int varid, const size_t start[], const size_t count[],  unsigned char *tp)
{
  return nc_get_vara_uchar(ncid, varid, start, count, tp);
}

template<>
int get_vara<signed char>(int ncid, int varid, const size_t start[], const size_t count[],  signed char *tp)
{
  return nc_get_vara_schar(ncid, varid, start, count, tp);
}


template<>
int get_vara<int>(int ncid, int varid, const size_t start[], const size_t count[],  int *tp)
{
  return nc_get_vara_int(ncid, varid, start, count, tp);
}

template<>
int get_vara<double>(int ncid, int varid, const size_t start[], const size_t count[],  double *tp)
{
  return nc_get_vara_double(ncid, varid, start, count, tp);
}

template<>
int get_vara<float>(int ncid, int varid, const size_t start[], const size_t count[],  float *tp)
{
  return nc_get_vara_float(ncid, varid, start, count, tp);
}



template<>
int put_vara<char>(int ncid, int varid, const size_t start[], const size_t count[],  const char *tp)
{
  return nc_put_vara_text(ncid, varid, start, count, tp);
}

template<>
int put_vara<unsigned char>(int ncid, int varid, const size_t start[], const size_t count[],  const unsigned char *tp)
{
  return nc_put_vara_uchar(ncid, varid, start, count, tp);
}

template<>
int put_vara<signed char>(int ncid, int varid, const size_t start[], const size_t count[],  const signed char *tp)
{
  return nc_put_vara_schar(ncid, varid, start, count, tp);
}

template<>
int put_vara<int>(int ncid, int varid, const size_t start[], const size_t count[],  const int *tp)
{
  return nc_put_vara_int(ncid, varid, start, count, tp);
}

template<>
int put_vara<double>(int ncid, int varid, const size_t start[], const size_t count[],  const double *tp)
{
  return nc_put_vara_double(ncid, varid, start, count, tp);
}

template<>
int put_vara<float>(int ncid, int varid, const size_t start[], const size_t count[],  const float *tp)
{
  return nc_put_vara_float(ncid, varid, start, count, tp);
}








