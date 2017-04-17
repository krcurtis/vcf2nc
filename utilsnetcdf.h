// Copyright 2017 Fred Hutchinson Cancer Research Center


#ifndef UTILSNETCDF_H
#define UTILSNETCDF_H


#include "sspt_array.h"
#include "sspt_tmatrix.h"
#include "netcdf.h"


#define FALSE_ON_NETCDF_ERROR(nret, ...)          \
  if (NC_NOERR != nret) {                         \
    fprintf(stderr, __VA_ARGS__);                 \
    fprintf(stderr, " -- netCDF error message %i : %s\n", nret, nc_strerror(nret)); \
    return false;                                 \
  }

#define CLOSE_ON_NETCDF_ERROR(nret, ncid, ...)    \
  if (NC_NOERR != nret) {                         \
    fprintf(stderr, __VA_ARGS__);                 \
    fprintf(stderr, " -- netCDF error message %i : %s\n", nret, nc_strerror(nret)); \
    ncclose(ncid);                                \
    return false;                                 \
  }

//! Assists in loading data from netCDF files
class UtilsNetcdf
{
 public:
  static bool inquireVariable(nc_type *xtype, int ncid, const char *variable);

  static bool inquireMatrixVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable, nc_type expectedXtype);
  static bool inquireMatrix8bitVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable);
  static bool inquireMatrixTextVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable);
  static bool inquireMatrixFloatVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable);
  static bool inquireMatrixDoubleVariable(int *varid, size_t *nRows, size_t *nCols, int ncid, const char *variable);


  static bool inquireVectorVariable(int *varid, size_t *nRows, int ncid, const char *variable, nc_type expectedXtype);
  static bool inquireVectorIntVariable(int *varid, size_t *nRows, int ncid, const char *variable);
  static bool inquireVector8bitVariable(int *varid, size_t *nRows, int ncid, const char *variable);
  static bool inquireVectorDoubleVariable(int *varid, size_t *nRows, int ncid, const char *variable);

  static bool inquireMatrixStorage(bool *SamplebySNPs, int varid, int ncid);

  static bool load(sspt_Array< int > *vec, int ncid, const char *variable);
  static bool load(sspt_Array<const char *> *vec, char **buffer, int ncid, const char *variable);
  static bool loadString(sspt_Array<const char *> *vec, char **buffer, int ncid, const char *variable);

#if 0
  void loadNetcdf(sspt_Array< std::string > *vec, const char *file, const char *variable);

  void loadNetcdf(sspt_TMatrix< int > *v, const char *file, const char *variable);
  void loadNetcdf(sspt_TMatrix< double > *v, const char *file, const char *variable);
#endif

};






template<typename T>
int get_var(int ncid, int varid, T *tp)
{
  return NC_EBADID;
}


template<typename T>
int put_var(int ncid, int varid, const T *tp)
{
  return NC_EBADID;
}



template<typename T>
int get_vara(int ncid, int varid, const size_t start[], const size_t count[],  T *tp)
{
  return NC_EBADID;
}


template<typename T>
int put_vara(int ncid, int varid, const size_t start[], const size_t count[],  const T *tp)
{
  return NC_EBADID;
}

//Declarations
template<> int get_var<char>(int ncid, int varid, char *tp);
template<> int get_var<unsigned char>(int ncid, int varid, unsigned char *tp);
template<> int get_var<signed char>(int ncid, int varid, signed char *tp);
template<> int get_var<int>(int ncid, int varid, int *tp);
template<> int get_var<double>(int ncid, int varid, double *tp);
template<> int get_var<float>(int ncid, int varid, float *tp);
template<> int put_var<char>(int ncid, int varid, const char *tp);
template<> int put_var<unsigned char>(int ncid, int varid, const unsigned char *tp);
template<> int put_var<signed char>(int ncid, int varid, const signed char *tp);
template<> int put_var<int>(int ncid, int varid, const int *tp);
template<> int put_var<double>(int ncid, int varid, const double *tp);
template<> int put_var<float>(int ncid, int varid, const float *tp);
template<> int get_vara<char>(int ncid, int varid, const size_t start[], const size_t count[],  char *tp);
template<> int get_vara<unsigned char>(int ncid, int varid, const size_t start[], const size_t count[],  unsigned char *tp);
template<> int get_vara<signed char>(int ncid, int varid, const size_t start[], const size_t count[],  signed char *tp);
template<> int get_vara<int>(int ncid, int varid, const size_t start[], const size_t count[],  int *tp);
template<> int get_vara<double>(int ncid, int varid, const size_t start[], const size_t count[],  double *tp);
template<> int get_vara<float>(int ncid, int varid, const size_t start[], const size_t count[],  float *tp);
template<> int put_vara<char>(int ncid, int varid, const size_t start[], const size_t count[],  const char *tp);
template<> int put_vara<unsigned char>(int ncid, int varid, const size_t start[], const size_t count[],  const unsigned char *tp);
template<> int put_vara<signed char>(int ncid, int varid, const size_t start[], const size_t count[],  const signed char *tp);
template<> int put_vara<int>(int ncid, int varid, const size_t start[], const size_t count[],  const int *tp);
template<> int put_vara<double>(int ncid, int varid, const size_t start[], const size_t count[],  const double *tp);
template<> int put_vara<float>(int ncid, int varid, const size_t start[], const size_t count[],  const float *tp);





#endif
