// Copyright 2017 Fred Hutchinson Cancer Research Center




#ifndef DATASETDESCRIPTION_H
#define DATASETDESCRIPTION_H

#include "sspt_cord.h"
#include "sspt_avltree.h"
#include "utilsnetcdf.h"


struct DimensionDesc;
struct VariableDesc;


class DataSetDescription {
public:
  DataSetDescription();
  ~DataSetDescription();


  //idempotent in terms of adding, if you need to change the size, use reviseDimensionSize
  bool addDimension(const char *dimname, size_t size);
  //add a variable with known dimensions
  bool addVariable(const char *varname, nc_type xtype, const char *dim1, const char *dim2=0, const char *dim3=0);


  bool reviseDimensionSize(const char *dimension, size_t revisedSize);
  bool createEmptyNetCDF(int *ncid, const char *file);

  bool dimensionSize(const char *dimension, size_t *size);


private:
  sspt_AVLTree<sspt_Cord, VariableDesc*> m_vars;
  sspt_AVLTree<sspt_Cord, DimensionDesc*> m_dims;

};

#endif 


