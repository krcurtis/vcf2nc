// Copyright 2017 Fred Hutchinson Cancer Research Center

#include <string.h>

#include "datasetdescription.h"
#include "variable.h"

#include "sspt_list.h"

struct DimensionDesc {
  char name[NC_MAX_NAME+1];
  size_t size;
  bool unlimited;
};


struct VariableDesc {
  nc_type xtype;
  char name[NC_MAX_NAME+1];
  sspt_Array<DimensionDesc*> dims;
};



DataSetDescription::DataSetDescription()
{
}


DataSetDescription::~DataSetDescription()
{
  //TODO free memory
}


bool DataSetDescription::addDimension(const char *dimname, size_t size)
{
  sspt_Cord dname(dimname);
  DimensionDesc *desc = 0;
  if (!m_dims.find(dname, &desc)) {
    desc = new DimensionDesc;
    strncpy(desc->name, dname.c_str(), NC_MAX_NAME);
    desc->unlimited = false;
    desc->size = size;
    m_dims.insert(dname, desc);
  }
  return true;
}


bool DataSetDescription::addVariable(const char *varname, nc_type xtype, const char *dim1, const char *dim2, const char *dim3)
{
  sspt_Cord key(varname);


  if (m_vars.find(key, 0)) {
    fprintf(stderr, "ERROR duplicate variable when creating description %s\n", varname);
    return false;
  }

  //todo switch to using variable arguments
  sspt_List<const char *> dimNamesList;
  dimNamesList.insertRear(dim1);
  dimNamesList.insertRear(dim2);
  dimNamesList.insertRear(dim3);

  sspt_Array<const char *> dimNames;
  dimNamesList.toArray(&dimNames);

  sspt_List<DimensionDesc*> dimList;

  for (size_t i = 0; i < dimNames.size(); i++) {
    const char *name = dimNames[i];
    if (0 == name)
      continue;

    sspt_Cord dname(name);
    DimensionDesc *desc = 0;
    if (!m_dims.find(dname, &desc)) {
      fprintf(stderr, "ERROR adding variable %s with unknown dimension %s\n", varname, name);
      return false;
    }
    dimList.insertRear(desc);
  }



  VariableDesc *vdesc = new VariableDesc;
  vdesc->xtype = xtype;
  strncpy(vdesc->name, varname, NC_MAX_NAME+1);
  dimList.toArray(&vdesc->dims);

  m_vars.insert(key, vdesc);
  return true;
}



bool DataSetDescription::reviseDimensionSize(const char *name, size_t revisedSize)
{
  sspt_Cord dname(name);

  DimensionDesc *desc = 0;
  if (!m_dims.find(dname, &desc)) {
    fprintf(stderr, "ERROR could not find dimension %s\n", name);
    return false;
  }
  desc->size = revisedSize;
  return true;
}




bool DataSetDescription::createEmptyNetCDF(int *ncid, const char *file)
{
  sspt_AVLTree<sspt_Cord,int> dimensionIndex;

  int nret;

  //open netcdf file
  nret = nc_create(file, NC_CLOBBER  | NC_NETCDF4, ncid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR failed to create netcdf output file %s", file);

  //two passes, first unlimited, second normal
  for (sspt_AVLIterator<sspt_Cord,DimensionDesc*> iter = m_dims.begin();
       !iter.atEnd(); iter.moveNext()) {
    DimensionDesc *d = iter.data();
    if (!d->unlimited)
      continue;
    int dim;
    nret = nc_def_dim(*ncid, d->name, NC_UNLIMITED, &dim);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR failed to create %s dimension", d->name);
    dimensionIndex.insert(iter.key(), dim);
  }
  for (sspt_AVLIterator<sspt_Cord,DimensionDesc*> iter = m_dims.begin();
       !iter.atEnd(); iter.moveNext()) {
    DimensionDesc *d = iter.data();
    if (d->unlimited)
      continue;
    int dim;
    nret = nc_def_dim(*ncid, d->name, d->size, &dim);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR failed to create %s dimension", d->name);
    dimensionIndex.insert(iter.key(), dim);
  }



  for (sspt_AVLIterator<sspt_Cord,VariableDesc*> iter = m_vars.begin();
       !iter.atEnd(); iter.moveNext()) {
    VariableDesc *v = iter.data();

    int *dims = new int[v->dims.size()];
    int var;

    for (size_t i = 0; i < v->dims.size(); i++) {
      DimensionDesc *d = v->dims[i];
      sspt_Cord name(d->name);
      if (!dimensionIndex.find(name, &dims[i])) {
        fprintf(stderr, "ERROR could not find dimension %s\n", d->name);
        return false;
      }
    }

    nret = nc_def_var(*ncid, v->name, v->xtype, v->dims.size(), dims, &var);
    FALSE_ON_NETCDF_ERROR(nret, "ERROR failed to create %s variable\n", v->name);
    delete[] dims;
  }

  ncendef(*ncid);
  //nc_close(*ncid);

  return true;
}



bool DataSetDescription::dimensionSize(const char *dimension, size_t *size)
{
  sspt_Cord dname(dimension);

  DimensionDesc *desc = 0;
  if (!m_dims.find(dname, &desc)) {
    fprintf(stderr, "ERROR could not find dimension %s\n", dimension);
    return false;
  }
  *size = desc->size;
  return true;

}
