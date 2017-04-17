// Copyright 2017 Fred Hutchinson Cancer Research Center


#ifndef VCF_VARIABLE_H
#define VCF_VARIABLE_H

#include "netcdf.h"

#include "sspt_cord.h"
#include "sspt_list.h"

// one idea is that this is a convient way of storing information about netcdf stuff to create
// plus convienent way of extract item from vcf data type

//assume at most one alternate sequence ...


#define VCF_SAMPLE_DIM "Samples"
#define VCF_SNP_DIM  "SNPs"
#define VCF_STRING_DIM  "string_position"


class DataSetDescription;
class VCF40;

class VCFVariable {
 public:
  enum VCFType {
    VCF_NULL,
    VCF_GENO,
    VCF_AB,
    VCF_CHAR,
    VCF_CHROMOSOME,
    VCF_BYTE,
    VCF_INT,
    VCF_DOUBLE,
    VCF_STRING,
    VCF_QCSTATUS,
    VCF_GENDER,
  };


  //may lead to creating multidimensional arrays with dimension name 'arb4', and similar
  //VCFVariable(enum VCFColumn column, const char *label,  nc_type xtype, int number);
  virtual ~VCFVariable() { } 
  virtual bool updateDescription( DataSetDescription *desc )=0;
  virtual bool populateNetCDF(int ncid,  VCF40 *vcf)=0;
  virtual void variableName(sspt_Cord *name)=0;

 private:

  //todo
};



//per-snp
class VCFVariableColumnInfo : public VCFVariable {
 public:
  //may lead to creating multidimensional arrays with dimension name 'arb4', and similar
  VCFVariableColumnInfo(const char *label,  enum VCFType vcftype, int number);
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;
  sspt_Cord m_field;
  enum VCFType m_vcftype;
  int m_number;
  size_t m_factor;

  //bool storeInteger(int ncid,  VCF40 *vcf);

};



//per-sample, per-snp
class VCFVariableColumnFormat : public VCFVariable {
 public:
  //may lead to creating multidimensional arrays with dimension name 'arb4', and similar
  VCFVariableColumnFormat(const char *label,  enum VCFType vcftype, int number);
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:

  sspt_Cord m_varname;
  sspt_Cord m_field;
  enum VCFType m_vcftype;
  int m_number;
  size_t m_factor;


};



//Special version for handling GT field like  '0/1' into two alleles/variables 
//per-sample, per-snp

class VCFVariableGenotype : public VCFVariable {
 public:
  //may lead to creating multidimensional arrays with dimension name 'arb4', and similar
  VCFVariableGenotype(const char *label); //,  enum VCFType vcftype, int number);
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;
  sspt_Cord m_field;
  enum VCFType m_vcftype;
  int m_number;
  size_t m_factor;


  bool storeAB(int ncid,  VCF40 *vcf);
};


//Special version to handle chromosome and position

class VCFVariableLocation : public VCFVariable {
 public:
  VCFVariableLocation();
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;


  bool storeChromosome(int ncid,  VCF40 *vcf);
  bool storePosition(int ncid,  VCF40 *vcf);
};


class VCFVariableID : public VCFVariable {
 public:
  VCFVariableID();
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;
  size_t m_stringWidth;

  bool store(int ncid,  VCF40 *vcf);
};


class VCFVariableQuality : public VCFVariable {
 public:
  VCFVariableQuality();
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;

  bool storeQuality(int ncid,  VCF40 *vcf);
};


class VCFVariableAutoFilter : public VCFVariable {
 public:
  VCFVariableAutoFilter(VCF40 *vcf);
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;
  sspt_List<sspt_Cord> m_filters;

  bool storeFlag(int ncid,  VCF40 *vcf,  const char *variableName, const sspt_Cord &filter);
};

class VCFVariableSimpleFilter : public VCFVariable {
 public:
  VCFVariableSimpleFilter();
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;
  size_t m_stringWidth;

  bool store(int ncid,  VCF40 *vcf);
};



// for ref allele
class VCFVariableRefAllele : public VCFVariable {
 public:
  VCFVariableRefAllele();
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;


  bool storeAllele(int ncid,  VCF40 *vcf);

};


// for alternate allele
class VCFVariableAltAllele : public VCFVariable {
 public:
  VCFVariableAltAllele(int whichAlternate);
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  size_t m_whichAlternate;
  sspt_Cord m_varname;

  bool storeAllele(int ncid,  VCF40 *vcf, size_t column);
};



//per-snp
class VCFVariablePlaceHolder : public VCFVariable {
 public:
  //may lead to creating multidimensional arrays with dimension name 'arb4', and similar
  VCFVariablePlaceHolder(const char *label,  enum VCFType vcftype, bool perSNP=true, bool perSample=false, bool perString=true);
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  enum DimCode {
    DIM_NONE = 0x0,
    DIM_STRING = 0x1,
    DIM_SAMPLE = 0x2,
    DIM_SAMPLE_STRING = 0x3,
    DIM_SNP    = 0x4,
    DIM_SNP_STRING = 0x5,
    DIM_SNP_SAMPLE = 0x6,
    DIM_ALL = 0x7
  };


  sspt_Cord m_varname;
  enum VCFType m_vcftype;
  enum DimCode m_dimCode;

  //bool storeInteger(int ncid,  VCF40 *vcf);

};




// for ref allele
class VCFVariableSample : public VCFVariable {
 public:
  VCFVariableSample();
  bool updateDescription( DataSetDescription *desc );
  bool populateNetCDF(int ncid,  VCF40 *vcf);
  void variableName(sspt_Cord *name) { *name = m_varname; }

 private:
  sspt_Cord m_varname;
  size_t m_stringWidth;
  enum VCFType m_vcftype;


  bool storeString(int ncid,  VCF40 *vcf);

};




#endif
