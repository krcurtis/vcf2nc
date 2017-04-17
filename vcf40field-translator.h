// Copyright 2017 Fred Hutchinson Cancer Research Center


#ifndef VCF40FIELD_TRANSLATOR_H
#define VCF40FIELD_TRANSLATOR_H

#include "vcf40.h"

#include <stdio.h>
#include <stdlib.h>

#include "sspt_hashtable.h"
#include "sspt_delimiterparse.h"
#include "sspt_cord.h"
#include "sspt_array.h"
#include "sspt_avltree.h"

#include "utilsnetcdf.h"
#include "netcdf.h"

#include "vcfvariable.h"


class VCF40FieldTranslator {
 public:

  VCF40FieldTranslator();
  ~VCF40FieldTranslator();

  bool process(const char *outputFile, 
               VCF40 *vcf, 
               VCF40 *alternateHeader,
               bool sortSNPs);
  

  void autofilter(bool flag) { m_autofilter = flag; }

 private:

  sspt_HashTable<sspt_Cord, VCFVariable*> m_variableTable;


  size_t m_nSamples;
  size_t m_nSNPs;
  int m_ncid;


  char *m_buffer;
  unsigned long m_bufferSize;

  bool m_autofilter;  //if true, expand filter column to boolean vectors
  //bool m_allowDuplicates;

  bool extractVariableInfo(std::string *label, std::string *vtype, std::string *number, const char *item);


  bool addInfoVar(const std::string &label, const std::string &vtype, const std::string &number);
  bool addFormatVar(const std::string &label, const std::string &vtype, const std::string &number);


  void selectVariables(VCF40 *vcf);
  bool createDescription(DataSetDescription **output, VCF40 *vcf);
  bool createNetcdfVariables(const char *outputFile, DataSetDescription *desc);


  bool sortByChromosome( sspt_Array<int> *mapping );
  bool parseGenotype(std::string *first, std::string *second, const std::string line);
  bool processVCF(VCF40 *vcf, bool sortSNPs);


  //  bool processR2(const char *file);
  //  bool indexSNPs(const char *netcdfFile);




};


#endif
