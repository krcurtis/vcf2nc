// Copyright 2017 Fred Hutchinson Cancer Research Center


#ifndef VCF40_H
#define VCF40_H

#include <string>
#include <map>
#include <vector>

#include "sspt_tmatrix.h"
#include "sspt_avltree.h"
#include "stringwrapper.h"



struct VCF40 {
  size_t nSNPs;
  size_t nSamples;

  std::multimap<std::string, std::string> headerPairs;


  //by snp
  std::vector<int> chromosome;
  std::vector<int> position;
  std::vector<std::string> snpName;  //maybe . if no known dbsnp mapping
  
  std::vector< std::string > referenceAllele;
  std::vector< std::vector<std::string> > alternateAllele;  //description of edits at the place
  std::vector< double > quality;
  std::vector< std::vector<std::string> > filters;
  std::vector< std::map<std::string, std::string> > info;
  std::vector< std::vector<std::string> > format;  //data types for coressponding sample columns
  
  //by sample
  std::vector< std::string > sampleID;

  //snps x samples
  //sspt_TMatrix< std::string > perSampleString;  // data
  sspt_TMatrix< const char * > perSampleString;  // data
  sspt_AVLTree< StringWrapper, const char *> uniqueStrings; // data store, potentially too slow, O(log2(N)) lookup time

  static bool loadVCF40(VCF40 *data, const char *file);


  std::vector<std::string> sampleGenotypeInfo(size_t i, size_t k);

};




#endif
