// Copyright 2017 Fred Hutchinson Cancer Research Center

#include <stdio.h>
#include <stdlib.h>

#include "sspt_ascription.h"
#include "vcf40.h"


int main(int argc, char *argv[]) 
{


  sspt_Ascription options;

  const char *inputFile;
  //char *outputFile;
  //bool sort = false;
  //bool duplicates = false;

  bool genotypeFlag = false;
  bool snpFlag = true;
  bool sampleFlag = false;
  bool filterFlag = false;

  options.quality("i", &inputFile, true, "input file name");
  //options.quality("o", &outputFile, true, "output file pathname");
  //options.quality("s", &sort, false, "<on|off> sort by chromosome,position");
  //options.quality("dup", &duplicates, false, "<on|off> allow duplicate positions when sorting");

  if (!options.evaluate(argc, argv)) {
    return -1;
  }

  VCF40 *vcf = new VCF40;
  if (!VCF40::loadVCF40(vcf, inputFile)) {
    fprintf(stderr, "ERROR could not load %s\n", inputFile);
    return -1;
  }

  printf("VCF snps %zu samples %zu\n", vcf->nSNPs, vcf->nSamples);
  printf("VCF key-value pairs %zu\n", vcf->headerPairs.size() );



  if (genotypeFlag) {
    for (size_t i = 0; i < vcf->nSNPs; i++) {
      std::vector<std::string> format = vcf->format[i];  //data types for coressponding sample columns

      int index = -1;
      for (size_t j = 0; j < format.size(); j++) {
        if (format[j] == "GT") {
          index = j;
          break;
        }
      }
      if (-1 == index) {
        fprintf(stderr, "ERROR no genotype format found at %zu\n", i);
        return -1;
      }

      for (size_t k=0; k < vcf->nSamples; k++) {
        std::vector<std::string> fields = vcf->sampleGenotypeInfo(i,k);
        printf("  %s", fields[index].c_str());
      }
      printf("\n");
    }
  }

  if (snpFlag) {
    for (size_t i = 0; i < vcf->nSNPs; i++) {
      printf("%s %i %i %s %zu %s %lf\n", 
             vcf->snpName[i].c_str(), 
             vcf->chromosome[i],
             vcf->position[i],
             vcf->referenceAllele[i].c_str(),
             vcf->alternateAllele[i].size(),
             vcf->alternateAllele[i][0].c_str(),
             vcf->quality[i]);
    }
  }

  if (sampleFlag) {
    for (size_t i = 0; i < vcf->nSamples; i++) {
      printf("%s\n", 
             vcf->sampleID[i].c_str());
    }

  }

  if (filterFlag) {

    for (size_t i = 0; i < vcf->nSNPs; i++) {
      printf("%s %zu", 
             vcf->snpName[i].c_str(), 
             vcf->filters[i].size());
      for (size_t k = 0; k < vcf->filters[i].size(); k++) {
        printf(" %s", vcf->filters[i][k].c_str());
      }
      printf("\n");
    }


  }

  return 0;
}
