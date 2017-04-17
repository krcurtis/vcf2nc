// Copyright 2017 Fred Hutchinson Cancer Research Center

#include <stdio.h>
#include <stdlib.h>

#include "sspt_ascription.h"
#include "vcf40field-translator.h"

int main(int argc, char *argv[]) 
{
  sspt_Ascription options;

  const char *inputFile;
  const char *outputFile;
  const char *alternateHeaderFile=0;
  bool sort = true;
  bool duplicates = false;

  options.quality("i", &inputFile, true, "input file names");
  options.quality("o", &outputFile, true, "output file pathname");
  options.quality("alt", &alternateHeaderFile, false, "alternate header file (if the original vcf has errors)");
  options.quality("s", &sort, false, "<on|off> sort by chromosome,position");
  //options.quality("dup", &duplicates, false, "<on|off> allow duplicate positions when sorting");

  if (!options.evaluate(argc, argv)) {
    return -1;
  }


  //todo if vcf33
#if 0
  VCF33 *vcf = new VCF33;
  if (!VCF33::loadVCF33(vcf, inputFile)) {
    fprintf(stderr, "ERROR could not load %s\n", inputFile);
    return -1;
  }

  printf("VCF snps %zu samples %zu\n", vcf->nSNPs, vcf->nSamples);
  printf("VCF key-value pairs %zu\n", vcf->headerPairs.size() );

  VCF33Translator vt;
  if (!vt.process(outputFile, vcf, sort)) {
    fprintf(stderr, "ERROR could not convert vcf info %s into netCDF\n", inputFile);
    return -1;
  }
#endif
 VCF40 *vcf = new VCF40;
  if (!VCF40::loadVCF40(vcf, inputFile)) {
    fprintf(stderr, "ERROR could not load %s\n", inputFile);
    return -1;
  }

  VCF40 *alt = 0;
  if (0 != alternateHeaderFile) {
    alt = new VCF40;
    if (!VCF40::loadVCF40(alt, alternateHeaderFile)) {
      fprintf(stderr, "ERROR could not load %s\n", alternateHeaderFile);
      return -1;
    }
  }

  printf("VCF snps %zu samples %zu\n", vcf->nSNPs, vcf->nSamples);
  printf("VCF key-value pairs %zu\n", vcf->headerPairs.size() );
  if (0 != alt)
    printf("VCF (alternate) key-value pairs %zu\n", alt->headerPairs.size() );


  //VCF40Translator vt;
  VCF40FieldTranslator vt;
  if (!vt.process(outputFile, vcf, alt, sort)) {
    fprintf(stderr, "ERROR could not convert vcf info %s into netCDF\n", inputFile);
    return -1;
  }
  return 0;
}
