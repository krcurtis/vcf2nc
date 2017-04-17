// Copyright 2017 Fred Hutchinson Cancer Research Center

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>

#include "vcf40.h"

#include "sspt_delimiterparse.h"

#include "utilstext.h"
#include "utilsfilter.h"

#define MAX_INFO_FIELD_WIDTH 128

static bool parseKeyValue(std::string *key, std::string *value, const std::string line, bool allowEmpty=false)
{
  size_t i;
  for (i = 0; '=' != line[i] && i < line.size(); i++);

  if (i == line.size()) {
    if (!allowEmpty) {
      fprintf(stderr, "ERROR could not find '=' token\n" );
      return false;
    }
    *key = line.substr(0,i);
    *value = "1";
    return true;
  }
  *key = line.substr(0, i);
  *value = line.substr(i+1, line.size()-(i+1));

  return true;
}

static bool parseKeyValue(std::string *key, std::string *value, const char *inputString)
{
  std::string line(inputString);
  return parseKeyValue(key, value, line);
}



static int searchColumns(const char *name, const sspt_DelimiterParse &columns)
{
  for (size_t i = 0; i < columns.values(); i++) {
    if (0 == strcmp(name, columns.value(i))) {
        return i;
    }
  }
  return -1;
}




bool VCF40::loadVCF40(VCF40 *vcf, const char *file)
{
  FILE *fptr = fopen(file, "rb");
  if (0 == fptr) {
    fprintf(stderr, "ERROR cannot open %s\n", file);
    return false;
  }

  size_t width=0;
  size_t rows=0;
  UtilsText::scanLineWidth(&width, file);
  UtilsText::scanRowCount(&rows, file);

  printf("width %zu rows %zu\n", width, rows);
  width += 256; //safety buffer

  int chromosomeColumn = -1;
  int positionColumn = -1;
  int snpColumn = -1;
  int referenceColumn = -1;
  int alternateColumn = -1;
  int qualityColumn = -1;
  int filterColumn = -1;
  int infoColumn = -1;
  int formatColumn = -1;
  int sample0Column = -1;

  char *line = new char[width];
  size_t lineCount = 0;
  size_t refSNPLine = 0;



  while (0 != fgets(line, width, fptr)) {
    int lineLength = strlen(line);
    lineCount++;
    //printf("line %zu %i\n", lineCount, lineLength);

    if (line[ lineLength-1 ] == '\n')
      line[ lineLength-1 ] = 0;
    if (lineLength < 2) {
      fprintf(stderr, "ERROR line length too small at line %zu\n", lineCount);
      return false;
    }

    if ('#' == line[0] && '#' == line[1]) {  //parse key value pair
      std::string key, value;
      if (!parseKeyValue(&key, &value, line+2))
        return false;
      vcf->headerPairs.insert(std::pair<std::string, std::string>(key, value));
    }
    else if ('#' == line[0] ) { //parse column header
      vcf->nSNPs = rows - lineCount;
      refSNPLine = lineCount+1;
      

      vcf->chromosome.resize( vcf->nSNPs );
      vcf->position.resize( vcf->nSNPs );
      vcf->snpName.resize( vcf->nSNPs );
      vcf->referenceAllele.resize( vcf->nSNPs );
      vcf->alternateAllele.resize( vcf->nSNPs );
      vcf->quality.resize( vcf->nSNPs );
      vcf->filters.resize( vcf->nSNPs );
      vcf->info.resize( vcf->nSNPs );
      vcf->format.resize( vcf->nSNPs );

      sspt_DelimiterParse columns( line, '\t', false);

      
      chromosomeColumn = searchColumns("#CHROM", columns);
      positionColumn   = searchColumns("POS", columns);
      snpColumn        = searchColumns("ID", columns);
      referenceColumn  = searchColumns("REF", columns);
      alternateColumn  = searchColumns("ALT", columns);
      qualityColumn    = searchColumns("QUAL", columns);
      filterColumn     = searchColumns("FILTER", columns);
      infoColumn       = searchColumns("INFO", columns);
      formatColumn     = searchColumns("FORMAT", columns);
      sample0Column = formatColumn+1;
      if (-1 == formatColumn) {
        fprintf(stderr, "WARNING no genotype data columns\n");
        vcf->nSamples = 0;
      }
      else {  // yes there is a format column which implies there is per-sample info
        vcf->nSamples = columns.values() - sample0Column;
        vcf->sampleID.resize( vcf->nSamples );

        //vcf->sampleGenotypeInfo = sspt_TMatrix< std::vector<std::string> > (vcf->nSNPs,
        //                                                                    vcf->nSamples );
        //vcf->perSampleString = sspt_TMatrix< std::string > (vcf->nSNPs, vcf->nSamples );
        vcf->perSampleString = sspt_TMatrix< const char * > (vcf->nSNPs, vcf->nSamples );

        for (size_t i = 0; i < vcf->nSamples; i++) {
          vcf->sampleID[i] = columns.value(i + sample0Column);
        }
      }

    }
    else {  //parse column data
      size_t snpIndex = lineCount - refSNPLine;
      
      sspt_DelimiterParse columns( line, '\t', false);


      if (-1 != chromosomeColumn) {
        if ('c' == columns.value(chromosomeColumn)[0])  //skip over 'chr'
          vcf->chromosome[ snpIndex ] = ChromosomeFilter::filter( columns.value(chromosomeColumn) + 3 );
        else
          vcf->chromosome[ snpIndex ] = ChromosomeFilter::filter( columns.value(chromosomeColumn) );
      }

      if (-1 != positionColumn) {
        vcf->position[ snpIndex ] =  atoi( columns.value(positionColumn) );
      }

      if (-1 != snpColumn) {
        vcf->snpName[ snpIndex ] =  columns.value(snpColumn);
      }
      if (-1 != referenceColumn) {
        //vcf->referenceAllele[ snpIndex ] =  AlleleFilter::filter( columns.value(referenceColumn) );
        vcf->referenceAllele[ snpIndex ] =   columns.value(referenceColumn);
      }

      if (-1 != alternateColumn) {
        sspt_DelimiterParse fields( columns.value(alternateColumn), ',', false );
        std::vector<std::string> edits( fields.values() );
        for (size_t i = 0; i <  fields.values(); i++)
          edits[i] = fields.value(i);
        vcf->alternateAllele[ snpIndex ] =  edits;
      }

      if (-1 != qualityColumn) {
        vcf->quality[ snpIndex ] =  atof( columns.value(qualityColumn) );
      }

      if (-1 != filterColumn) {
        sspt_DelimiterParse fields( columns.value(filterColumn), ';', false);
        std::vector<std::string> filters( fields.values() );
        for (size_t i = 0; i <  fields.values(); i++)
          filters[i] = fields.value(i);
        vcf->filters[ snpIndex ] =  filters;
      }


      if (-1 != infoColumn) {
        sspt_DelimiterParse fields( columns.value(infoColumn), ';', false);
        std::map<std::string, std::string> pairs;
        for (size_t i = 0; i <  fields.values(); i++) {
          std::string key, value, group=fields.value(i);
          //from vcf 40 spec Keys without corresponding values are allowed in order to indicate group membership 
          //so just set to one
          if (!parseKeyValue(&key, &value, group, true))
            return false;
          //WORKAROUND large ANNO field size (greater that 2048) which causes problem in creation of netCDF
          if (value.size() > MAX_INFO_FIELD_WIDTH) {
            if (key != "ANNO")
              fprintf(stderr, "WARNING at SNP %zu, %s field is too big\n", snpIndex, key.c_str());
            value.resize(MAX_INFO_FIELD_WIDTH-1);
          }
          pairs.insert( std::pair<std::string, std::string>(key, value) );
        }
        vcf->info[ snpIndex ] =  pairs;
      }


      if (-1 != formatColumn) {
        sspt_DelimiterParse fields( columns.value(formatColumn), ':', false);
        std::vector<std::string> datatypes( fields.values() );
        for (size_t i = 0; i <  fields.values(); i++)
          datatypes[i] = fields.value(i);
        vcf->format[ snpIndex ] =  datatypes;
      }
#if 0
      //old way, parses only once, but the standard template library seems to be using a lot of memory, ~277 bytes per-sample, per-snp
      for (size_t i = 0; i < vcf->nSamples; i++) {
        sspt_DelimiterParse fields( columns.value(i + sample0Column), ':', false);
        std::vector<std::string> datafields( fields.values() );
        for (int k = 0; k <  fields.values(); k++)
          datafields[k] = fields.value(k);
        vcf->sampleGenotypeInfo(snpIndex, i) =  datafields;
      }
#endif

      //supposedly for short strings std::string is not so good.
      //http://jovislab.com/blog/?p=76

#if 0
      //try using strdup instead
      for (size_t i = 0; i < vcf->nSamples; i++) {
        //vcf->perSampleString(snpIndex, i) = columns.value(i + sample0Column);
        vcf->perSampleString(snpIndex, i) = strdup( columns.value(i + sample0Column) );
      }
#endif

#if 1
      //try finding unique strings, should be plenty based on previous experimenting
      for (size_t i = 0; i < vcf->nSamples; i++) {
        StringWrapper key( columns.value(i + sample0Column) );

        const char *unique = 0;
        if (vcf->uniqueStrings.find(key, &unique)) {
          //printf("dup: %s\n", unique);
        }
        else {
          unique = strdup( columns.value(i + sample0Column)  );
          StringWrapper insertKey(unique);

          //printf("before insert: %s\n", unique);
          if (!vcf->uniqueStrings.insert(insertKey, unique)) {
            fprintf(stderr, "ERROR could not insert %s\n", unique);
            return false;
          }
          //printf("after insert : %s\n", unique);
        }

        vcf->perSampleString(snpIndex, i) = unique;
      }
#endif

    }

  }
  fclose(fptr);



  return true;
}


std::vector<std::string> VCF40::sampleGenotypeInfo(size_t i, size_t k)
{
  //sspt_DelimiterParse fields( perSampleString(i,k).c_str(), ':', false);
  //printf("retrieve %s\n", perSampleString(i,k));
  sspt_DelimiterParse fields( perSampleString(i,k), ':', false);

  std::vector<std::string> datafields( fields.values() );
  for (int k = 0; k <  fields.values(); k++)
    datafields[k] = fields.value(k);
  return datafields;
}
