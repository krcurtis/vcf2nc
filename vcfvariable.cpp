// Copyright 2017 Fred Hutchinson Cancer Research Center

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "sspt_delimiterparse.h"
#include "sspt_avltree.h"
#include "sspt_symbolstore.h"

#include "vcfvariable.h"

#include "datasetdescription.h"
#include "vcf40.h"

#include "stringtranslator.h"



#include "vcf_names.h"


static nc_type mapVCFType(enum VCFVariable::VCFType vcftype)
{
  switch (vcftype) {
  case VCFVariable::VCF_GENO:
    return NC_BYTE;

  case VCFVariable::VCF_AB:
    return NC_BYTE;

  case VCFVariable::VCF_CHAR:
    return NC_CHAR;

  case VCFVariable::VCF_CHROMOSOME:
    return NC_BYTE;

  case VCFVariable::VCF_BYTE:
    return NC_BYTE;

  case VCFVariable::VCF_INT:
    return NC_INT;

  case VCFVariable::VCF_DOUBLE:
    return NC_DOUBLE;

  case VCFVariable::VCF_STRING:
    return NC_CHAR;

  case VCFVariable::VCF_QCSTATUS:
    return NC_BYTE;
    break;

  case VCFVariable::VCF_GENDER:
    return NC_BYTE;

  default:
    break;
  }
  fprintf(stderr, "ERROR invalid internal mapping of VCFType enums to netcdf type %i\n", vcftype);
  exit(-1);
}


static inline bool getColumnField(std::string *out, VCF40 *vcf, size_t index, const char *field)
{
  std::map<std::string, std::string> info = vcf->info[index];
  *out = info[ field ];
  return true;
}



template <typename T>
bool storeColumn(int ncid,  VCF40 *vcf, size_t factor, const char *varname, const char *field, StringTranslator<T> *translator)
{
  int nret;
  int varid;
  size_t N = vcf->nSNPs * factor;

  // will know that it uses integer filter ...
  //build integer filter for specific m_vcftype ?


  nret = nc_inq_varid(ncid, varname, &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", varname);



  T *buffer = new T[N];
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    std::string s;
    if (!getColumnField(&s, vcf, i, field) ) {
      fprintf(stderr, "ERROR could not find %s at snp index %zu\n", varname, i);
      return false;
    }

    if (factor > 1) { //parse ...
      sspt_DelimiterParse p(s.c_str(), ',', false);
      size_t found = p.values();
      if (found > factor) {
        fprintf(stderr, "ERROR expected %zu values, found %zu at snp index %zu\n", factor, p.values(), i);
        return false;
      }
      for (size_t k = 0; k < factor; k++) {
        if (k < found)
          buffer[factor*i+k] = translator->translate( p.value(k) );
        else
          buffer[factor*i+k] = translator->translate( "" );
      }
    }
    else {
      buffer[i] = translator->translate(s.c_str());
    }
  }
  //end for loop

  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", varname );

  return true;
}





template <typename T>
bool storeMatrix(int ncid,  VCF40 *vcf, size_t factor, const char *varname, const char *field, StringTranslator<T> *translator)
{
  int nret;
  int varid;
  size_t N = vcf->nSamples * vcf->nSNPs * factor;

  // will know that it uses integer filter ...
  //build integer filter for specific m_vcftype ?


  nret = nc_inq_varid(ncid, varname, &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", varname);



  T *buffer = new T[N];
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    std::vector<std::string> formatIDs = vcf->format[i];

    //sspt_DelimiterParse formatParse(formatIDs.c_str(), ':', false);
    int subfield = -1;
    for (size_t j = 0; j < formatIDs.size() && subfield < 0; j++) {
      //printf("-%s\n", formatIDs[j].c_str());
      if (formatIDs[j] == field)  {
        subfield = j;
      }
    }
    if (subfield == -1) {
      fprintf(stderr, "WARNING could not find field %s in format at snp index %zu\n", field, i);
      //return true;
      for (size_t k = 0; k < vcf->nSamples; k++) {
        for (size_t m = 0; m < factor; m++) {
          buffer[ k*(vcf->nSNPs * factor) + i*factor + m] = -1;
        }
      }
      continue;
    }

    for (size_t k = 0; k < vcf->nSamples; k++) {
      std::vector<std::string> fields = vcf->sampleGenotypeInfo(i, k);

      //check for special case empty data, v3.3 apparently used empty string...
      if (1 == fields.size() && 
          (fields[0] == "./." || fields[0] == "") ) {
        for (size_t m = 0; m < factor; m++) {
          buffer[ k*(vcf->nSNPs * factor) + i*factor + m] = -1;
        }
        continue;
      }
      


      std::string entry = fields[subfield];

      if (factor > 1) { //parse ...
        sspt_DelimiterParse p(entry.c_str(), ',', false);
        size_t found = p.values();
        if (p.values() > factor) {
          fprintf(stderr, "ERROR (in %s) expected %zu values, found %zu at snp index %zu\n", __FUNCTION__,factor, p.values(), i);
          return false;
        }
        for (size_t m = 0; m < factor; m++) {
          if (m < found)
            buffer[ k*(vcf->nSNPs * factor) + i*factor + m] = translator->translate( p.value(m) );
          else
            buffer[ k*(vcf->nSNPs * factor) + i*factor + m] = translator->translate( "" );
        }
      }
      else {
        buffer[ k*(vcf->nSNPs * factor) + i*factor  ] = translator->translate( entry.c_str() );
      }

    }     //end sample loop
  } // end snp loop

  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", varname );

  return true;
}






VCFVariableColumnInfo::VCFVariableColumnInfo(const char *label,  enum VCFType vcftype, int number)
{
  m_varname = "info_";
  m_varname.append(label);
  m_field = label;


  m_vcftype = vcftype;
  m_number = number;
  m_factor = 0;
  if (m_number == 0 || m_number == 1) {
    m_factor = 1;
  }
  else if (m_number > 1) {
    m_factor = m_number;
  }
}


bool VCFVariableColumnInfo::updateDescription( DataSetDescription *desc )
{
  if (m_number == 0 || m_number == 1) {
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype), VCF_SNP_DIM);
  }
  else if (m_number > 1) {
    char dim2[128];
    snprintf(dim2, 128, "arb%i", m_number);
    desc->addDimension(dim2, m_number);
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype), VCF_SNP_DIM, dim2);
  }
  return false;
}


bool VCFVariableColumnInfo::populateNetCDF(int ncid,  VCF40 *vcf)
{
  nc_type xtype = mapVCFType(m_vcftype);

  switch (xtype) {
  case NC_INT: {
    PlainTranslator<int> translator;
    return storeColumn(ncid, vcf, m_factor, m_varname.c_str(), m_field.c_str(), &translator);
  }
    break;

  case NC_DOUBLE: {
    PlainTranslator<double> translator;
    return storeColumn(ncid, vcf, m_factor, m_varname.c_str(), m_field.c_str(), &translator);
  }
    break;

  default:
    return false;
  };

  return true;
}






VCFVariableColumnFormat::VCFVariableColumnFormat(const char *label,  enum VCFType vcftype, int number)
{
  m_varname = "array_";
  m_varname.append(label);
  m_field = label;


  m_vcftype = vcftype;
  m_number = number;
  m_factor = 0;
  if (m_number == 0 || m_number == 1) {
    m_factor = 1;
  }
  else if (m_number > 1) {
    m_factor = m_number;
  }
}


bool VCFVariableColumnFormat::updateDescription( DataSetDescription *desc )
{
  if (m_number == 0 || m_number == 1) {
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype), VCF_SAMPLE_DIM, VCF_SNP_DIM);
  }
  else if (m_number > 1) {
    char dim2[128];
    snprintf(dim2, 128, "arb%i", m_number);
    desc->addDimension(dim2, m_number);
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype), VCF_SAMPLE_DIM, VCF_SNP_DIM, dim2);
  }
  return false;
}


bool VCFVariableColumnFormat::populateNetCDF(int ncid,  VCF40 *vcf)
{
  nc_type xtype = mapVCFType(m_vcftype);

  switch (xtype) {
  case NC_INT: {
    PlainTranslator<int> translator;
    return storeMatrix(ncid, vcf, m_factor, m_varname.c_str(), m_field.c_str(), &translator);
  }
    break;

  case NC_DOUBLE: {
    PlainTranslator<double> translator;
    return storeMatrix(ncid, vcf, m_factor, m_varname.c_str(), m_field.c_str(), &translator);
  }
    break;

  default:
    return false;
  };

  return true;
}




VCFVariableGenotype::VCFVariableGenotype(const char *label)
{
  m_varname = "array_";
  m_varname.append(label);
  m_field = label;

  m_vcftype = VCF_AB;
  m_number = 3;
  m_factor = 3;
}


bool VCFVariableGenotype::updateDescription( DataSetDescription *desc )
{
  char dim2[128];
  snprintf(dim2, 128, "arb%i", m_number);
  desc->addDimension(dim2, m_number);
  return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype), VCF_SAMPLE_DIM, VCF_SNP_DIM, dim2);
}



bool VCFVariableGenotype::populateNetCDF(int ncid,  VCF40 *vcf)
{
  return storeAB(ncid, vcf);
}





bool VCFVariableGenotype::storeAB(int ncid,  VCF40 *vcf)
{
  GTTranslator<signed char> translator;

  int nret;
  int varid;
  size_t N = vcf->nSamples * vcf->nSNPs * m_factor;

  // will know that it uses integer filter ...
  //build integer filter for specific m_vcftype ?


  nret = nc_inq_varid(ncid, m_varname.c_str(), &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", m_varname.c_str());



  signed char *buffer = new signed char[N];
  size_t shortFieldCount = 0;
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    std::vector<std::string> formatIDs = vcf->format[i];

    //sspt_DelimiterParse formatParse(formatIDs.c_str(), ':', false);
    int subfield = -1;
    for (size_t j = 0; j < formatIDs.size() && subfield < 0; j++) {
      //printf("-%s\n", formatIDs[j].c_str());
      if (formatIDs[j] == m_field.c_str())  {
        subfield = j;
      }
    }
    if (subfield == -1) {
      fprintf(stderr, "ERROR could not find field %s in format at snp index %zu\n", m_field.c_str(), i);
      return false;
    }

    for (size_t k = 0; k < vcf->nSamples; k++) {
      std::vector<std::string> fields = vcf->sampleGenotypeInfo(i, k);
      std::string entry = fields[subfield];
      const char *string = entry.c_str();


      if (0 == strlen(string)) {
        //fprintf(stderr, "ERROR field %s is empty\n", m_field.c_str());
        //return false;
        //glfMultiples VCF file output seems to produce empty strings
        shortFieldCount++;
        const char *slash = "/";
        const char *dot = ".";
        buffer[ k*(vcf->nSNPs * m_factor) + i*m_factor + 0] = translator.translate( dot );
        buffer[ k*(vcf->nSNPs * m_factor) + i*m_factor + 1] = translator.translate( slash );
        buffer[ k*(vcf->nSNPs * m_factor) + i*m_factor + 2] = translator.translate( dot );
      }
      else if (strlen(string) < m_factor) {  // this seems common in vcf files 
        shortFieldCount++;
        const char *r1 = "/";
        const char *r2 = ".";

        buffer[ k*(vcf->nSNPs * m_factor) + i*m_factor + 0] = translator.translate( string );
        buffer[ k*(vcf->nSNPs * m_factor) + i*m_factor + 1] = translator.translate( r1 );
        buffer[ k*(vcf->nSNPs * m_factor) + i*m_factor + 2] = translator.translate( r2 );
        /* TODO put in back in check for 0,1,2,3 ?
           if ('0' == string[0])
           else if ('1' == string[0]) {
           else if ('2' == string[0]) {
           else if ('3' == string[0]) {
           fprintf(stderr, "ERROR field %s entry is too short, %s \n", m_field.c_str(), string);
           fprintf(stderr, "cont'd: i %zu k %zu %s\n", i, k, string);
           return false;
           }
        */
      }
      else { //normal expected case
        for (size_t m = 0; m < m_factor; m++) {
          buffer[ k*(vcf->nSNPs * m_factor) + i*m_factor + m] = translator.translate( string+m );
        }
      }

    }     //end sample loop
  } // end snp loop

  nret = put_var(ncid, varid, buffer);
  delete[] buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", m_varname.c_str() );

  if (shortFieldCount > 0)
    printf("Short field count %zu\n", shortFieldCount);

  return true;

}




VCFVariableLocation::VCFVariableLocation()
{
  m_varname = "ChromosomePosition";
}

bool  VCFVariableLocation::updateDescription( DataSetDescription *desc )
{
  if (!desc->addVariable(CHROMOSOME, NC_BYTE, VCF_SNP_DIM))
    return false;
  return desc->addVariable(POSITION, NC_INT, VCF_SNP_DIM);
}


bool  VCFVariableLocation::populateNetCDF(int ncid,  VCF40 *vcf)
{
  return storeChromosome(ncid, vcf) && storePosition(ncid, vcf);
}


bool  VCFVariableLocation::storeChromosome(int ncid,  VCF40 *vcf)
{
  int nret;
  int varid;
  size_t N = vcf->nSNPs;

  // will know that it uses integer filter ...
  //build integer filter for specific m_vcftype ?


  nret = nc_inq_varid(ncid, CHROMOSOME, &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", CHROMOSOME);

  signed char *buffer = new signed char[N];
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    buffer[i] = vcf->chromosome[i];
  }
  //end for loop

  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", CHROMOSOME );
  return true;
}


bool  VCFVariableLocation::storePosition(int ncid,  VCF40 *vcf)
{
  int nret;
  int varid;
  size_t N = vcf->nSNPs;

  // will know that it uses integer filter ...
  //build integer filter for specific m_vcftype ?


  nret = nc_inq_varid(ncid, POSITION, &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", POSITION);

  int *buffer = new int[N];
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    buffer[i] = vcf->position[i];
  }
  //end for loop

  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", POSITION );
  return true;
}






VCFVariableID::VCFVariableID()
{
  m_varname = "ID";
}

bool  VCFVariableID::updateDescription( DataSetDescription *desc )
{
  return desc->dimensionSize(VCF_STRING_DIM, &m_stringWidth)
    && desc->addVariable(m_varname.c_str(), NC_CHAR,  VCF_SNP_DIM, VCF_STRING_DIM);
}


bool  VCFVariableID::populateNetCDF(int ncid,  VCF40 *vcf)
{
  return store(ncid, vcf);
}


bool  VCFVariableID::store(int ncid,  VCF40 *vcf)
{
  assert(m_stringWidth > 1);
  size_t N = vcf->nSNPs * m_stringWidth;


  int varid;
  int nret = nc_inq_varid(ncid, m_varname.c_str(), &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", m_varname.c_str());

  char *buffer = new char[N];
  memset(buffer, 0, N);

  for (size_t i = 0; i < vcf->nSNPs; i++) {
    const char *ptr = vcf->snpName[i].c_str();
    size_t n = strlen(ptr);
    n = (n <= (m_stringWidth-1)) ? n : m_stringWidth-1;
    memcpy(buffer + i*m_stringWidth, ptr, n);
  }
  
  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", m_varname.c_str() );
  return true;

}








VCFVariableQuality::VCFVariableQuality()
{
  m_varname = "QUAL";
}

bool  VCFVariableQuality::updateDescription( DataSetDescription *desc )
{
  return desc->addVariable(m_varname.c_str(), NC_DOUBLE, VCF_SNP_DIM);
}


bool  VCFVariableQuality::populateNetCDF(int ncid,  VCF40 *vcf)
{
  return storeQuality(ncid, vcf);
}


bool  VCFVariableQuality::storeQuality(int ncid,  VCF40 *vcf)
{
  int nret;
  int varid;
  size_t N = vcf->nSNPs;

  // will know that it uses integer filter ...
  //build integer filter for specific m_vcftype ?


  nret = nc_inq_varid(ncid, m_varname.c_str(), &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", m_varname.c_str());

  double *buffer = new double[N];
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    buffer[i] = vcf->quality[i];
  }
  //end for loop

  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", m_varname.c_str() );
  return true;
}






VCFVariableAutoFilter::VCFVariableAutoFilter(VCF40 *vcf)
{
  m_varname = "AutoFilter";

  int count = 0;
  sspt_AVLTree<sspt_Cord,int> table;
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    std::vector<std::string> filters = vcf->filters[i];
    for (size_t k = 0; k < filters.size(); k++) {
      std::string value = filters[k];
      sspt_Cord key( value.c_str() );
      if (key == ".") {
        //nothing
      }
      else if (table.find(key, &count)) {
        //nothing
      }
      else {
        m_filters.insertRear(key);


        table.insert(key, count);
      }
    } //end for loop over filters for one position
  } //end for loop over filters for multiple positions


}



bool  VCFVariableAutoFilter::updateDescription( DataSetDescription *desc )
{
  bool result = true;
  for (sspt_ListIterator<sspt_Cord> iter = m_filters.begin(); !iter.atEnd() && result; iter.moveNext()) {
    sspt_Cord filter = iter.current();
    sspt_Cord name("flag_");
    name.append(&filter);
    result = desc->addVariable(name.c_str(), NC_BYTE, VCF_SNP_DIM);
  }
  return result;
}


bool  VCFVariableAutoFilter::populateNetCDF(int ncid,  VCF40 *vcf)
{
  bool result = true;
  for (sspt_ListIterator<sspt_Cord> iter = m_filters.begin(); !iter.atEnd() && result; iter.moveNext()) {
    sspt_Cord filter = iter.current();
    sspt_Cord varname("flag_");
    varname.append(filter.c_str());
    result = storeFlag(ncid, vcf, varname.c_str(), filter);
  }
  return result;

}


bool  VCFVariableAutoFilter::storeFlag(int ncid,  VCF40 *vcf,  const char *variableName, const sspt_Cord &filter)
{
  int nret;
  int varid;
  size_t N = vcf->nSNPs;

  // will know that it uses integer filter ...
  //build integer filter for specific m_vcftype ?

  sspt_Cord test = filter;
  nret = nc_inq_varid(ncid, variableName, &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", variableName);

  signed char *buffer = new signed char[N];

  for (size_t i = 0; i < vcf->nSNPs; i++) {
    std::vector<std::string> filters = vcf->filters[i];
    signed char flag = 0;

    for (size_t k = 0; k < filters.size(); k++) {
      std::string value = filters[k];
      sspt_Cord key( value.c_str() );
      //printf("-%s-%s-\n", key.c_str(), test.c_str());
      if (key == filter) {
        flag = 1;
      }
    } //end for loop over filters for one position

    buffer[i] = flag;
  } //end for loop over filters for multiple positions


  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", variableName );
  return true;
}





VCFVariableSimpleFilter::VCFVariableSimpleFilter()
{
  m_varname = "FILTER";
}

bool  VCFVariableSimpleFilter::updateDescription( DataSetDescription *desc )
{
  return desc->dimensionSize(VCF_STRING_DIM, &m_stringWidth)
    && desc->addVariable(m_varname.c_str(), NC_CHAR,  VCF_SNP_DIM, VCF_STRING_DIM);
}


bool  VCFVariableSimpleFilter::populateNetCDF(int ncid,  VCF40 *vcf)
{
  return store(ncid, vcf);
}


bool  VCFVariableSimpleFilter::store(int ncid,  VCF40 *vcf)
{
  assert(m_stringWidth > 1);
  size_t N = vcf->nSNPs * m_stringWidth;


  int varid;
  int nret = nc_inq_varid(ncid, m_varname.c_str(), &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", m_varname.c_str());

  char *buffer = new char[N];
  memset(buffer, 0, N);

  for (size_t i = 0; i < vcf->nSNPs; i++) {
    std::vector<std::string> filters = vcf->filters[i];

    sspt_Cord rebuild("");

    
    for (size_t k = 0; k < filters.size(); k++) {
      std::string value = filters[k];
      if (0 != k)
        rebuild.append(";");
      rebuild.append(value.c_str());
    }
    const char *ptr = rebuild.c_str();
    size_t n = strlen(ptr);
    n = (n <= (m_stringWidth-1)) ? n : m_stringWidth-1;
    memcpy(buffer + i*m_stringWidth, ptr, n);
  }
  
  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", m_varname.c_str() );
  return true;

}





VCFVariableRefAllele::VCFVariableRefAllele()
{
  m_varname = REF_ALLELE;
}

bool  VCFVariableRefAllele::updateDescription( DataSetDescription *desc )
{
  return desc->addVariable(m_varname.c_str(), NC_STRING, VCF_SNP_DIM);
}


bool  VCFVariableRefAllele::populateNetCDF(int ncid,  VCF40 *vcf)
{
  return storeAllele(ncid, vcf);
}


bool VCFVariableRefAllele::storeAllele(int ncid, VCF40 *vcf)
{
  int nret;
  int varid;
  // size_t N = vcf->nSNPs;

  nret = nc_inq_varid(ncid, m_varname.c_str(), &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", m_varname.c_str());

  sspt_Array<std::string> m_array;

  // determine buffer size
  size_t n = 0;

  for (size_t i = 0; i < vcf->nSNPs; i++) {
    n += vcf->referenceAllele[i].size();
  }

  char **arrayOfStrings = new char*[ vcf->nSNPs ];
  sspt_SymbolStore buf(n + vcf->nSNPs);  //size for end of string character too.
  for (size_t i = 0; i < vcf->nSNPs; i++)
    arrayOfStrings[i] = (char *) buf.addSymbol( vcf->referenceAllele[i].c_str() );

    
  //manual is a bit unclear on meaning of startp, and countp
  size_t start[] = {0};
  size_t count[] = { vcf->nSNPs };
  nret = nc_put_vara_string(ncid, varid, start, count, (const char**) arrayOfStrings);

  if (NC_NOERR != nret) {
    fprintf(stderr, "Could not write array of strings %i\n", varid);
    return false;
  }

  delete[] arrayOfStrings;

  return true;
}







VCFVariableAltAllele::VCFVariableAltAllele(int whichAlternate)
{
  m_whichAlternate = whichAlternate;
  switch( m_whichAlternate ) {
  case 1:
    m_varname = ALT1_ALLELE;
    break;
  case 2:
    m_varname = ALT2_ALLELE;
    break;
  case 3:
    m_varname = ALT3_ALLELE;
    break;
  default:
    m_varname = "Unknown_Allele";
    m_whichAlternate = 0;  //this will lead to 0 being stored
    break;

  }
}

bool  VCFVariableAltAllele::updateDescription( DataSetDescription *desc )
{
  return desc->addVariable(m_varname.c_str(), NC_STRING, VCF_SNP_DIM);
}


bool  VCFVariableAltAllele::populateNetCDF(int ncid,  VCF40 *vcf)
{
  return storeAllele(ncid, vcf, m_whichAlternate-1);
}


bool  VCFVariableAltAllele::storeAllele(int ncid,  VCF40 *vcf, size_t column)
{
  int nret;
  int varid;
  // size_t N = vcf->nSNPs;

  nret = nc_inq_varid(ncid, m_varname.c_str(), &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", m_varname.c_str());

  // determine buffer size
  size_t n = 0;
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    std::vector<std::string> alleles = vcf->alternateAllele[i];
    if (alleles.size() <= column) {
      n += 0;
    }
    else {
      n += alleles[column].size();
    }
  }

  char **arrayOfStrings = new char*[ vcf->nSNPs ];
  sspt_SymbolStore buf(n + vcf->nSNPs);  //size for end of string character too.
  for (size_t i = 0; i < vcf->nSNPs; i++) {
    std::vector<std::string> alleles = vcf->alternateAllele[i];
    if (alleles.size() <= column) {
      arrayOfStrings[i] = (char *) buf.addSymbol( "" );
    }
    else {
      arrayOfStrings[i] = (char *) buf.addSymbol( alleles[column].c_str() );
    }

  }

    
  //manual is a bit unclear on meaning of startp, and countp
  size_t start[] = {0};
  size_t count[] = { vcf->nSNPs };
  nret = nc_put_vara_string(ncid, varid, start, count, (const char**) arrayOfStrings);

  if (NC_NOERR != nret) {
    fprintf(stderr, "Could not write array of strings %i\n", varid);
    return false;
  }

  delete[] arrayOfStrings;

  return true;
}





VCFVariablePlaceHolder::VCFVariablePlaceHolder(const char *label,  enum VCFType vcftype, bool perSNP, bool perSample, bool perString)
{
  m_varname.append(label);
  m_vcftype = vcftype;

  m_dimCode = (DimCode)(  (perSNP << 2) | (perSample << 1) | (perString) );
}




bool VCFVariablePlaceHolder::updateDescription( DataSetDescription *desc )
{
  switch (m_dimCode) {
  case DIM_NONE:
    return false;
  case DIM_STRING:
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype), VCF_STRING_DIM);
  case DIM_SAMPLE:
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype), VCF_SAMPLE_DIM);
  case DIM_SAMPLE_STRING:
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype),  VCF_SAMPLE_DIM, VCF_STRING_DIM);
  case DIM_SNP:
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype),  VCF_SNP_DIM);
  case DIM_SNP_STRING:
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype),  VCF_SNP_DIM, VCF_STRING_DIM);
  case DIM_SNP_SAMPLE:
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype),  VCF_SAMPLE_DIM, VCF_SNP_DIM);
  case DIM_ALL:
    return desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype),  VCF_SAMPLE_DIM, VCF_SNP_DIM, VCF_STRING_DIM);
  };

  return false;
}


bool VCFVariablePlaceHolder::populateNetCDF(int ncid,  VCF40 *vcf)
{
  return true;
}




VCFVariableSample::VCFVariableSample()
{
  m_vcftype = VCFVariable::VCF_CHAR;
  m_varname = SAMPLE_ID;
}


bool VCFVariableSample::updateDescription( DataSetDescription *desc )
{
  return desc->dimensionSize(VCF_STRING_DIM, &m_stringWidth)
    && desc->addVariable(m_varname.c_str(), mapVCFType(m_vcftype),  VCF_SAMPLE_DIM, VCF_STRING_DIM);
}


bool VCFVariableSample::populateNetCDF(int ncid,  VCF40 *vcf)
{
  //TODO maybe it would be better to query the netCDF file for the string width
  return storeString(ncid, vcf);
}


bool VCFVariableSample::storeString(int ncid,  VCF40 *vcf)
{
  assert(m_stringWidth > 1);
  size_t N = vcf->nSamples * m_stringWidth;

  // will know that it uses integer filter ...
  //build integer filter for specific m_vcftype ?

  int varid;
  int nret = nc_inq_varid(ncid, m_varname.c_str(), &varid);
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not find variable %s\n", m_varname.c_str());

  char *buffer = new char[N];
  memset(buffer, 0, N);

  for (size_t sampleIndex = 0; sampleIndex < vcf->nSamples; sampleIndex++) {
    const char *ptr = vcf->sampleID[sampleIndex].c_str();
    size_t n = strlen(ptr);
    n = (n <= (m_stringWidth-1)) ? n : m_stringWidth-1;
    memcpy(buffer + sampleIndex*m_stringWidth, ptr, n);
  }
  
  nret = put_var(ncid, varid, buffer);
  delete buffer;
  FALSE_ON_NETCDF_ERROR(nret, "ERROR could not write %s\n", m_varname.c_str() );
  return true;
  
}


