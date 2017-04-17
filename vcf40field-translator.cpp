// Copyright 2017 Fred Hutchinson Cancer Research Center



#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <stdlib.h>


#include "sspt_cord.h"
#include "sspt_skewheap.h"
#include "sspt_delimiterparse.h"
#include "sspt_avltree.h"

#include "vcf40field-translator.h"

#include "genovariable.h"
#include "chromosomeposition.h"
#include "utilsnetcdf.h"
#include "utilstext.h"

#include "datasetdescription.h"

#include "vcf_names.h"



//if duplicates are allowed ...
struct ChromosomePriority {
  int chromosome;
  int position;
  sspt_Cord name;
};

class PositionCompare {
public:


  //  bool operator()(const ChromosomePriority* &a, const  ChromosomePriority* &b) const {
  bool operator()(const ChromosomePriority *a, const ChromosomePriority *b) const {
    return ( (a->chromosome < b->chromosome) 
             || ((a->chromosome == b->chromosome) &&  (a->position < b->position))
             || ((a->chromosome == b->chromosome) &&  (a->position == b->position) && (a->name < b->name))
             );
  }


};


//SNP Name to index data structure
//Sample ID to index data structure

//column labels




#define TABLE_DEFAULT 64
VCF40FieldTranslator::VCF40FieldTranslator() : m_variableTable(TABLE_DEFAULT)
{

  m_buffer = 0;
  m_autofilter = false;
}


// input looks like:   <ID=PV4,Number=4,Type=Float,Descri ... >
bool VCF40FieldTranslator::extractVariableInfo(std::string *label, std::string *vtype, std::string *number, const char *item)
{
  sspt_DelimiterParse a(item, '<', false);
  sspt_DelimiterParse b( a.value(a.values()-1), '>', false);
  sspt_DelimiterParse keypair(b.value(0), ',', false);

  //printf("%s\n", b.value(0));

  int found = 0;
  for (size_t i = 0; i < keypair.values(); i++) {
    sspt_DelimiterParse p( keypair.value(i), '=', false);
    //printf("%s\n", keypair.value(i));
    if (p.values() != 2)
      continue;
    else if (0 == strcmp(p.value(0), "ID")) {
      *label = p.value(1);
      found++;
    }
    else if (0 == strcmp(p.value(0), "Number")) {
      *number = p.value(1);
      found++;
    }
    else if (0 == strcmp(p.value(0), "Type")) {
      *vtype = p.value(1);
      found++;
    }
  }

  return 3==found;
}


void VCF40FieldTranslator::selectVariables(VCF40 *vcf)
{

  //is this needed??
  //read filters from key-value pairs
  for (std::multimap<std::string, std::string>::iterator iter = vcf->headerPairs.find("FILTER");
       iter != vcf->headerPairs.upper_bound("FILTER") 
       && iter != vcf->headerPairs.end(); ++iter) {
    std::pair<std::string, std::string> pair = *iter;
    sspt_DelimiterParse p(pair.second.c_str(), ',', false );
    //printf("FILTER -- %s\n", p.value(0));
  }


  //read info from key-value pairs
  for (std::multimap<std::string, std::string>::iterator iter = vcf->headerPairs.find("INFO");
       iter != vcf->headerPairs.upper_bound("INFO") 
       && iter != vcf->headerPairs.end(); ++iter) {
    std::pair<std::string, std::string> pair = *iter;
    std::string label, vtype, number;
    extractVariableInfo(&label, &vtype, &number, pair.second.c_str());
    //printf("INFO -- %s,%s,%s\n", label.c_str(), vtype.c_str(), number.c_str());
    addInfoVar(label, vtype, number);
  }


  //read format from key-value pairs
  for (std::multimap<std::string, std::string>::iterator iter = vcf->headerPairs.find("FORMAT");
       iter != vcf->headerPairs.upper_bound("FORMAT") 
       && iter != vcf->headerPairs.end(); ++iter) {
    std::pair<std::string, std::string> pair = *iter;
    std::string label, vtype, number;
    extractVariableInfo(&label, &vtype, &number, pair.second.c_str());
    printf("FORMAT -- %s,%s,%s\n", label.c_str(), vtype.c_str(), number.c_str());
    addFormatVar(label, vtype, number);
  }



  for (std::multimap<std::string, std::string>::iterator iter = vcf->headerPairs.begin();
       iter != vcf->headerPairs.end(); ++iter) {
    std::pair<std::string, std::string> pair = *iter;
    //printf("header -- %s  =  %s\n", pair.first.c_str(), pair.second.c_str());
  }
 

  //Chromosome and Position
  {
    VCFVariable *var = new VCFVariableLocation();
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }

  {
    VCFVariable *var = new VCFVariableID();
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }


  {
    VCFVariable *var = new VCFVariableQuality();
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }

  if (m_autofilter) {
    VCFVariable *var = new VCFVariableAutoFilter(vcf);
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }
  else {
    VCFVariable *var = new VCFVariableSimpleFilter();
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }

  //add support for ref and alt alleles
  {
    VCFVariable *var = new VCFVariableRefAllele();
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }

  for (int i = 1; i <= 3; i++) {
    VCFVariable *var = new VCFVariableAltAllele(i);
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }

  {
    VCFVariable *var = new VCFVariablePlaceHolder(SNP_NAME, VCFVariable::VCF_STRING, true, false, true);
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }

  {
    VCFVariable *var = new VCFVariableSample(); //SAMPLE_ID, VCFVariable::VCF_STRING, false, true, true);
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }


  {
    VCFVariable *var = new VCFVariablePlaceHolder(GENOTYPE, VCFVariable::VCF_BYTE, true, true, false);
    sspt_Cord name;
    var->variableName(&name);
    m_variableTable.insert(name, var);
  }


  //hack to ignore some variables like PL3
  {
    sspt_Cord name( "array_PL3" );
    if (m_variableTable.isKey(name))
      m_variableTable.remove(name);
  }


}


VCF40FieldTranslator::~VCF40FieldTranslator()
{
  if (!m_buffer)
    delete[] m_buffer;

}



bool VCF40FieldTranslator::addInfoVar(const std::string &label, const std::string &vtype, const std::string &number)
{

  VCFVariable::VCFType vcftype;

  if (vtype == "Integer") {
    vcftype = VCFVariable::VCF_INT;
  }
  else if (vtype == "Float") {
    vcftype = VCFVariable::VCF_DOUBLE;
  }
  else
    return false;

  VCFVariable *var = new VCFVariableColumnInfo(label.c_str(), vcftype, atoi(number.c_str()) );

  sspt_Cord name;
  var->variableName(&name);

  m_variableTable.insert(name, var);
  return true;
}

bool VCF40FieldTranslator::addFormatVar(const std::string &label, const std::string &vtype, const std::string &number)
{

  VCFVariable::VCFType vcftype;
  int n = atoi(number.c_str());

  //workaround
  if (-1 == n && "PL" == label) {
    //there seems to be only one alternate sequence, so that's three numbers for the likelihood pairs
    n = 3;
  }
  //else if ("GL" == label) {
  //  return true; //for the initial indel vcf, it doesn't actually have GL in the data!
  //}

  //end workaround

  VCFVariable *var = 0;

  if (vtype == "Integer") {
    vcftype = VCFVariable::VCF_INT;
    var = new VCFVariableColumnFormat(label.c_str(), vcftype, n );
  }
  else if (vtype == "Float") {
    vcftype = VCFVariable::VCF_DOUBLE;
    var = new VCFVariableColumnFormat(label.c_str(), vcftype, n );
  }
  else if (vtype == "String" && label == "GT") {
    //vcftype = VCFVariable::VCF_AB;
    //n = 2;
    //var = new VCFVariableColumnFormat(label.c_str(), vcftype, n, '/');
    var = new VCFVariableGenotype(label.c_str());
  }
  else
    return false;

  //VCFVariable *var = new VCFVariableColumnInfo(label.c_str(), vcftype, n );
  sspt_Cord name;
  var->variableName(&name);

  m_variableTable.insert(name, var);
  return true;
}





bool VCF40FieldTranslator::process(const char *outputFile, 
                                   VCF40 *vcf, VCF40 *alternateHeader, bool sortSNPs)
{
  //m_sortingByChromosome = sortingByChromosome;
  //m_allowDuplicates = allowDuplicatePositions;

  //phase zero - determine number of snps and samples
  m_nSamples = vcf->nSamples;
  m_nSNPs = vcf->nSNPs;
  printf("Samples %zu SNPs %zu\n", m_nSamples, m_nSNPs);



  //phase one - allocate variables, with sample as unlimited dimension
  if (0 != alternateHeader) {
    selectVariables(alternateHeader);
  }
  else {
    selectVariables(vcf);
  }

  printf("parsed variables from VCF header\n");
  
  DataSetDescription *desc;

  if (!createDescription(&desc, vcf)) {
    return false;
  }

  printf("created description\n");

  if (!createNetcdfVariables(outputFile, desc)) {
    return false;
  }


  printf("created netCDF variables\n");



  //phase two - iterate over vcf data
  if (!processVCF(vcf, sortSNPs)) {
    return false;
  }


  nc_close(m_ncid);
  printf("finished files\n");

  //additional phases happen with other programs, transpose, merge, rs-write, genotype-write, ?position mapping
  //position-sort?


  return true;
}



bool VCF40FieldTranslator::createDescription(DataSetDescription **output, VCF40 *vcf)
{
  

  sspt_List<VCFVariable*> list;
  m_variableTable.contents(&list);

  DataSetDescription *desc = new DataSetDescription;

  if (!desc->addDimension(VCF_SAMPLE_DIM, vcf->nSamples))
    return false;

  if (!desc->addDimension(VCF_SNP_DIM, vcf->nSNPs))
    return false;

  if (!desc->addDimension(VCF_STRING_DIM, MAX_STRING))
    return false;



  //create other 2D arrays
  for (sspt_ListIterator<VCFVariable*> iter = list.begin(); !iter.atEnd(); iter.moveNext()) {
    VCFVariable *v = iter.current();
    if (!v->updateDescription(desc))
      return false;
  }

  *output = desc;
  return true;
}



bool VCF40FieldTranslator::createNetcdfVariables(const char *outputFile, DataSetDescription *desc)
{

  if (!desc->createEmptyNetCDF(&m_ncid, outputFile)) {
    return false;
  }
  return true;

}




#define LOOKUP(var, a) {                             \
    sspt_Cord name(a);                               \
    if (!m_variableTable.find(name, var)) {          \
      fprintf(stderr, "ERROR could not look up %s\n", name.c_str()); \
      return false;                                  \
    }                                                \
  }



bool VCF40FieldTranslator::sortByChromosome( sspt_Array<int> *mapping )
{
  GenoVariable *snpName;
  GenoVariable *chromosome;
  GenoVariable *position;

  //LOOKUP(&snpName, SNP_NAME);
  //LOOKUP(&chromosome, CHROMOSOME);
  //LOOKUP(&position, POSITION);

  printf("Creating mapping\n");
  sspt_Array<int> c(m_nSNPs);
  sspt_Array<int> p(m_nSNPs);
  *mapping = sspt_Array<int>(m_nSNPs);

  chromosome->readArray(&c);
  position->readArray(&p);

  //sort
  bool allowDuplicates = false;
  if (!allowDuplicates) {
    sspt_AVLTree<ChromosomePosition, int> table;
    for (size_t k = 0; k < c.size(); k++) {
      if (!table.insert( ChromosomePosition(c[k], p[k]), k) ) {
        fprintf(stderr, "ERROR duplicate chromosome position %i %i\n", c[k], p[k]);
        return false;
      }
    }

    int index = 0;
    for (sspt_AVLIterator<ChromosomePosition, int> iter = table.begin();
         !iter.atEnd(); iter.moveNext(), index++) {
      (*mapping)[ iter.data() ] = index;
    }
  }
  else {  //need to do something different if allowing duplicates
    sspt_SkewHeap<ChromosomePriority*, int, PositionCompare> heap;

    //get strings for names ...
    for (size_t k = 0; k < c.size(); k++) {
      ChromosomePriority *P = new ChromosomePriority;
      P->chromosome = c[k];
      P->position = p[k];
      P->name = snpName->readString(k);
      heap.insert(P, k);
    }

    int index = 0;
    while (!heap.isEmpty()) {
      ChromosomePriority *P;
      int originalIndex;
      heap.remove(&P, &originalIndex);
      (*mapping)[ originalIndex ] = index;
      delete P;
      index++;
    }
  }

  printf("Done creating mapping\n");
  return true;
}





bool VCF40FieldTranslator::parseGenotype(std::string *first, std::string *second, const std::string line)
{
  //find separator token
  size_t i;
  for (i = 0; i < line.size(); i++) {
    if ('/' == line[i] ||
        '|' == line[i] ||
        '\\' == line[i]  ) {
      break;
    }
  }


  if (i == line.size()) {
    fprintf(stderr, "ERROR could not find '=' token\n" );
    return false;
  }
  *first = line.substr(0, i);
  *second = line.substr(i+1, line.size()-(i+1));

  return true;
}





bool VCF40FieldTranslator::processVCF(VCF40 *vcf, bool sortSNPs)
{

  sspt_List<VCFVariable*> list;
  m_variableTable.contents(&list);

  for (sspt_ListIterator<VCFVariable*> iter = list.begin(); !iter.atEnd(); iter.moveNext()) {
    VCFVariable *v = iter.current();
    sspt_Cord name;
    v->variableName(&name);

    printf("processing %s ...\n", name.c_str());

    if ( !v->populateNetCDF(m_ncid, vcf) )
      return false;
  }


  return true;

}



  





