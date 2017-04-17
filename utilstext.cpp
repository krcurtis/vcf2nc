// Copyright 2017 Fred Hutchinson Cancer Research Center



#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "sspt_list.h"
#include "sspt_datachunk.h"
#include "sspt_delimiterparse.h"
#include "utilstext.h"

bool UtilsText::readVector(sspt_Array<int> *vec, const char *file)
{

  FILE *fptr = fopen(file, "rb");
  if (0 == fptr) {
    fprintf(stderr, "ERROR Could not open %s\n", file);
    return false;
  }

  const int n = 128;
  char line[n];
  sspt_List<int> list;
  while (0 != fgets(line, n, fptr)) {
    list.insertRear( atoi(line) );
  }
  fclose(fptr);

  int index = 0;
  *vec = sspt_Array<int>( list.size() );
  for (sspt_ListIterator<int> iter = list.begin(); !iter.atEnd(); iter.moveNext(), index++) {
    (*vec)[index] = iter.current();
  }

  return true;
}


bool UtilsText::readVector(sspt_Array<const char *> *vec, sspt_SymbolStore **buffer, const char *file)
{

  FILE *fptr = fopen(file, "rb");
  if (0 == fptr) {
    fprintf(stderr, "ERROR Could not open %s\n", file);
    return false;
  }
  fclose(fptr);

  sspt_DataChunk block(file);

  if (0 == block.size()) {
    *vec =  sspt_Array<const char*>(0);
    *buffer = new sspt_SymbolStore(0);
    return true;
  }

  *buffer = new sspt_SymbolStore( 1024 + (int) block.size() );
  sspt_DelimiterParse p( (const char *) block.contents(), '\n', true  );

  sspt_List<const char *> list;
  for (size_t i = 0; i < p.values(); i++) {
    if (0 == strlen( p.value(i) ) )
      continue;
    list.insertRear( (*buffer)-> addSymbol( p.value(i) )  );
  }

  int index = 0;
  *vec = sspt_Array<const char *>( list.size() );
  for (sspt_ListIterator<const char *> iter = list.begin(); !iter.atEnd(); iter.moveNext(), index++) {
    (*vec)[index] = iter.current();
  }

  return true;
}

bool UtilsText::readMatrix(sspt_TMatrix<const char *> *matrix, sspt_SymbolStore **buffer, const char *file, char delimiter, int skipLines)
{
  sspt_DataChunk block(file);

  if (0 == block.size()) {
    (*matrix) = sspt_TMatrix<const char*>(0,0);
    *buffer = new sspt_SymbolStore(0);
    return true;
  }

  *buffer = new sspt_SymbolStore( 1024 + (int) block.size() );
  sspt_DelimiterParse lines( (const char *) block.contents(), '\n', false );

  size_t nCols=0;
  {
    sspt_DelimiterParse columns( lines.value(skipLines), delimiter, false);
    nCols = columns.values();
  }

  size_t nRows = 0;
  for (size_t i = skipLines; i < lines.values(); i++) {
    if (0 == strlen( lines.value(i) ) )
      continue;
    nRows++;
  }

  *matrix = sspt_TMatrix<const char *>( nRows, nCols );
  size_t index = 0;
  for (size_t i = skipLines; i < lines.values(); i++) {
    if (0 == strlen( lines.value(i) ) )
      continue;
    sspt_DelimiterParse columns( lines.value(i), delimiter, false);
    if (nCols != columns.values()) {
      fprintf(stderr, "ERROR inconsistent number of columns (expected %zu, found %zu) at line %zu\n", nCols, columns.values(), i);
      return false;
    }
    
    for (size_t k = 0; k < nCols; k++) {
      (*matrix)(index,k) =  (*buffer)-> addSymbol( columns.value(k) );
    }
    index++;
  }


  return true;
}




bool UtilsText::scanColumnCount(size_t *count, const char *inputFile, char delimiter)
{
  FILE *fptr = fopen(inputFile, "rb");
  if (0 == fptr) {
    fprintf(stderr, "ERROR could not open %s\n", inputFile);
    return false;
  }

  size_t byteCount = 0;
  int value;
  while ( (value = fgetc(fptr)) != EOF) {
    byteCount++;
    if ('\n' == (char) value) {
      break;
    }
  }

  if (EOF == value) {  //invalid file format
    fclose(fptr);
    return false;
  }

  fseeko(fptr, 0, SEEK_SET);
  char *buffer = new char[byteCount+1];  //add one for null terminator
    
  fread(buffer, byteCount, 1, fptr);
  fclose(fptr);

  buffer[byteCount]=0;
  assert( '\n' == buffer[byteCount-1] ); //should have newline at end



  sspt_DelimiterParse line(buffer, delimiter, false);  
  *count = (size_t) line.values();
  delete[] buffer;

  return true;
}


bool UtilsText::scanRowCount(size_t *count, const char *inputFile)
{

  sspt_DataChunk block(inputFile);

  if (0 == block.size()) {
    *count = 0;
    return true;
  }

  *count = 0;
  const char *content = (const char *) block.contents();
  for (unsigned long i = 0; i < block.size(); i++) {
    if ('\n' == content[i]) {
      (*count)++;
    }
  }

  return true;
}


bool UtilsText::scanLineWidth(size_t *value, const char *inputFile)
{
  sspt_DataChunk block(inputFile);

  if (0 == block.size()) {
    *value = 0;
    return true;
  }

  size_t max = 0;
  size_t count = 0;
  const char *content = (const char *) block.contents();
  for (size_t i = 0; i < block.size(); i++) {
    if ('\n' != content[i]) {
      count++;
    }
    else {
      count++;
      if (count > max)
        max = count;
      count = 0;
    }
  }

  *value = max;
  return true;
}

