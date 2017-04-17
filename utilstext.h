// Copyright 2017 Fred Hutchinson Cancer Research Center


#ifndef UTILSTEXT_H
#define UTILSTEXT_H

#include "sspt_array.h"
#include "sspt_tmatrix.h"
#include "sspt_symbolstore.h"

//! Loads text into vectors or matrices for easier access
class UtilsText {
 public:

  static bool readVector(sspt_Array<int> *vec, const char *file);
  static bool readVector(sspt_Array<const char *> *vec, sspt_SymbolStore **buffer, const char *file);
  static bool readMatrix(sspt_TMatrix<const char *> *matrix, sspt_SymbolStore **buffer, const char *file, char delimiter, int skipLines=0);
  static bool scanColumnCount(size_t *count, const char *inputFile, char delimiter=' ');
  static bool scanRowCount(size_t *count, const char *inputFile);
  static bool scanLineWidth(size_t *value, const char *inputFile);

};

#endif
