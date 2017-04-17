// Copyright 2017 Fred Hutchinson Cancer Research Center

#ifndef STRINGTRANSLATOR_H
#define STRINGTRANSLATOR_H





#include <stdlib.h>
#include <string.h>
#include <ctype.h>


template <typename T>
class StringTranslator {
 public:
  virtual ~StringTranslator() {}
  virtual T translate(const char *string)=0;
};




template <typename T>
class PlainTranslator : public StringTranslator<T> {
 public:
  T translate(const char *string);
 private:
};


template <>
class PlainTranslator<int> : public StringTranslator<int> {
 public:
  int translate(const char *string) { return atoi(string); }
 private:
};


template <>
class PlainTranslator<double> : public StringTranslator<double> {
 public:
  double translate(const char *string) { return atof(string); }
 private:
};





template <typename T>
class GTTranslator : public StringTranslator<T> {
 public:
  T translate(const char *string);
 private:
};


template <>
class GTTranslator<signed char> : public StringTranslator<signed char> {
 public:
  signed char translate(const char *string) {
    if (isdigit(string[0])) {
      return string[0] - '0' + 1;
    }
    else if ('.' == string[0]) {
      return 0;
    }
    else if ('|' == string[0]) {
      return 1;
    }
    else if ('\\' == string[0]) {
      return 2;
    }
    else if ('/' == string[0]) {
      return 3;
    }
    return 0;
  }
 private:
};


#define LOOKUP_MAX 256
template <typename T>
class GenoTranslator : public StringTranslator<T> {
 public:
  GenoTranslator();
  T translate(const char *string);
 private:
};


template <>
class GenoTranslator<signed char> : public StringTranslator<signed char> {
 public:

  GenoTranslator() {
   for (int i = 0; i < LOOKUP_MAX; i++)
      m_lookup[i] = 42;
    m_lookup['a'] = 1;
    m_lookup['c'] = 2;
    m_lookup['g'] = 3;
    m_lookup['t'] = 4;
    m_lookup['A'] = 1;
    m_lookup['C'] = 2;
    m_lookup['G'] = 3;
    m_lookup['T'] = 4;
    m_lookup['-'] = 0;
    m_lookup['0'] = 0;  //for plink translator
    m_lookup[0]   = 0;  //for empty strings, for prostate data set
  }

  signed char translate(const char *string) {
    return m_lookup[ (unsigned char) string[0] ];
  }
 private:
  signed char m_lookup[LOOKUP_MAX];
};



#if 0

template <typename T>
class ChromosomeTranslator : public StringTranslator<T> {
 public:
  T translate(const char *string);
 private:
};


//specialized variant

template <unsigned char>
class ChromosomeTranslator : public StringTranslator<unsigned char> {
 public:
  unsigned char translate(const char *string) {
    unsigned char value;
    if ( isalpha(string[0])) {
      if (0 == strcmp("MT", string))
        value = 26;
      else if (0 == strcmp("X", string))
        value = 23;
      else if (0 == strcmp("Y", string))
        value = 24;
      else if (0 == strcmp("XY", string))
        value = 25;
      else if (0 == strcmp("U", string)) //unknown?
        value = 0;
      else {
        fprintf(stderr, "ERROR unknown value %s for chromosome\n", string);
        exit(-1);
      }
    }
    else if (isdigit(string[0])) {
      value = atoi(string);
    }
    else {
      fprintf(stderr, "ERROR unknown value %s for chromosome\n", string);
      exit(-1);
      //return -1;
    }
    return value;
  }


 private:
};


#endif




#endif
