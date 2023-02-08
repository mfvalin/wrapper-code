// Hopefully useful code for C
// Copyright (C) 2022  Recherche en Prevision Numerique
//
// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
#if ! defined(MISC_TYPES)
#define MISC_TYPES

#include <stdint.h>
#include <stddef.h>

typedef union{     // float | (un)signed 32 bit integer
  uint32_t u ;
  int32_t i ;
  float   f ;
} FloatInt;

typedef union{     // float | (un)signed 64 bit integer
  uint64_t ul ;
  int64_t l ;
  double  d ;
} DoubleLong;

typedef struct{    // pair of signed 32 bit integers
  int32_t t[2] ;
} IntPair ;

typedef struct{    // pair of floats
  float t[2] ;
} FloatPair ;

// some properties of a float array, 32 bits total
typedef struct{
  uint32_t emax:8 ,  // largest float (absolute value) exponent
           emin:8 ,  // smallest non zero float (absolute value) exponent
           allp:1 ,  // 1 if all non negative numbers
           allm:1 ,  // 1 if all negative numbers
           zero:1 ,  // 1 if zero values detected
           mima:1 ,  // 1 if same exponent and no zero
           errf:1 ,  // error flag
           n8x8:1 ,  // full 8x8 block
           xtra:2 ,  // reserved for future use
           npti:4 ,  // nuber of points in row (1-8), 0 means unknown
           nptj:4 ;  // number of rows (1-8), 0 means unknown
} ieee_prop ;

#endif
