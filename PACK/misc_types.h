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

typedef union{
  uint32_t i ;
  float    f ;
} FloatUint ;

typedef union{
  int32_t i ;
  float   f ;
} FloatInt;

typedef union{
  uint64_t l ;
  double   d ;
} DoubleUlong;

typedef union{
  int64_t l ;
  double  d ;
} DoubleLong;

#endif
