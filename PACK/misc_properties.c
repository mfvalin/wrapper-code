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

#if ! defined(__INTEL_COMPILER_UPDATE)
#pragma GCC optimize "tree-vectorize"
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// SIMD code seems useful for get_ieee_properties when using the Intel icc compiler
// but no longer true when using icx
#if ! defined(__INTEL_COMPILER_UPDATE)
#define NO_SIMD
#endif
#include "with_simd.h"

#include <rmn/misc_properties.h>

static inline ieee_prop_64 encode_ieee_properties(uint32_t isand,uint32_t isor, int32_t ismax, int32_t ismin, int ni, int nj ){
  ieee_prop_64 prop ;

  prop.u = 0 ;                                          // initialize all fields to 0
  prop.p.allm = isand >> 31 ;                           // all values are negative if 1
  prop.p.allp = (isor >> 31) == 0 ;                     // all values non negative if 0
  prop.p.emax = (ismax >> 23) & 0xFF ;                  // exponent of maximum value
  prop.p.emin = (ismin >> 23) & 0xFF ;                  // exponent of minimum value
  prop.p.mima = prop.p.emax == prop.p.emin ;            // same exponent for all values in block
  prop.p.mima &= (prop.p.allm | prop.p.allp) ;          // only if all values have the same sign
  prop.p.npti = ni ;
  prop.p.nptj = nj ;
  if(ni ==  8 && nj ==  8) prop.p.n_08 = 1 ;
  if(ni == 64 && nj == 64) prop.p.n_64 = 1 ;
  prop.p.errf = 0 ;                                     // no error
  return prop ;
}

// get some properties of a ni x nj contiguous block of IEEE 32 bit float values
ieee_prop_64 get_ieee_properties(void *restrict blk, int ni, int nj){
  uint32_t *restrict s = (uint32_t *) blk ;  // data will be treated as if the were unsigned integers
  int i0, i, nij = ni*nj, nij7 ;
  uint32_t ivand[8], ivor[8], isand, isor, tu ;
  int32_t  ivmax[8], ivmin[8], ismin, ismax, ts ;
  ieee_prop_64 prop = { .u = 0}  ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vo0, vmin, vmax, vor, vand, vm ;
  vm = _mm256_set1_epi32(0x7FFFFFFF) ;     // mask to suppress IEEE sign
  vor = _mm256_xor_si256(vor, vor) ;       // 0s
  vand = _mm256_cmpeq_epi32(vor, vor) ;    // 1s
  vmin = vm ;                              // largest positive number
  vmax = _mm256_sub_epi32(vor, vmin) ;     // -vmin (large negative value)
#else
  for(i=0 ; i<8 ; i++){
    ivor[i]  = 0 ;                         // 0s
    ivand[i] = ~0 ;                        // 1s
    ivmin[i] = 0x7FFFFFFF ;                // largest positive value
    ivmax[i] = -ivmin[i] ;                 // -ivmin (large negative value)
  }
#endif

  for(i0 = 0 ; i0 < nij-7 ; i0 += 8){                   // blocks of 8 values
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
    __m256i vt ;
    vo0 = _mm256_loadu_si256((__m256i *)(s+i0)) ;
    vand = _mm256_and_si256(vand, vo0) ;  // wired AND
    vor  = _mm256_or_si256(vor, vo0) ;    // wired OR
    vt   = _mm256_and_si256(vm, vo0) ;    // abs(vo0), get rid of sign bit
    vo0  = _mm256_srai_epi32(vo0, 31) ;   // propagate sign bit of vo0
    vo0  = _mm256_xor_si256(vo0, vt) ;    // 1's complement fake integer representing float
    vmax = _mm256_max_epi32(vmax, vo0) ;  // max value of fake integer
    vmin = _mm256_min_epi32(vmin, vo0) ;  // min value of fake integer
#else
    for(i=0 ; i<8 ; i++){
      tu           = s[i0+i] ;                          // get a value
      ivand[i]    &= tu ;                               // running AND
      ivor[i]     |= tu ;                               // running OR
      ts           = (int32_t) tu ;                     // make ts a signed value representation of tu
      tu          &= 0x7FFFFFFFu ;                      // abs(tu), get rid of sign bit
      ts           = tu ^ (ts >> 31) ;                  // 1's complement fake integer representing float
      ivmax[i]     = (ts > ivmax[i]) ? ts : ivmax[i] ;  // max value of fake integer
      ivmin[i]     = (ts < ivmin[i]) ? ts : ivmin[i] ;  // min value of fake integer
    }
#endif
  }
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
   _mm256_storeu_si256((__m256i *) ivmax, vmax) ;       // store SIMD registers
   _mm256_storeu_si256((__m256i *) ivmin, vmin) ;
   _mm256_storeu_si256((__m256i *) ivand, vand) ;
   _mm256_storeu_si256((__m256i *) ivor, vor) ;
#endif
  nij7 = nij & 7 ;                                      // nij modulo 7
  for(i = 0 ; i < nij7 ; i++){                          // shorter last slice (if nij not a multiple of 8)
    tu           = s[i0+i] ;                            // get a value
    ivand[i]    &= tu ;                                 // running AND
    ivor[i]     |= tu ;                                 // running OR
    ts           = (int32_t) tu ;                       // make tu a signed value
    tu          &= 0x7FFFFFFFu ;                        // get rid of sign bit
    ts           = tu ^ (ts >> 31) ;                    // 1's complement fake integer representing float
    ivmax[i]     = (ts > ivmax[i]) ? ts : ivmax[i] ;    // max value of fake integer
    ivmin[i]     = (ts < ivmin[i]) ? ts : ivmin[i] ;    // min value of fake integer
  }

  isor = ivor[0] ; isand = ivand[0] ; ismax = ivmax[0] ; ismin = ivmin[0] ;
  for(i=0 ; i<8 ; i++){                                 // fold to single values
    isand &= ivand[i] ;                                 // running AND
    isor  |= ivor[i] ;                                  // running OR
    ismax  = ivmax[i] > ismax ? ivmax[i] : ismax ;      // max value
    ismin  = ivmin[i] < ismin ? ivmin[i] : ismin ;      // min value
  }
  ismax = ismax ^ (ismax >> 31) ;                       // restore IEEE bit pattern from 1s complement fake integer
  ismin = ismin ^ (ismin >> 31) ;
  return encode_ieee_properties(isand, isor, ismax, ismin, ni, nj) ;
//   prop.u = 0 ;                                          // initialize all fields to 0
//   prop.p.allm = isand >> 31 ;                           // all values are negative if 1
//   prop.p.allp = (isor >> 31) == 0 ;                     // all values non negative if 0
//   prop.p.emax = (ismax >> 23) & 0xFF ;                  // exponent of maximum value
//   prop.p.emin = (ismin >> 23) & 0xFF ;                  // exponent of minimum value
//   prop.p.mima = prop.p.emax == prop.p.emin ;            // same exponent for all values in block
//   prop.p.mima &= (prop.p.allm | prop.p.allp) ;          // only if all values have the same sign
//   prop.p.npti = ni ;
//   prop.p.nptj = nj ;
//   if(ni ==  8 && nj ==  8) prop.p.n_08 = 1 ;
//   if(ni == 64 && nj == 64) prop.p.n_64 = 1 ;
//   prop.p.errf = 0 ;                                     // no error
//   return prop ;
}
