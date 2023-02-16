// Hopefully useful code for C (memory block movers)
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

#include <stdint.h>

#include "with_simd.h"

#if ! defined(__INTEL_COMPILER_UPDATE)
#pragma GCC optimize "tree-vectorize"
#endif

// special case for row shorter than 8 elements
// insert a block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block_07(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i, ni7 ;

  if(ni > 7) return 1 ;
  ni7 = (ni & 7) ;
  while(nj--){
    for(i=0 ; i < ni7 ; i++) d[i] = s[i] ;
    s += ni ; d += lni ;
  }
  return 0 ;
}

// special case for row shorter than 8 elements
// extract a block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block_07(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i, ni7 ;

  if(ni > 7) return 1 ;
  ni7 = (ni & 7) ;
  while(nj--){
    for(i=0 ; i < ni7 ; i++) d[i] = s[i] ;
    s += lni ; d += ni ;
  }
  return 0 ;
}

// special case for row length 32
// insert a block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block_32(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni != 32) return 1 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__)
      _mm256_storeu_si256((__m256i *)(d+ 0), _mm256_loadu_si256((__m256i *)(s+ 0))) ;
      _mm256_storeu_si256((__m256i *)(d+ 8), _mm256_loadu_si256((__m256i *)(s+ 8))) ;
      _mm256_storeu_si256((__m256i *)(d+16), _mm256_loadu_si256((__m256i *)(s+16))) ;
      _mm256_storeu_si256((__m256i *)(d+24), _mm256_loadu_si256((__m256i *)(s+24))) ;
#else
    for(i=0  ; i<32 ; i++) d[i] = s[i] ;
#endif
    s += 32 ; d += lni ;
  }
  return 0 ;
}

// special case for row length 32
// extract a block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block_32(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni != 32) return 1 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__)
      _mm256_storeu_si256((__m256i *)(d+ 0), _mm256_loadu_si256((__m256i *)(s+ 0))) ;
      _mm256_storeu_si256((__m256i *)(d+ 8), _mm256_loadu_si256((__m256i *)(s+ 8))) ;
      _mm256_storeu_si256((__m256i *)(d+16), _mm256_loadu_si256((__m256i *)(s+16))) ;
      _mm256_storeu_si256((__m256i *)(d+24), _mm256_loadu_si256((__m256i *)(s+24))) ;
#else
    for(i=0  ; i<32 ; i++) d[i] = s[i] ;
#endif
    s += lni ; d += 32 ;
  }
  return 0 ;
}

// special case for row length 64
// insert a block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block_64(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni != 64) return 1 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__)
      _mm256_storeu_si256((__m256i *)(d+ 0), _mm256_loadu_si256((__m256i *)(s+ 0))) ;
      _mm256_storeu_si256((__m256i *)(d+ 8), _mm256_loadu_si256((__m256i *)(s+ 8))) ;
      _mm256_storeu_si256((__m256i *)(d+16), _mm256_loadu_si256((__m256i *)(s+16))) ;
      _mm256_storeu_si256((__m256i *)(d+24), _mm256_loadu_si256((__m256i *)(s+24))) ;
      _mm256_storeu_si256((__m256i *)(d+32), _mm256_loadu_si256((__m256i *)(s+32))) ;
      _mm256_storeu_si256((__m256i *)(d+40), _mm256_loadu_si256((__m256i *)(s+40))) ;
      _mm256_storeu_si256((__m256i *)(d+48), _mm256_loadu_si256((__m256i *)(s+48))) ;
      _mm256_storeu_si256((__m256i *)(d+56), _mm256_loadu_si256((__m256i *)(s+56))) ;
#else
    for(i=0  ; i<32 ; i++) d[i] = s[i] ;
    for(i=32 ; i<64 ; i++) d[i] = s[i] ;
#endif
    s += 64 ; d += lni ;
  }
  return 0 ;
}

// special case for row length 64
// extract a block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block_64(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni != 64) return 1 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__)
      _mm256_storeu_si256((__m256i *)(d+ 0), _mm256_loadu_si256((__m256i *)(s+ 0))) ;
      _mm256_storeu_si256((__m256i *)(d+ 8), _mm256_loadu_si256((__m256i *)(s+ 8))) ;
      _mm256_storeu_si256((__m256i *)(d+16), _mm256_loadu_si256((__m256i *)(s+16))) ;
      _mm256_storeu_si256((__m256i *)(d+24), _mm256_loadu_si256((__m256i *)(s+24))) ;
      _mm256_storeu_si256((__m256i *)(d+32), _mm256_loadu_si256((__m256i *)(s+32))) ;
      _mm256_storeu_si256((__m256i *)(d+40), _mm256_loadu_si256((__m256i *)(s+40))) ;
      _mm256_storeu_si256((__m256i *)(d+48), _mm256_loadu_si256((__m256i *)(s+48))) ;
      _mm256_storeu_si256((__m256i *)(d+56), _mm256_loadu_si256((__m256i *)(s+56))) ;
#else
    for(i=0  ; i<32 ; i++) d[i] = s[i] ;
    for(i=32 ; i<64 ; i++) d[i] = s[i] ;
#endif
    s += lni ; d += 64 ;
  }
  return 0 ;
}

// insert a block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni  <  8) return put_word_block_07(f, blk, ni, lni, nj) ;
  if(ni == 32) return put_word_block_32(f, blk, ni, lni, nj) ;
  if(ni == 64) return put_word_block_64(f, blk, ni, lni, nj) ;

  ni7 = (ni & 7) ;
  ni7 = ni7 ? ni7 : 8 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__)
      _mm256_storeu_si256((__m256i *)(d), _mm256_loadu_si256((__m256i *)(s))) ;
    for(i0 = ni7 ; i0 < ni-7 ; i0 += 8 ){
      _mm256_storeu_si256((__m256i *)(d+i0), _mm256_loadu_si256((__m256i *)(s+i0))) ;
    }
#else
    for(i=0 ; i<8 ; i++) d[i] = s[i] ;
    for(i0 = ni7 ; i0 < ni-7 ; i0 += 8 ){
      for(i=0 ; i<8 ; i++) d[i0+i] = s[i0+i] ;
    }
#endif
    s += ni ; d += lni ;
  }
  return 0 ;
}

// extract a block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni  <  8) return get_word_block_07(f, blk, ni, lni, nj) ;
  if(ni == 32) return get_word_block_32(f, blk, ni, lni, nj) ;
  if(ni == 64) return get_word_block_64(f, blk, ni, lni, nj) ;

  ni7 = (ni & 7) ;
  ni7 = ni7 ? ni7 : 8 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__)
      _mm256_storeu_si256((__m256i *)(d), _mm256_loadu_si256((__m256i *)(s))) ;
    for(i0 = ni7 ; i0 < ni-7 ; i0 += 8 ){
      _mm256_storeu_si256((__m256i *)(d+i0), _mm256_loadu_si256((__m256i *)(s+i0))) ;
    }
#else
    for(i=0 ; i<8 ; i++) d[i] = s[i] ;
    for(i0 = ni7 ; i0 < ni-7 ; i0 += 8 ){
      for(i=0 ; i<8 ; i++) d[i0+i] = s[i0+i] ;
    }
#endif
    s += lni ; d += ni ;
  }
  return 0 ;
}
