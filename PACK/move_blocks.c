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
#include <rmn/misc_operators.h>

// SIMD does not seem to be useful any more for these funtions
#undef WITH_SIMD

#if ! defined(__INTEL_COMPILER_UPDATE)
#pragma GCC optimize "tree-vectorize"
#endif

// special case for rows shorter than 8 elements
// insert a contiguous block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block_07(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i, ni7 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vm = _mm256_memmask_si256(ni) ;  // mask for load and store operations
#endif

  if(ni > 7) return 1 ;
  ni7 = (ni & 7) ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ;
#else
    for(i=0 ; i < ni7 ; i++) d[i] = s[i] ;
#endif
    s += ni ; d += lni ;
  }
  return 0 ;
}

// special case for rows shorter than 8 elements
// extract a contiguous block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block_07(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i, ni7 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vm = _mm256_memmask_si256(ni) ;  // mask for load and store operations
#endif

  if(ni > 7) return 1 ;
  ni7 = (ni & 7) ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ;
#else
    for(i=0 ; i < ni7 ; i++) d[i] = s[i] ;
#endif
    s += lni ; d += ni ;
  }
  return 0 ;
}

// specialized versions for ni = 8 / 32 / 64 (not needed for now)
#if 0
// special case for row length 8
// insert a contiguous block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block_08(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i0, i, j ;

  if(ni != 8) return 1 ;
  if(nj == 8){
    for(j=0 ; j<8 ; j++){
      for(i=0  ; i<8 ; i++) d[i] = s[i] ;
      s += 8 ; d += lni ;
    }
  }else{
    while(nj--){
      for(i=0  ; i<8 ; i++) d[i] = s[i] ;
      s += 8 ; d += lni ;
    }
  }
  return 0 ;
}

// special case for row length 8
// extract a contiguous block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block_08(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i0, i, j ;

  if(ni != 8) return 1 ;
  if(nj == 8){
    for(j=0 ; j<8 ; j++){
      for(i=0  ; i<8 ; i++) d[i] = s[i] ;
      s += lni ; d += 8 ;
    }
  }else{
    while(nj--){
      for(i=0  ; i<8 ; i++) d[i] = s[i] ;
      s += lni ; d += 8 ;
    }
  }
  return 0 ;
}

// special case for row length 32
// insert a contiguous block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block_32(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i0, i ;

  if(ni != 32) return 1 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
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
// extract a contiguous block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block_32(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i0, i ;

  if(ni != 32) return 1 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
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
// insert a contiguous block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block_64(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni != 64) return 1 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
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
// extract a contiguous block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block_64(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni != 64) return 1 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
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
#endif

// insert a contiguous block (ni x nj) of 32 bit words into f from blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int put_word_block(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni  <  8) return put_word_block_07(f, blk, ni, lni, nj) ;
//   if(ni ==  8) return put_word_block_08(f, blk, ni, lni, nj) ;
//   if(ni == 32) return put_word_block_32(f, blk, ni, lni, nj) ;
//   if(ni == 64) return put_word_block_64(f, blk, ni, lni, nj) ;

  ni7 = (ni & 7) ;
  ni7 = ni7 ? ni7 : 8 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
      _mm256_storeu_si256((__m256i *)(d), _mm256_loadu_si256((__m256i *)(s))) ;
    for(i0 = ni7 ; i0 < ni-7 ; i0 += 8 ){
      _mm256_storeu_si256((__m256i *)(d+i0), _mm256_loadu_si256((__m256i *)(s+i0))) ;
    }
#else
    for(i=0 ; i<8 ; i++) d[i] = s[i] ;     // first and second chunk may overlap if ni not a multiple of 8
    for(i0 = ni7 ; i0 < ni-7 ; i0 += 8 ){
      for(i=0 ; i<8 ; i++) d[i0+i] = s[i0+i] ;
    }
#endif
    s += ni ; d += lni ;
  }
  return 0 ;
}

// extract a contiguous block (ni x nj) of 32 bit words from f into blk
// ni    : row size (row storage size in blk)
// lni   : row storage size in f
// nj    : number of rows
int get_word_block(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i0, i, ni7 ;

  if(ni  <  8) return get_word_block_07(f, blk, ni, lni, nj) ;
//   if(ni ==  8) return get_word_block_08(f, blk, ni, lni, nj) ;
//   if(ni == 32) return get_word_block_32(f, blk, ni, lni, nj) ;
//   if(ni == 64) return get_word_block_64(f, blk, ni, lni, nj) ;

  ni7 = (ni & 7) ;
  ni7 = ni7 ? ni7 : 8 ;
  while(nj--){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
      _mm256_storeu_si256((__m256i *)(d), _mm256_loadu_si256((__m256i *)(s))) ;
    for(i0 = ni7 ; i0 < ni-7 ; i0 += 8 ){
      _mm256_storeu_si256((__m256i *)(d+i0), _mm256_loadu_si256((__m256i *)(s+i0))) ;
    }
#else
    for(i=0 ; i<8 ; i++) d[i] = s[i] ;     // first and second chunk may overlap if ni not a multiple of 8
    for(i0 = ni7 ; i0 < ni-7 ; i0 += 8 ){
      for(i=0 ; i<8 ; i++) d[i0+i] = s[i0+i] ;
    }
#endif
    s += lni ; d += ni ;
  }
  return 0 ;
}
