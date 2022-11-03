/*
 * Hopefully useful code for C and Fortran
 * Copyright (C) 2022  Recherche en Prevision Numerique
 *
 * This code is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 */
#include <stdint.h>
#include <string.h>

#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
#include <immintrin.h>
#endif

#if ! defined(STATIC)
#define STATIC static
#endif

#include <lorenzo.h>

// plain C version for cases where ni < 9
static void LorenzoPredictShort(int32_t * restrict orig, int32_t * restrict diff, int ni, int lnio, int lnid, int nj){
  int i ;
  diff[0] = orig[0] ;
  for(i=1 ; i<ni ; i++) diff[i] = orig[i] - orig[i-1] ;
  while(--nj > 0){
    diff += lnid ; 
    orig += lnio ;
    diff[0] = orig[0] - orig[0-lnio] ;
    for(i=1 ; i<ni ; i++) diff[i] = orig[i] - (orig[i-1] + orig[i-lnio] - orig[i-1-lnio]) ;
  }
}

// target is AVX2 256 bit SIMD
// predict the bottom row (1D prediction)
// row  : bottom row
// diff : prediction for bottom row
// n    : number of points in row
STATIC inline void LorenzoPredictRow0(int32_t * restrict row, int32_t * restrict diff, int n){
#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
    __m256i vi, vi1, vj1, vij1, t ;
    int i0, ii0 ;
#endif
  int i ;
  diff[0] = row[0] ;
#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
  for(ii0 = 1 ; ii0 < n ; ii0 += 8) {
    i0 = (ii0 > n-8) ? (n-8) : ii0 ;
    vi   = _mm256_loadu_si256((__m256i *) (row+i0)  ) ;
    vi1  = _mm256_loadu_si256((__m256i *) (row+i0-1)) ;
    _mm256_storeu_si256( (__m256i *) (diff+i0), _mm256_sub_epi32( vi, vi1 ) ) ;
  }
#else
  for(i=1 ; i<n ; i++) diff[i] = diff[i] = row[i] - row[i-1] ;
#endif
}

// predict row j where j > 0 (2D prediction)
// top  : row j
// bot  : row (j - 1)
// diff : prediction error for row j
// n    : number of points in row
// this function WILL NOT WORK IN-PLACE (i.e. if diff == top)
STATIC inline void LorenzoPredictRowJ(int32_t * restrict top, int32_t * restrict bot, int32_t * restrict diff, int n){
#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
    __m256i vi, vi1, vj1, vij1, t ;
    int i0, ii0 ;
#endif
  int i ;
  diff[0] = top[0] - bot[0] ;         // first point in row, 1D prediction using row below
#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
  for(ii0 = 1 ; ii0 < n ; ii0 += 8) {
    i0 = (ii0 > n-8) ? (n-8) : ii0 ;
    vi   = _mm256_loadu_si256((__m256i *) (top+i0)  ) ;   // top[i]
    vi1  = _mm256_loadu_si256((__m256i *) (top+i0-1)) ;   // top[i-1]
    vj1  = _mm256_loadu_si256((__m256i *) (bot+i0)  ) ;   // bot[i]
    vij1 = _mm256_loadu_si256((__m256i *) (bot+i0-1)) ;   // bot[i-1]
    // predicted[i,j] = z[i-1,j] + z[i,j-1] - z[i-1] = top[i-1] + bot[i] - bot[i-1]
    // diff[i,j] = orig[i,j] - predicted[i,j]
    _mm256_storeu_si256( (__m256i *) (diff+i0), _mm256_sub_epi32( vi, _mm256_sub_epi32( _mm256_add_epi32(vi1, vj1) , vij1 ) ) );
  }
#else
  for(i=1 ; i<n ; i++) diff[i] = top[i] - ( top[i-1] + bot[i] - bot[i-1] ) ;
#endif
}

// 2D lorenzo prediction (32 bit signed integers)
// orig : input : original values (32 bit signed integers)
// diff : output : original value - predicted value (using 2D Lorenzo predictor) (32 bit signed integers)
// ni   : number of useful points in row
// lnio : row storage dimension for orig
// lnid : row storage dimension for diff
// nj   : number of rows
// the SIMD version tends to be 1.5-4 times faster than the non SIMD version
// non SIMD version : Fortran anc C performances roughly equivalent (Fortran slightly faster with some compilers)
void LorenzoPredict_c(int32_t * restrict orig, int32_t * restrict diff, int ni, int lnio, int lnid, int nj){
  if(ni < 9){             // less than 9 points, SIMD version will not give correct results
    LorenzoPredictShort(orig, diff, ni, lnio, lnid, nj) ;
    return ;
  }
  LorenzoPredictRow0(orig, diff, ni) ;                // bottom row
  while(--nj > 0){
    diff += lnid ; 
    orig += lnio ;
    LorenzoPredictRowJ(orig, orig-lnio, diff, ni) ;   // all other rows
  }
}

// in place version of above function
// in order to operate in place, prediction is done backwards from top row to bottom row
void LorenzoPredictInplace_c(int32_t * restrict orig, int ni, int lnio, int nj){
  int32_t diff[ni] ;
  orig += (lnio * (nj - 1)) ;
  while(--nj > 0){                                    // all rows other than bottom row
    LorenzoPredictRowJ(orig, orig-lnio, diff, ni) ;   // predict upper row in row pair -> diff
    memcpy(orig, diff, sizeof(diff)) ;                // copy predicted row back into orig
    orig -= lnio ;                                    // next row
  }
  LorenzoPredictRow0(orig, diff, ni) ;                // bottom row
  memcpy(orig, diff, sizeof(diff)) ;
}

// restore ogiginal from 2D lorenzo prediction (32 bit signed integers)
// diff : input : original value - predicted value (32 bit signed integers)
// orig : output : restored original values from predicted differences
// ni   : number of useful points in row
// lnio : row storage dimension for orig
// lnid : row storage dimension for diff
// nj   : number of rows
// NOTE : no SIMD version exists, as the process is fully recursive
// NOTE : with the Intel compiler, it seems slower than the Fortran version
void LorenzoUnpredict_c(int32_t * restrict orig, int32_t * restrict diff, int ni, int lnio, int lnid, int nj){
  int i ;

  orig[0] = diff[0] ;                                    // restore first point of bottom row
  for(i=1 ; i<ni ; i++) orig[i] = diff[i] + orig[i-1] ;  // restore bottom row

  while(--nj > 0){
    orig += lnio ; 
    diff += lnid ;
    orig[0] = diff[0] + orig[0-lnio] ;                   // first point in row (1D prediction)
    // (original - predicted) + predicted
    for(i=1 ; i<ni ; i++) orig[i] = diff[i] + (orig[i-1] + orig[i-lnio] - orig[i-1-lnio]) ;
  }
}

void LorenzoUnpredictInplace_c(int32_t * restrict orig, int ni, int lnio, int nj){
  int i ;

  for(i=1 ; i<ni ; i++) orig[i] = orig[i] + orig[i-1] ;  // restore bottom row

  while(--nj > 0){
    orig += lnio ; 
    orig[0] = orig[0] + orig[0-lnio] ;                   // first point in row (1D prediction)
    // (original - predicted) + predicted
    for(i=1 ; i<ni ; i++) orig[i] = orig[i] + (orig[i-1] + orig[i-lnio] - orig[i-1-lnio]) ;
  }
}

