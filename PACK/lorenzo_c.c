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

#include <with_simd.h>

#if ! defined(STATIC)
#define STATIC static
#endif

// copy n 32 bit words from array s0 to array d0 with strong SIMD hint
// s0 and d0 MUST NOT OVERLAP
// should be faster than memcpy for small to medium lengths
// n MUST BE >= 8
// static inline void mem_cpy_w32(void * restrict d0, void * restrict s0, int n){
//   int32_t * restrict s = (int32_t *)s0, * restrict d = (int32_t *)d0 ;
//   int i, ni7, i0 ;
//   if(n < 8) {   // less than SIMD vector length
//     for(i = 0 ; i < (n & 7) ; i++) d[i] = s[i] ;
//   }else{
//     ni7 = (n & 7) ? (n & 7) : 8 ;
//     for(i0=0 ; i0<n ; i0+=ni7 , ni7 = 8){
//       for(i=0 ; i<8 ; i++) d[i0+i] = s[i0+i] ;  // 8 element SIMD hint
//     }
//   }
// }

#include <rmn/lorenzo.h>

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
  if(n < 9){   // the SIMD version will not work for n < 9
    for(i=1 ; i<n ; i++) diff[i] = diff[i] = row[i] - row[i-1] ;
    return ;
  }
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

// bottom row, n < 8, in place transform
static inline void LorenzoPredictRow0_inplace_07(int32_t * restrict row, int n){
  int i, n7 = (n & 7) ;
  int32_t r[7] ;
  r[0] = row[0] ;
  for(i=1 ; i<n7 ; i++) r[i] = row[i] - row[i-1] ;
  for(i=0 ; i<n7 ; i++) row[i] = r[i] ;
}

// bottom row, in place transform
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)

STATIC inline void LorenzoPredictRow0_inplace(int32_t * restrict row, int n){
  int i, i0, j0, n7 = (n & 7) ;
  __m256i vi, vi1, vt, vr, v0, vs ;

  if(n < 8) {
    LorenzoPredictRow0_inplace_07(row, n) ;
    return ;
  }

  n7 = n7 ? n7 : 8 ;
  v0   = _mm256_setzero_si256() ;                         // 0s
  vs   = _mm256_cvtepi8_epi32( _mm_set1_epi64x(0x0605040302010000lu) ) ;  // shuffle patterm
  // first 8 elements
  vi   = _mm256_loadu_si256((__m256i *) (row)  ) ;        // row[0:7]
  vi1  = _mm256_permutevar8x32_epi32(vi, vs) ;            // 0, 0, 1, 2, 3, 4, 5, 6  permutation
  vi1  = _mm256_blend_epi32(vi1, v0, 1) ;                 // first element set to 0 (row[-1])
  // row[i] - row[i-1]
  vt   = _mm256_sub_epi32( vi, vi1 ) ;
  vr   = vt ;
  j0 = 0 ;
  // chunks of 8 elements (second chunk may overlap first chunk)
  for(i0=n7 ; i0<n ; i0+=8){
    vi   = _mm256_loadu_si256((__m256i *) (row+i0)  ) ;   // row[i0:i0+7]
#if defined(WITH_FETCH)
    vi1  = _mm256_loadu_si256((__m256i *) (row+i0-1)) ;   // row[i0-1:i0+6]
#else
    v0   = _mm256_set1_epi32(row[i0-1]) ;                 // row[i0-1]
    vi1  = _mm256_permutevar8x32_epi32(vi, vs) ;          // 0, 0, 1, 2, 3, 4, 5, 6  permutation
    vi1  = _mm256_blend_epi32(vi1, v0, 1) ;               // first element set to row[i0-1]
#endif
    // row[i0+i] - row[i0+i-1]
    vt   = _mm256_sub_epi32( vi, vi1 ) ;
    _mm256_storeu_si256( (__m256i *) (row+j0), vr) ;      // delayed store, result of previous pass
    vr = vt ;
    j0 = i0 ;
  }
  _mm256_storeu_si256( (__m256i *) (row+j0), vr) ;        // delayed store, result of previous pass
}

#else

STATIC inline void LorenzoPredictRow0_inplace(int32_t * restrict row, int n){
  int i, i0, j0, n7 = (n & 7) ;
  int32_t t[8], r[8], r0 ;

  if(n < 8) {
    LorenzoPredictRow0_inplace_07(row, n) ;
    return ;
  }
  r0 = row[0] ;
  n7 = n7 ? n7 : 8 ;
  // first 8 elements
  for(i=0 ; i<8 ; i++) r[i] = row[i] - row[i-1] ;
  j0 = 0 ;
  // chunks of 8 elements (second chunk may overlap first chunk)
  for(i0=n7 ; i0<n ; i0+=8){
    for(i=0 ; i<8 ; i++) t[i] = row[i0+i] - row[i0+i-1] ;
    for(i=0 ; i<8 ; i++) row[j0+i] = r[i] ;
    for(i=0 ; i<8 ; i++) r[i] = t[i] ;
    j0 = i0 ;
  }
  for(i=0 ; i<8 ; i++) row[j0+i] = r[i] ;
  row[0] = r0 ;
}

#endif

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
  if(n < 9){   // the SIMD version will not work for n < 9
    for(i=1 ; i<n ; i++) diff[i] = top[i] - ( top[i-1] + bot[i] - bot[i-1] ) ;
    return ;
  }
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

// all rows but bottom row, n < 8, in place transform
static void LorenzoPredictRowJ_inplace_07(int32_t * restrict top, int32_t * restrict bot, int n){
  int i, n7 = (n & 7) ;
  int32_t r[7] ;
  r[0] = top[0] - bot[0] ;
  for(i=1 ; i<n7 ; i++) r[i] = top[i] - ( top[i-1] + bot[i] - bot[i-1] ) ;
  for(i=0 ; i<n7 ; i++) top[i] = r[i] ;
}

// all rows but bottom row, in place transform
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)

STATIC inline void LorenzoPredictRowJ_inplace(int32_t * restrict top, int32_t * restrict bot, int n){
  int i, i0, j0, n7 = (n & 7) ;
  __m256i vi, vi1, vj1, vij1, vt, vr, v0, vs ;

  if(n < 8) {
    LorenzoPredictRowJ_inplace_07(top, bot, n) ;
    return ;
  }
  n7 = n7 ? n7 : 8 ;
  v0   = _mm256_setzero_si256() ;                         // 0s
  vs   = _mm256_cvtepi8_epi32( _mm_set1_epi64x(0x0605040302010000lu) ) ;  // shuffle patterm
  // first 8 elements
  vi   = _mm256_loadu_si256((__m256i *) (top)  ) ;        // top[0:7]
  vj1  = _mm256_loadu_si256((__m256i *) (bot)  ) ;        // bot[0:7]
//   vi1  = _mm256_loadu_si256((__m256i *) (top-1)) ;        // top[-1:6]
  vi1  = _mm256_permutevar8x32_epi32(vi, vs) ;            // 0, 0, 1, 2, 3, 4, 5, 6  permutation
  vi1  = _mm256_blend_epi32(vi1, v0, 1) ;                 // first element set to 0 (top[-1])
//   vij1 = _mm256_loadu_si256((__m256i *) (bot-1)) ;        // bot[-1:6]
  vij1 = _mm256_permutevar8x32_epi32(vj1, vs) ;            // 0, 0, 1, 2, 3, 4, 5, 6  permutation
  vij1 = _mm256_blend_epi32(vij1, v0, 1) ;                 // first element set to 0 (bot[-1])
  // top[i] - ( top[i-1] + bot[i] - bot[i-1] )
  vt   = _mm256_sub_epi32( vi, _mm256_sub_epi32( _mm256_add_epi32(vi1, vj1) , vij1 ) ) ;
//   vr   =  _mm256_xor_si256( _mm256_slli_epi32(vt, 1) , _mm256_srai_epi32(vt, 31) ) ;
  vr   = vt ;
  j0 = 0 ;
  // chunks of 8 elements (second chunk may overlap first chunk)
  for(i0=n7 ; i0<n ; i0+=8){
    vi   = _mm256_loadu_si256((__m256i *) (top+i0)  ) ;   // top[i0:i0+7]
    vj1  = _mm256_loadu_si256((__m256i *) (bot+i0)  ) ;   // bot[i0:i0+7]
#if defined(WITH_FETCH)
    vi1  = _mm256_loadu_si256((__m256i *) (top+i0-1)) ;   // top[i0-1:i0+6]
    vij1 = _mm256_loadu_si256((__m256i *) (bot+i0-1)) ;   // bot[i0-1:i0+6]
#else
    v0   = _mm256_set1_epi32(top[i0-1]) ;                 // top[i0-1]
    vi1  = _mm256_permutevar8x32_epi32(vi, vs) ;          // 0, 0, 1, 2, 3, 4, 5, 6  permutation
    vi1  = _mm256_blend_epi32(vi1, v0, 1) ;               // first element set to top[i0-1]
    v0   = _mm256_set1_epi32(bot[i0-1]) ;                 // bot[i0-1]
    vij1 = _mm256_permutevar8x32_epi32(vj1, vs) ;         // 0, 0, 1, 2, 3, 4, 5, 6  permutation
    vij1 = _mm256_blend_epi32(vij1, v0, 1) ;              // first element set to bot[i0-1]
#endif
    // top[i0+i] - ( top[i0+i-1] + bot[i0+i] - bot[i0+i-1] )
    vt   = _mm256_sub_epi32( vi, _mm256_sub_epi32( _mm256_add_epi32(vi1, vj1) , vij1 ) ) ;
    _mm256_storeu_si256( (__m256i *) (top+j0), vr) ;      // delayed store, result of previous pass
//     vr =  _mm256_xor_si256( _mm256_slli_epi32(vt, 1) , _mm256_srai_epi32(vt, 31) ) ;
    vr = vt ;
    j0 = i0 ;
  }
  _mm256_storeu_si256( (__m256i *) (top+j0), vr) ;        // delayed store, result of previous pass
}
#else

STATIC inline void LorenzoPredictRowJ_inplace(int32_t * restrict top, int32_t * restrict bot, int n){
  int i, i0, j0, n7 = (n & 7) ;
  int32_t t[8], r[8], r0 ;

  if(n < 8) {
    LorenzoPredictRowJ_inplace_07(top, bot, n) ;
    return ;
  }
  r0 = top[0] - bot[0] ;
  n7 = n7 ? n7 : 8 ;
  // first 8 elements
  for(i=0 ; i<8 ; i++) r[i] = top[i] - ( top[i-1] + bot[i] - bot[i-1] ) ;
  j0 = 0 ;
  // chunks of 8 elements (second chunk may overlap first chunk)
  for(i0=n7 ; i0<n ; i0+=8){
    for(i=0 ; i<8 ; i++) t[i] = top[i0+i] - ( top[i0+i-1] + bot[i0+i] - bot[i0+i-1] ) ;
    for(i=0 ; i<8 ; i++) top[j0+i] = r[i] ;
    for(i=0 ; i<8 ; i++) r[i] = t[i] ;
    j0 = i0 ;
  }
  for(i=0 ; i<8 ; i++) top[j0+i] = r[i] ;
  top[0] = r0 ;
}
#endif

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
    LorenzoPredictRowJ_inplace(orig, orig-lnio, ni) ;
//     LorenzoPredictRowJ(orig, orig-lnio, diff, ni) ;   // predict upper row in row pair -> diff
//     mem_cpy_w32(orig, diff, ni) ;                     // copy predicted row back into orig
//     memcpy(orig, diff, sizeof(diff)) ;                // copy predicted row back into orig
    orig -= lnio ;                                    // next row
  }
  LorenzoPredictRow0_inplace(orig, ni) ;
//   LorenzoPredictRow0(orig, diff, ni) ;                // bottom row
//   mem_cpy_w32(orig, diff, ni) ;                       // copy predicted row back into orig
//   memcpy(orig, diff, sizeof(diff)) ;                  // copy predicted row back into orig
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
// void LorenzoUnpredict_c(int32_t * restrict orig, int32_t * restrict diff, int ni, int lnio, int lnid, int nj){
//   int i ;
// 
//   orig[0] = diff[0] ;                                    // restore first point of bottom row
//   for(i=1 ; i<ni ; i++) orig[i] = diff[i] + orig[i-1] ;  // restore bottom row
// 
//   while(--nj > 0){
//     orig += lnio ; 
//     diff += lnid ;
//     orig[0] = diff[0] + orig[0-lnio] ;                   // first point in row (1D prediction)
//     // (original - predicted) + predicted
//     for(i=1 ; i<ni ; i++) orig[i] = diff[i] + (orig[i-1] + orig[i-lnio] - orig[i-1-lnio]) ;
//   }
// }

void LorenzoUnpredict_c(int32_t * restrict orig, int32_t * restrict diff, int ni, int lnio, int lnid, int nj){
  int i ;
  int32_t *top, *bot ;
  int32_t d00, d01, d10, d11 ;
//   d01 d11   unpredict : d11 = d11 + d01 + d10 - d00
//   d00 d10

  d00 = orig[0] = diff[0] ;         // restore first point of bottom row
  for(i=1 ; i<ni ; i++) {           // restore rest of bottom row
    d00 = orig[i] = diff[i] + d00 ;
  }

  while(--nj > 0){
    bot = orig ;
    orig += lnio ; 
    diff += lnid ;
    top = diff ;
    d01 = orig[0] = top[0] + bot[0] ;                // first point in row (1D prediction)
    d00 = bot[0] ;
    for(i=1 ; i<ni ; i++) {
      d10 = bot[i] ;
      d11 = top[i] ;
      d01 = orig[i] = (d11 + d01) + (d10 - d00) ;    // (original - predicted) + predicted
      d00 = d10 ;
    }
  }
}

// void LorenzoUnpredictInplace_c(int32_t * restrict orig, int ni, int lnio, int nj){
//   int i ;
// 
//   for(i=1 ; i<ni ; i++) orig[i] = orig[i] + orig[i-1] ;  // restore bottom row
// 
//   while(--nj > 0){
//     orig += lnio ; 
//     orig[0] = orig[0] + orig[0-lnio] ;                   // first point in row (1D prediction)
//     // (original - predicted) + predicted
//     for(i=1 ; i<ni ; i++) orig[i] = orig[i] + (orig[i-1] + orig[i-lnio] - orig[i-1-lnio]) ;
//   }
// }

void LorenzoUnpredictInplace_c(int32_t *orig, int ni, int lnio, int nj){
  int i ;
  int32_t *top, *bot ;
  int32_t d00, d01, d10, d11 ;
//   d01 d11   unpredict : d11 = d11 + d01 + d10 - d00
//   d00 d10

  d00 = orig[0] ;
  for(i=1 ; i<ni ; i++) { d00 = orig[i] = orig[i] + d00 ; } // restore rest of bottom row

  while(--nj > 0){
    bot   = orig ;
    orig += lnio ; 
    top   = orig ;
    d01 = orig[0] = top[0] + bot[0] ;                // first point in row (1D prediction)
    d00 = bot[0] ;

    for(i=1 ; i<ni ; i++) {
      d10 = bot[i] ;
      d11 = top[i] ;
      d01 = orig[i] = (d11 + d01) + (d10 - d00) ;    // (original - predicted) + predicted
      d00 = d10 ;
    }
  }
}

