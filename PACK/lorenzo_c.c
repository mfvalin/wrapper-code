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
#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
#include <immintrin.h>
#endif

// target is AVX2 256 long SIMD
#define VL 8

// 2D lorenzo prediction (32 bit signed integers)
// orig : input : original values (32 bit signed integers)
// diff : output : original value - predicted value (using 2D Lorenzo predictor)
// ni   : number of useful points in row
// lnio : row storage dimension for orig
// lnid : row storage dimension for diff
// nj   : number of rows
// the SIMD version is 2-4x faster then the non SIMD version
void LorenzoPredict2D(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj){
  int i, i0, ii0 ;
#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
    __m256i vi, vi1, vj1, vij1, t ;
#endif
    int tmp[VL], j ;

  diff[0] = orig[0] ;  // first point is left untouched
  // predicted value = orig(i-1,1) (lower row)
  if(ni < 9) {
    for(i=1 ; i<ni ; i++) diff[i] = orig[i] - orig[i-1] ;
  }else{
    for(ii0 = 1 ; ii0 < ni ; ii0 += VL){
      i0 = (ii0 > ni-VL) ? (ni-VL) : ii0 ;
#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
      vi   = (__m256i) _mm256_loadu_ps((float *) orig+i0) ;
      vi1  = (__m256i) _mm256_loadu_ps((float *) orig+i0-1) ;
      t    = _mm256_sub_epi32(vi, vi1) ;
      _mm256_storeu_ps( (float *) diff+i0,  (__m256) t ) ;
#else
      for(i=0 ; i<VL ; i++) diff[i0+i] = orig[i0+i] - orig[i0+i-1] ;
#endif
    }
  }

  while(--nj > 0){
    diff += lnid ; 
    orig += lnio ;
    diff[0] = orig[0] - orig[0-lnio] ;    // first point in row, predicted value = orig(i,j-1)
    // predicted value = orig(i,j) + orig(i,j-1) - orig(i-1,j-1)
    if(ni < 9) {
      for(i=1 ; i<ni ; i++) diff[i] = orig[i] - (orig[i-1] + orig[i-lnio] - orig[i-1-lnio]) ;
    }else{
      for(ii0 = 1 ; ii0 < ni ; ii0 += VL){
        i0 = (ii0 > ni-VL) ? (ni-VL) : ii0 ;
#if defined(WITH_SIMD) && defined(__AVX2__) && defined(__x86_64__)
        vi   = (__m256i) _mm256_loadu_ps((float *) orig+i0) ;
        vi1  = (__m256i) _mm256_loadu_ps((float *) orig+i0-1) ;
        vj1  = (__m256i) _mm256_loadu_ps((float *) orig+i0-lnio) ;
        vij1 = (__m256i) _mm256_loadu_ps((float *) orig+i0-lnio-1) ;
        _mm256_storeu_ps( (float *) diff+i0, (__m256) _mm256_sub_epi32( vi, _mm256_sub_epi32( _mm256_add_epi32(vi1, vj1) , vij1 ) ) );
#else
        for(i=0 ; i<VL ; i++) diff[i0+i] = orig[i0+i] - (orig[i0+i-1] + orig[i0+i-lnio] - orig[i0+i-1-lnio]) ;
#endif
      }
    }
  }
}

// restore ogiginal from 2D lorenzo prediction (32 bit signed integers)
// diff : input : original value - predicted value (32 bit signed integers)
// orig : output : restored original values from predicted differences
// ni   : number of useful points in row
// lnio : row storage dimension for orig
// lnid : row storage dimension for diff
// nj   : number of rows
// NOTE : no SIMD version exists, as the process is fully recursive
void LorenzoUnpredict2D(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj){
  int i, i0, ii0 ;

  orig[0] = diff[0] ;                                    // restore first point
  for(i=1 ; i<ni ; i++) orig[i] = diff[i] + orig[i-1] ;  // restore lower row

  while(--nj > 0){
    orig += lnio ; 
    diff += lnid ;
    orig[0] = diff[0] + orig[0-lnio] ;
    // (original - predicted) + predicted
    for(i=1 ; i<ni ; i++) orig[i] = diff[i] + (orig[i-1] + orig[i-lnio] - orig[i-1-lnio]) ;
  }
}
