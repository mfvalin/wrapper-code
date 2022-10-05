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

// specific optimization option for gcc
#if defined(__GNUC__)
#if ! defined(__INTEL_COMPILER_UPDATE)
#pragma GCC optimize "tree-vectorize"
#endif
#endif

#include <stdint.h>
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
#include <immintrin.h>
#endif

// static int b_conts = 0b01010101 ;

#define STATIC static
#include <misc_operators.h>
#undef STATIC
#define STATIC extern

// 2x2 averaging of 8 or less pairs of points on 2 lines (averaging along line and across lines)
// src1[n*2] : first line
// src2[n*2] : second line
// avg[n]    : result (2x2 average)
// n         : number of pairs of points
static inline void average_2x2_16_I32(int32_t * restrict src1, int32_t * restrict src2, int32_t * restrict avg, int n){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  if(n == 8){
    __m256i vs0, vs1, va1, va2, vb1, vb2, vk2 ;
    va1 = _mm256_loadu_si256( (__m256i *)  src1   ) ;
    vk2 = _mm256_cmpeq_epi32(va1, va1) ;           // -1
    vk2 = _mm256_add_epi32(vk2, vk2) ;             // -2
    va2 = _mm256_loadu_si256( (__m256i *) (src1+8)) ;
    vb1 = _mm256_loadu_si256( (__m256i *)  src2   ) ;
    vb2 = _mm256_loadu_si256( (__m256i *) (src2+8)) ;
    va1 = _mm256_add_epi32(va1, vb1) ;             // add rows terms(0-7)
    va2 = _mm256_add_epi32(va2, vb2) ;             // add rows terma(8-15)
    vs1 = _mm256_hadd_epi32(va1, va2) ;            // add pairs of terms
    vs1 = _mm256_permute4x64_epi64(vs1, 0b11011000) ;
    vs0 = _mm256_sub_epi32(vs1, vk2) ;             // add 2 (subtract -2)
    vs1 = _mm256_srai_epi32(vs1, 31) ;             // -1 or 0 according to sign
    vs0 = _mm256_add_epi32(vs1, vs0) ;             // add to vs0
    vs0 = _mm256_srai_epi32(vs0, 2) ;              // finalize divide by 4 with rounding
    _mm256_storeu_si256( (__m256i *) avg, vs0) ;   // store result
    return ;
  }
#endif
  int i, ii ;
  for(i=0, ii=0 ; ii<n ; ii++, i+=2){
    avg[ii]  = src1[i] + src1[i+1] ;
    avg[ii] += src2[i] + src2[i+1] ;
    avg[ii]  = IDIV4R(avg[ii]) ;
  }
}

static inline void average_2x2_16_F32(float * restrict src1, float * restrict src2, float * restrict avg, int n){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  if(n == 8){
    __m256 vs0, vs1, va1, va2, vb1, vb2, vk2 ;
    vk2 = _mm256_set1_ps(.25f) ;
    va1 = _mm256_loadu_ps(  src1   ) ;
    va2 = _mm256_loadu_ps( (src1+8)) ;
    vb1 = _mm256_loadu_ps(  src2   ) ;
    vb2 = _mm256_loadu_ps( (src2+8)) ;
    va1 = _mm256_add_ps(va1, vb1) ;             // add rows terms(0-7)
    va2 = _mm256_add_ps(va2, vb2) ;             // add rows terma(8-15)
    vs1 = _mm256_hadd_ps(va1, va2) ;            // add pairs of terms
    vs1 = (__m256) _mm256_permute4x64_epi64((__m256i) vs1, 0b11011000) ;
    vs0 = _mm256_mul_ps(vs0, vk2) ;            // divide by 4
    _mm256_storeu_ps(avg, vs0) ;                // store result
    return ;
  }
#endif
  int i, ii ;
  for(i=0, ii=0 ; ii<n ; ii++, i+=2){
    avg[ii]  = src1[i] + src1[i+1] ;
    avg[ii] += src2[i] + src2[i+1] ;
    avg[ii] *= .25f ;
  }
}

// 2x2 averaging of n points on 2 lines (averaging along line and across lines)
// src1[n]      : first line
// src2[n]      : second line
// avg[(n+1)/2] : result
// n            : number of points (may be odd)
STATIC inline void average_2x2_I32(int32_t * restrict src1, int32_t * restrict src2, int32_t * restrict avg, uint32_t n){
  int i, ii ;
  int n2 = n >> 1 ;                         // number of averaged pairs
  int n1 = (n2 & 7) ? (n2 & 7) : 8 ;        // size of first chunk
  int n8 = (n>16) ? 8 : n1 ;                // 8 if more than 16 points

//   if(n1 <= 4) average_2x2_8_I32(src1, src2, avg, n) ;
//   else        average_2x2_16_I32(src1, src2, avg, 8) ;
  average_2x2_16_I32(src1, src2, avg, n8) ;          // first chunk

  for(i=n1+n1, ii=n1 ; ii<n2 ; ii+=8, i+=16){        // 16x2 -> 8x1
    average_2x2_16_I32(src1+i, src2+i, avg+ii, 8) ;  // next chunks (8 pairs)
  }
  if(n & 1) {    // odd number of points in lines, 2 points are averaged instead of 4
    avg[ii] = src1[i] + src2[i] ;
    avg[ii] = IDIV2R(avg[ii]) ;
  }
}

STATIC inline void average_2x2_F32(float * restrict src1, float * restrict src2, float * restrict avg, uint32_t n){
  int i, ii ;
  int n2 = n >> 1 ;                         // number of averaged pairs
  int n1 = (n2 & 7) ? (n2 & 7) : 8 ;        // size of first chunk
  int n8 = (n>16) ? 8 : n1 ;                // 8 if more than 16 points

//   if(n1 <= 4) average_2x2_8_I32(src1, src2, avg, n) ;
//   else        average_2x2_16_I32(src1, src2, avg, 8) ;
  average_2x2_16_F32(src1, src2, avg, n8) ;          // first chunk

  for(i=n1+n1, ii=n1 ; ii<n2 ; ii+=8, i+=16){        // 16x2 -> 8x1
    average_2x2_16_F32(src1+i, src2+i, avg+ii, 8) ;  // next chunks (8 pairs)
  }
  if(n & 1) {    // odd number of points in lines, 2 points are averaged instead of 4
    avg[ii]  = src1[i] + src2[i] ;
    avg[ii] *= .5f ;
  }
}

// average one row
// src1[n]      : row to be averaged
// avg[(n+1)/2] : result
// n            : number of points (may be odd) (MUST BE >= 8)
// NOTE : average_2x2_I32 is used, simulating 2 identical rows
STATIC inline void average_2x1_I32(int32_t * restrict src, int32_t * restrict avg, uint32_t n){
  int i, ii ;
  int n2 = n>>1 ;
  average_2x2_I32(src, src, avg, n) ;
}

STATIC inline void average_2x1_F32(float * restrict src, float * restrict avg, uint32_t n){
  int i, ii ;
  int n2 = n>>1 ;
  average_2x2_F32(src, src, avg, n) ;
}

// 2x2 averaging of a 2 dimensional array
// src[lni       *        nj] : source array
// avg[((ni+1)/2)*((nj+1)/2)] : result
// lni  : storage dimension of rows in src
// ni   : number of points to be averaged in rows (may be odd) (MUST BE >= 8)
// nj   : number of rows to be averaged (may be odd) (MUST BE >= 4)
void average_2x2_2D_I32(int32_t * restrict src, int32_t * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj){
  int i, j ;
  int ni2 = (ni+1)/2 ;
  int lni2 = lni+lni ;
  for(j=0 ; j<nj/2 ; j++){
    average_2x2_I32(src, src+lni, avg, ni) ;
    src += lni2 ;
    avg += ni2 ;
  }
  if(nj & 1) average_2x2_I32(src, src, avg, ni) ;
}

void average_2x2_2D_F32(float * restrict src, float * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj){
  int i, j ;
  int ni2 = (ni+1)/2 ;
  int lni2 = lni+lni ;
  for(j=0 ; j<nj/2 ; j++){
    average_2x2_F32(src, src+lni, avg, ni) ;
    src += lni2 ;
    avg += ni2 ;
  }
  if(nj & 1) average_2x2_F32(src, src, avg, ni) ;
}

// compute half resolution average and full resolution average error

// compute half resolution average and full resolution average error
void avgres_2x2_2D_I32(int32_t * restrict src, int32_t * restrict avg, int32_t * restrict res, int ni, int lni, int nj){
  int i, ii, j ;
  int ni2 = (ni+1)/2 ;
  int lni2 = lni+lni ;
  for(j=0 ; j<nj ; j+=2){
    average_2x2_I32(src, src+lni, avg, ni) ;
    for(i=0, ii=0 ; i<ni-1 ; i+=2, ii++){
      res[i]       = avg[ii] ;
      res[i+1]     = avg[ii] ;
    }
    if(i < ni) res[ni-1] = avg[ni2-1] ;
    for(i=0 ; i<ni ; i++){
      res[i+lni]   = src[i+lni] - res[i] ;
      res[i]       = src[i]     - res[i] ;
    }
    src += lni2 ;
    res += lni2 ;
    avg += (ni+1)/2 ;
  }
  if(j < nj){   // odd number of rows
    average_2x2_I32(src, src, avg, ni) ;
    for(i=0, ii=0 ; i<ni-1 ; i+=2, ii++){
      res[i]       = src[i]       - avg[ii] ;
      res[i+1]     = src[i+1]     - avg[ii] ;
    }
    if(i < ni) {   // odd number of points in last row
      res[i]  = 0 ;
    }
  }
}

// expand a row along i
// SIMD version WILL NOT WORK if n < 16, in that case the non SIMD code will be used
STATIC void expand_2x1_row_along_i(int32_t * restrict dst, int32_t * restrict avg, int n){
  int i, ii ;
  int n2 = n>>1 ;
  int t1[n2], t2[n2] ;

  dst[0] = 5*avg[0] - avg[1] ; dst[0] = IDIV4R(dst[0]) ;
  dst[n2+n2-1] = 5*avg[n2-1] - avg[n2-2] ; dst[n2+n2-1] = IDIV4R(dst[n2+n2-1]) ;
  if(n & 1) dst[n2+n2] = avg[n2] ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  if(n >= 16){                          // the SIMD version needs at least 16 elements
    int nm2 = n2-1, i2 = 1 ;
    int n7 = (nm2 & 7) ? (nm2 & 7) : 8 ;
    __m256i va1, va2, va3, vb1, vb2, vb3, v31, v13, q31, q13, vlo, vhi, vk2, v07, v8f ;
    for(i=0 ; i<nm2 ; i+= n7, n7 = 8) {                   // expand 8 averages -> 16x1 values
      va1 = _mm256_loadu_si256( (__m256i *) (avg + i   ) );   // avg[ii  :ii+7]
      vb1 = _mm256_loadu_si256( (__m256i *) (avg + i +1) );   // avg[ii+1:ii+8]
      vk2 = _mm256_cmpeq_epi32(va1, va1) ;                // -1
      vk2 = _mm256_add_epi32(vk2, vk2) ;                  // -2
      va2 = _mm256_add_epi32(va1, va1) ;                  // va1 * 2
      vb2 = _mm256_add_epi32(vb1, vb1) ;                  // vb1 * 2
      va3 = _mm256_add_epi32(va2, va1) ;                  // va1 * 3
      vb3 = _mm256_add_epi32(vb2, vb1) ;                  // vb1 * 3
      v31 = _mm256_add_epi32(va3, vb1) ;                  // 3*row0[i] + 1*row1[i]
      v13 = _mm256_add_epi32(va1, vb3) ;                  // 1*row0[i] + 3*row1[i]
      q31 = _mm256_srai_epi32(v31, 31) ;                  // -1 or 0 according to sign
      q13 = _mm256_srai_epi32(v13, 31) ;                  // -1 or 0 according to sign
      q31 = _mm256_sub_epi32(q31, vk2) ;                  // q31 + 2
      q13 = _mm256_sub_epi32(q13, vk2) ;                  // q13 + 2
      q31 = _mm256_add_epi32(v31, q31) ;                  // + v31
      q13 = _mm256_add_epi32(v13, q13) ;                  // + v13
      q31 = _mm256_srai_epi32(q31, 2) ;                   // IDIV4R(v31)  (dst[1, 3, 5, ...])
      q13 = _mm256_srai_epi32(q13, 2) ;                   // IDIV4R(v13)  (dst[2, 4, 6, ...])
      vlo = _mm256_unpacklo_epi32(q31, q13) ;             // 0 1 2 3 8 9 A B
      vhi = _mm256_unpackhi_epi32(q31, q13) ;             // 4 5 6 7 C D E F
      v07 = _mm256_permute2x128_si256(vlo, vhi, 0x20) ;   // 0 1 2 3 4 5 6 7 (dst[1, 2, 3, ...])
      v8f = _mm256_permute2x128_si256(vlo, vhi, 0x31) ;   // 8 9 A B C D E F (dst(9, A, B, ...])
      _mm256_storeu_si256( (__m256i *) (dst + i2    ), v07 );  // dst[i  :i+ 7]
      _mm256_storeu_si256( (__m256i *) (dst + i2 + 8), v8f );  // dst[i+8:i+15]
      i2 = i2 + n7 + n7 ;
    }
    return ;
  }
#endif
  for(ii=0 ; ii<n2-1 ; ii++) {
    t1[ii] = 3*avg[ii] + 1*avg[ii+1] ;
    t1[ii] = IDIV4R(t1[ii]) ;
    t2[ii] = 1*avg[ii] + 3*avg[ii+1] ;
    t2[ii] = IDIV4R(t2[ii]) ;
  }
  for(i=1, ii=0 ; i<n-2 ; i+=2, ii++){
    dst[i  ] = t1[ii] ; dst[i+1] = t2[ii] ;
  }
}

// bottom row
STATIC void expand_2x2_row_0(int32_t * restrict row0, int32_t * restrict row1, int32_t * src0, int n){
  int i ;
  int ni2 = (n+1)>>1 ;
  int32_t a51[ni2] ;
  for(i=0 ; i<ni2 ; i++) {
    a51[i] = 5*row0[i] - 1*row1[i] ;  // unaverage along j to restore src0
    a51[i] = IDIV4R(a51[i]) ;
    expand_2x1_row_along_i(src0, a51, n) ;  // restore src0 along i
  }
}

// pairs of middle rows
STATIC void expand_2x2_2_rows(int32_t * restrict row0, int32_t * restrict row1, int32_t * src0, int32_t * src1, int n){
  int i ;
  int ni2 = (n+1)>>1 ;
  int32_t a31[ni2], a13[ni2] ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  int n7 = (ni2 & 7) ? (ni2 & 7) : 8 ;
  __m256i va1, va2, va3, vb1, vb2, vb3, v31, v13, q31, q13, vk2 ;
  for(i=0 ; i<ni2 ; i+= n7, n7 = 8) {
    va1 = _mm256_loadu_si256( (__m256i *) (row0 + i) );
    vb1 = _mm256_loadu_si256( (__m256i *) (row1 + i) );
    vk2 = _mm256_cmpeq_epi32(va1, va1) ;    // -1
    vk2 = _mm256_add_epi32(vk2, vk2) ;      // -2
    va2 = _mm256_add_epi32(va1, va1) ;      // va1 * 2
    vb2 = _mm256_add_epi32(vb1, vb1) ;
    va3 = _mm256_add_epi32(va2, va1) ;      // va1 * 3
    vb3 = _mm256_add_epi32(vb2, vb1) ;
    v31 = _mm256_add_epi32(va3, vb1) ;      // 3*row0[i] + 1*row1[i]
    v13 = _mm256_add_epi32(va1, vb3) ;      // 1*row0[i] + 3*row1[i]
    q31 = _mm256_srai_epi32(v31, 31) ;      // -1 or 0 accroding to sign
    q13 = _mm256_srai_epi32(v13, 31) ;      // -1 or 0 accroding to sign
    q31 = _mm256_sub_epi32(q31, vk2) ;      // q31 + 2
    q13 = _mm256_sub_epi32(q13, vk2) ;      // q13 + 2
    q31 = _mm256_add_epi32(v31, q31) ;      // + v31
    q13 = _mm256_add_epi32(v13, q13) ;      // + v13
    q31 = _mm256_srai_epi32(q31, 2) ;       // IDIV4R(v31)
    q13 = _mm256_srai_epi32(q13, 2) ;       // IDIV4R(v13)
    _mm256_storeu_si256( (__m256i *) (a31 + i), q31 );
    _mm256_storeu_si256( (__m256i *) (a13 + i), q13 );
  }
#else
  for(i=0 ; i<ni2 ; i++) {
    a31[i] = 3*row0[i] + 1*row1[i] ;  // unaverage along j to restore src0
    a31[i] = IDIV4R(a31[i]) ;
    a13[i] = 1*row0[i] + 3*row1[i] ;  // unaverage along j to restore src1
    a13[i] = IDIV4R(a13[i]) ;
  }
#endif
  expand_2x1_row_along_i(src0, a31, n) ;  // restore src0 along i
  expand_2x1_row_along_i(src1, a13, n) ;  // restore src1 along i
}

// top row if even number of rows
STATIC void expand_2x2_row_n(int32_t * restrict row0, int32_t * restrict row1, int32_t * src0, int n){
  int i ;
  int ni2 = (n+1)>>1 ;
  int32_t a15[ni2] ;
  for(i=0 ; i<ni2 ; i++) {
    a15[i] = 5*row1[i] - 1*row0[i] ;  // unaverage along j to restore src0
    a15[i] = IDIV4R(a15[i]) ;
    expand_2x1_row_along_i(src0, a15, n) ;  // restore src0 along i
  }
}

// restore a previously averaged array using linear interpolation/extrapolation
// dst[lni       *        nj] : result
// avg[((ni+1)/2)*((nj+1)/2)] : previously averaged array
// lni  : storage dimension of rows in dst
// ni   : number of averaged points in rows (may be odd)
// nj   : number of averaged rows (may be odd)
void expand_2x2_2D_I32(int32_t * restrict dst, int32_t * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj){
  int i, j ;
  int ni2 = (ni+1)/2 ;
  int lni2 = lni+lni ;
  int32_t t[ni] ;

  // first row
  expand_2x2_row_0(avg, avg+ni2, dst, ni) ;
  dst += lni ;
  for(j=1 ; j<nj/2 ; j++){        // pairs of rows
    expand_2x2_2_rows(avg, avg+ni2, dst, dst+lni, ni) ;
    avg += ni2 ;
    dst += lni2 ;
  }
  // last row
  expand_2x2_row_n(avg-ni2, avg, dst, ni) ;
  if(nj & 1) expand_2x1_row_along_i(dst+lni, avg+ni2, ni) ;
}

