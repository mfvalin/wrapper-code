// 
// Copyright (C) 2022  Environnement Canada
// 
// This is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
// 
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// Author:
//     M. Valin,   Recherche en Prevision Numerique, august 2022
// 
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <misc_types.h>

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
#include <emmintrin.h>
#include <immintrin.h>

// shuffle table to pack the upper 24 bits from 4 32 bit words into 3 32 bit words
static uint8_t _32_from_upper24[] = { 7, 1, 2, 3, 10, 11, 5, 6, 13, 14, 15, 9, 128, 128, 128, 128 } ;
// shuffle table to pack the lower 24 bits from 4 32 bit words into 3 32 bit words
static uint8_t _32_from_lower24[] = { 6, 0, 1, 2,9, 10, 4, 5, 12, 13, 14, 8, 128, 128, 128, 128 } ;
// shuffle table to pack the lower 16 bits from 4 32 bit words into 2 32 bit words
static uint8_t _32_from_lower16[] = { 4, 5, 0, 1, 12, 13, 8, 9, 128, 128, 128, 128, 128, 128, 128, 128 } ;

// shuffle table to extract 4 24 bit items from 3 32 bit words into the upper part of a 32 bit word
static uint8_t _32_to_upper24[] = {128, 1, 2, 3, 128, 6, 7, 0, 128, 11, 4, 5, 128, 8, 9, 10 } ;
// shuffle table to extract 4 24 bit items from 3 32 bit words into the lower part of a 32 bit word
static uint8_t _32_to_lower24[] = { 1, 2, 3, 128, 6, 7, 0, 128, 11, 4, 5, 128, 8, 9, 10, 128 } ;
// shuffle table to extract 8 16 bit items from 4 32 bit words into the lower part of a 32 bit word
static uint8_t _32_to_lower16[] = {  2,  3, 128, 128,  0,  1, 128, 128,  6,  7, 128, 128,  4,  5, 128, 128,    // lower 2 words
                                    10, 11, 128, 128,  8,  9, 128, 128, 14, 15, 128, 128, 12, 13, 128, 128 } ; // upper 2 words

// mask for masked store of first 3 elements in vector
static uint32_t mask24[4] = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0 } ;

#endif

// pack upper 16 bits from 4 32 bit words into 2 32 bit unsigned integers
// normally used to pack floats 
// sign / 8 bit exponent / 7 bit mantissa : 16 bit "brain float"
void fp32_bf16(void *f32, uint32_t *u16, int32_t n){
  int32_t i ;
  uint32_t *u32 = f32 ;
  uint32_t round = 0x8000 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vr, v32 ;
  __m128i v16a, v16b, v16c, vx ;
  vr = _mm256_set1_epi32(round) ;  // rounding term to be added before clipping mantissa
  vx = _mm_loadu_si128((__m128i*) _32_from_lower16) ;
#endif
  i=0 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  for(    ; i < n-15 ; i += 16){
    v32  = _mm256_loadu_si256((__m256i*)  u32) ;
    v32  = _mm256_add_epi32(v32, vr) ;
    v32  = _mm256_srli_epi32(v32,16) ;
    v16a = _mm256_extracti128_si256(v32, 0) ;
    v16b = _mm256_extracti128_si256(v32, 1) ;
    v16a = _mm_shuffle_epi8(v16a, vx) ;
    v16b = _mm_shuffle_epi8(v16b, vx) ;
    v16c = _mm_unpacklo_epi64(v16a, v16b) ;
    _mm_storeu_si128((__m128i *) u16  , v16c) ;
//     _mm_storeu_si64(u16  , v16a) ;
//     _mm_storeu_si64(u16+2, v16b) ;
//     _mm_storeu_si128((__m128i*) u16, _mm_packus_epi32(v16a, v16b)) ;
    v32  = _mm256_loadu_si256((__m256i*)  (u32+8)) ;
    v32  = _mm256_add_epi32(v32, vr) ;
    v32  = _mm256_srli_epi32(v32,16) ;
    v16a = _mm256_extracti128_si256(v32, 0) ;
    v16b = _mm256_extracti128_si256(v32, 1) ;
    v16a = _mm_shuffle_epi8(v16a, vx) ;
    v16b = _mm_shuffle_epi8(v16b, vx) ;
    v16c = _mm_unpacklo_epi64(v16a, v16b) ;
    _mm_storeu_si128((__m128i *) (u16+4), v16c) ;
//     _mm_storeu_si64(u16+4, v16a) ;
//     _mm_storeu_si64(u16+6, v16b) ;
//     _mm_storeu_si128((__m128i*) (u16+4), _mm_packus_epi32(v16a, v16b)) ;
    u32 += 16 ;
    u16 +=  8 ;
  }
#endif
  // process leftovers
  for(    ; i < n-1 ; i += 2){
    u16[0]  =  (u32[0]+round) & (~0xFFFF) ;  // upper 16 bits of u32[0] into upper 16 bits of u16[0]
    u16[0] |= ((u32[1]+round) >> 16) ;       // upper 16 bits of u32[1] into lower 16 bits of u16[0]
    u32 += 2 ;
    u16++ ;
  }
  if(i < n) u16[0]  =  (u32[0]+round) & (~0xFFFF) ;
}

// pack upper 16 bits from 4 32 bit words into 2 32 bit unsigned integers
// normally used to pack floats 
// 1 bit : sign / nexp bits : exponent / (16 - nexp - 1) bits : mantissa
void fp32_nf16(void *f32, uint32_t *u16, int32_t n, int32_t nexp){
  int32_t i ;
  FloatUint *u32 = f32 ;
  uint32_t round, top, fexp, sign ;
  uint32_t s8 ;
  FloatUint nfac, z ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vi, v8, vs, vt, vr, vn, v1, v32, vm ;
  __m256  vf ;
  __m128i v16a, v16b, v16c, vx ;
#endif
  i=0 ;
  s8 = (8 - nexp) ;                           // shift count for exp / mantissa
  round = (0x8000 >> s8) ;                    // rounding term to be added before clipping mantissa
  fexp = (0x7F >> s8) ;                       // new offset for exponent
  nfac.u = (fexp  << 23) ;                    // normalisation factor
  top = nfac.u | 0x7FFFFF ;                   // maximum allowable value after normalisation
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  vr = _mm256_set1_epi32(round) ;
  vs = _mm256_set1_epi32(s8) ;
  vt = _mm256_set1_epi32(top) ;
  vn = _mm256_set1_epi32(nfac.u) ;
  vm = _mm256_set1_epi32(0x80000000u) ;    // sign mask
  vx = _mm_loadu_si128((__m128i*) _32_from_lower16) ;
  for(    ; i < n-15 ; i += 16){
    vf = _mm256_loadu_ps((float *) u32 ) ;          // get float
    vi = (__m256i) _mm256_mul_ps(vf, (__m256) vn) ; // normalize
    v1 = _mm256_and_si256(vm, vi) ;                 // get sign
    vi = _mm256_andnot_si256(vm, vi) ;              // absolute value
    vi = _mm256_add_epi32(vr, vi) ;                 // add rounding factor
    vi = _mm256_min_epi32(vt, vi) ;                 // "infinity" clip
    vi = _mm256_sllv_epi32(vi, vs) ;                // shift
    vi = _mm256_or_si256(v1, vi) ;                  // restore sign
    vi  = _mm256_srli_epi32(vi,16) ;                // shift right by 16 bits
    v16a = _mm256_extracti128_si256(vi, 0) ;
    v16b = _mm256_extracti128_si256(vi, 1) ;
    v16a = _mm_shuffle_epi8(v16a, vx) ;
    v16b = _mm_shuffle_epi8(v16b, vx) ;
    v16c = _mm_unpacklo_epi64(v16a, v16b) ;
    _mm_storeu_si128((__m128i *) u16, v16c) ;

//     _mm_storeu_si128((__m128i*) u16, _mm_packus_epi32(v16a, v16b)) ;
    vf = _mm256_loadu_ps((float *) (u32+8) ) ;          // get float
    vi = (__m256i) _mm256_mul_ps(vf, (__m256) vn) ; // normalize
    v1 = _mm256_and_si256(vm, vi) ;                 // get sign
    vi = _mm256_andnot_si256(vm, vi) ;              // absolute value
    vi = _mm256_add_epi32(vr, vi) ;                 // add rounding factor
    vi = _mm256_min_epi32(vt, vi) ;                 // "infinity" clip
    vi = _mm256_sllv_epi32(vi, vs) ;                // shift
    vi = _mm256_or_si256(v1, vi) ;                  // restore sign
    vi  = _mm256_srli_epi32(vi,16) ;                // shift right by 16 bits
    v16a = _mm256_extracti128_si256(vi, 0) ;
    v16b = _mm256_extracti128_si256(vi, 1) ;
    v16a = _mm_shuffle_epi8(v16a, vx) ;
    v16b = _mm_shuffle_epi8(v16b, vx) ;
    v16c = _mm_unpacklo_epi64(v16a, v16b) ;
    _mm_storeu_si128((__m128i *) (u16+4), v16c) ;
//     _mm_storeu_si128((__m128i*) (u16+4), _mm_packus_epi32(v16a, v16b)) ;
    u32 += 16 ;
    u16 +=  8 ;
  }
#endif
  // process leftovers
  for(    ; i < n-1 ; i += 2){
    z.f = u32[0].f ;
    sign = z.u & 0x80000000u ;           // keep sign
    z.u &= 0x7FFFFFFFu ;                 // absolute value
    z.f *= nfac.f ;                      // normalize
    z.u += round ;                       // add rounding
    z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
    z.u <<= s8 ;                         // shift
    z.u |= sign ;                        // restore sign
    u16[0]  =  z.u & 0xFFFF0000u ;       // upper 16 bits to upper 16 bits of u16
    z.f = u32[1].f ;
    sign = z.u & 0x80000000u ;           // keep sign
    z.u &= 0x7FFFFFFFu ;                 // absolute value
    z.f *= nfac.f ;                      // normalize
    z.u += round ;                       // add rounding
    z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
    z.u <<= s8 ;                         // shift
    z.u |= sign ;                        // restore sign
    u16[0] |= (z.u >> 16) ;              // upper 16 bits to lower 16 bits of u16
    u32 += 2 ;
    u16++ ;
  }
  if(i < n) {
    z.f = u32[0].f ;
    sign = z.u & 0x80000000u ;           // keep sign
    z.u &= 0x7FFFFFFFu ;                 // absolute value
    z.f *= nfac.f ;                      // normalize
    z.u += round ;                       // add rounding
    z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
    z.u <<= s8 ;                         // shift
    z.u |= sign ;                        // restore sign
    u16[0]  =  z.u & 0xFFFF0000u ;       // upper 16 bits to upper 16 bits of u16
  }
}

// pack upper 24 bits from 4 32 bit words into 3 32 bit unsigned integers
// normally used to pack floats 
// exponent field is kept at 8 bits like 16bits "brain float"
void fp32_bf24(void *f32, uint32_t *u24, int32_t n){
  int32_t i ;
  uint32_t *u32 = f32 ;
  uint32_t n3 = (n & 3) ;
  uint32_t round = 0x80 ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32a, v24a, v32b, v24b, vm, vr ;
  vr = _mm_set1_epi32(round) ;
  vx = _mm_loadu_si128((__m128i*) _32_from_upper24) ;
  vm = _mm_loadu_si128((__m128i*) mask24) ;
  for(i=0 ; i < n-8 ; i+=8){   // 8 words into 6 words
    v32a = _mm_loadu_si128((__m128i*)  u32) ;
    v32b = _mm_loadu_si128((__m128i*) (u32+4)) ;
    v32a = _mm_add_epi32(v32a, vr) ;                  // add rounding term
    v32b = _mm_add_epi32(v32b, vr) ;                  // add rounding term
    v24a = _mm_shuffle_epi8(v32a, vx) ;
    v24b = _mm_shuffle_epi8(v32b, vx) ;
    _mm_storeu_si128((__m128i*)  u24,    v24a) ;      // plenty of extra room
    _mm_storeu_si128((__m128i*) (u24+3), v24b) ;      // there is enough room for one extra word
//     _mm_storeu_si64((__m128i*) (u24+3), v24b) ;
//     v24b = _mm_bsrli_si128 (v24b, 8) ;
//     _mm_storeu_si32((__m128i*) (u24+5), v24b) ;
    u32 += 8 ;
    u24 += 6 ;
  }
  while(i < n-3){   // 4 words into 3 words
    v32a = _mm_loadu_si128((__m128i*) u32) ;
    v32a = _mm_add_epi32(v32a, vr) ;                  // add rounding term
    v24a = _mm_shuffle_epi8(v32a, vx) ;
    _mm_maskstore_epi32((int *) u24, vm, v24a) ;  // mask used to store only 3 words
    u32 += 4 ;
    u24 += 3 ;
    i   += 4 ;
  }
#else
  for(i=0 ; i<n-3 ; i+=4){
    u24[0] = (u32[0]+round) & (~0xFF) ;             // upper 24 bits of u32[0]
    u24[0] |= (u32[1]+round) >> 24 ;                // upper 8 bits of u32[1]
    u24[1]  = (u32[1]+round & (~0xFF)) << 8 ;     // lower 16 of upper 24 bits of u32[1]
    u24[1] |= ((u32[2]+round) >> 16) ;              // upper 16 of upper 24 bits of u32[2]
    u24[2]  = (u32[2]+round & (~0xFF)) << 16 ;    // lower 8 of upper 24 bits of u32[2]
    u24[2] |= ((u32[3]+round) >> 8) ;               // upper 24 bits of u32[3]
    u32 += 4 ;
    u24 += 3 ;
  }
#endif
  // process leftovers
  if(n3 > 0) u24[0] = u32[0]+round & (~0xFF) ;    // upper 24 bits of u32[0]
  if(n3 > 1) {
    u24[0] |= (u32[1]+round) >> 24 ;                // upper 8 bits of u32[1]
    u24[1]  = (u32[1]+round & (~0xFF)) << 8 ;     // lower 16 of upper 24 bits of u32[1]
  }
  if(n3 > 2) {
    u24[1] |= ((u32[2]+round) >> 16) ;              // upper 16 of upper 24 bits of u32[2]
    u24[2]  = (u32[2]+round & (~0xFF)) << 16 ;    // lower 8 of upper 24 bits of u32[2]
  }
// printf("fp32_bf24 leftovers = %d %8.8x %8.8x %8.8x\n", n3, u24[0], u24[1], u24[2]);
}

// pack lower 24 bits from 4 32 bit unsigned integers into 3 32 bit unsigned integers
void u32_u24(uint32_t *u32, uint32_t *u24, int32_t n){
  int i ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32a, v24a, v32b, v24b, vm ;
  vx = _mm_loadu_si128((__m128i*) _32_from_lower24) ;
  vm = _mm_loadu_si128((__m128i*) mask24) ;
  for(i=0 ; i<n-8 ; i+=8){
    v32a = _mm_loadu_si128((__m128i*)  u32) ;
    v32b = _mm_loadu_si128((__m128i*) (u32+4)) ;
    v24a = _mm_shuffle_epi8(v32a, vx) ;
    v24b = _mm_shuffle_epi8(v32b, vx) ;
    _mm_storeu_si128((__m128i*)  u24,    v24a) ;      // plenty of extra room
    _mm_storeu_si128((__m128i*) (u24+3), v24b) ;      // plenty of extra room
    u32 += 8 ;
    u24 += 6 ;
  }
  while(i < n-3){   // 4 words into 3 words
    v32a = _mm_loadu_si128((__m128i*) u32) ;
    v24a = _mm_shuffle_epi8(v32a, vx) ;
    _mm_maskstore_epi32((int *) u24, vm, v24a) ;  // mask used to store only 3 words
    u32 += 4 ;
    u24 += 3 ;
    i += 4 ;
  }
#else
  for(i=0 ; i<n-3 ; i+=4){
    u24[0] = u32[0] << 8 ;                   // lower 24 bits of u32[0]
    u24[0] |= ((u32[1] >> 16) & 0xFF) ;      // upper  8 bits of lower 24 bits of u32[1]
    u24[1]  = u32[1] << 16 ;                 // lower 16 bits of lower 24 bits of u32[1]
    u24[1] |= ((u32[2] >> 8) & 0xFFFF) ;     // upper 16 bits of lower 24 bits of u32[2]
    u24[2]  = u32[2] << 24 ;                 // lower  8 bits of lower 24 bits of u32[2]
    u24[2] |= (u32[3] & 0xFFFFFF) ;          // lower 24 bits of u32[3]
    u32 += 4 ;
    u24 += 3 ;
  }
#endif
  // process leftovers
  if(n3 > 0) u24[0] = u32[0] << 8 ;          // lower 24 bits of u32[0]
  if(n3 > 1) {
    u24[0] |= ((u32[1] >> 16) & 0xFF) ;      // upper  8 bits of lower 24 bits of u32[1]
    u24[1]  = u32[1] << 16 ;                 // lower 16 bits of lower 24 bits of u32[1]
  }
  if(n3 > 2) {
    u24[1] |= ((u32[2] >> 8) & 0xFFFF) ;     // upper 16 bits of lower 24 bits of u32[2]
    u24[2]  = u32[2] << 24 ;                 // lower  8 bits of lower 24 bits of u32[2]
  }
}

// pack lower 24 bits from 4 32 bit signed integers into 3 32 bit unsigned integers
void i32_u24(int32_t *i32, uint32_t *u24, int32_t n){
  u32_u24((uint32_t *) i32, u24, n) ;
}

// restore the upper 16 bits of 4 32 bit words from 16 bit tokens in 32 bit words
// normally used to restore floats
void bf16_fp32(void *f32, uint32_t *u16, int32_t n){
  int i ;
  uint32_t *u32 = f32 ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i v32 ;
  __m128i v16 ;
#endif
  i=0 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  for(    ; i < n-15 ; i += 16){
    v16 = _mm_loadu_si128((__m128i *) u16) ;
    v32 = _mm256_cvtepu16_epi32(v16) ;
    v32 = _mm256_slli_epi32(v32, 16) ;
    v32 = _mm256_shuffle_epi32(v32, 0xB1) ;
    _mm256_storeu_si256((__m256i *) u32, v32) ;
    v16 = _mm_loadu_si128((__m128i *) (u16+4)) ;
    v32 = _mm256_cvtepu16_epi32(v16) ;
    v32 = _mm256_slli_epi32(v32, 16) ;
    v32 = _mm256_shuffle_epi32(v32, 0xB1) ;
    _mm256_storeu_si256((__m256i *) (u32+8), v32) ;
    u32 += 16 ;
    u16 +=  8 ;
  }
#endif
  // process leftovers
  for(    ; i < n-1 ; i += 2){
    u32[0] = u16[0] & (~0xFFFF) ;
    u32[1] = u16[0] << 16 ;
    u32 += 2 ;
    u16 ++ ;
  }
  if(i == n-1) u32[0] = u16[0] & (~0xFFFF) ;
}

// restore the upper 16 bits of 4 32 bit words from 16 bit tokens in 2 32 bit words
// normally used to restore floats with a variable exponent width
void nf16_fp32(void *f32, uint32_t *u16, int32_t n, int nexp){
  int32_t i ;
  FloatUint *u32 = f32 ;
  uint32_t top, fexp, sign ;
  uint32_t s8 ;
  FloatUint nfac, z ;
//   bf16_fp32(f32, u16, n) ;
//   return ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vi, v8, vs, vt, vn, v1, v32, vm ;
  __m256  vf ;
  __m128i v16 ;
#endif
  i = 0 ;
  s8 = (8 - nexp) ;                           // shift count for exp / mantissa
  fexp = (0x7F >> s8) ;                       // new offset for exponent
  nfac.u = (fexp  << 23) ;                    // normalisation factor
//   top = nfac.u | 0x7FFFFF ;                   // maximum allowable value after normalisation
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  vs = _mm256_set1_epi32(s8) ;
  vt = _mm256_set1_epi32(top) ;
  vn = _mm256_set1_epi32(nfac.u) ;
  vm = _mm256_set1_epi32(0x80000000u) ;    // sign mask
  for(    ; i < n-7 ; i += 8){
    v16 = _mm_loadu_si128((__m128i *) u16) ;
    v32 = _mm256_cvtepu16_epi32(v16) ;
    v32 = _mm256_slli_epi32(v32, 16) ;
    v32 = _mm256_shuffle_epi32(v32, 0xB1) ;
    v32 = _mm256_srlv_epi32(v32, vs) ;
    vf  = _mm256_mul_ps((__m256) v32, (__m256) vn) ;
    _mm256_storeu_ps((float *) u32, vf) ;
    u32 +=  8 ;
    u16 +=  4 ;
  }
#endif
  // process leftovers
  for(    ; i < n-1 ; i += 2){
    z.u = u16[0] & 0xFFFF0000u ;     // upper 16 bits of u16
    sign = z.u & 0x80000000u ;       // get sign
    z.u &= 0x7FFFFFFFu ;             // abs value
    z.u >>= s8 ;                     // shift back to proper position
    z.u |= sign ;                    // re apply sign
    u32[0].f = z.f * nfac.f ;        // normalize

    z.u = u16[0] << 16 ;             // lower 16 bits of u16
    sign = z.u & 0x80000000u ;       // get sign
    z.u &= 0x7FFFFFFFu ;             // abs value
    z.u >>= s8 ;                     // shift back to proper position
    z.u |= sign ;                    // re apply sign
    u32[1].f = z.f * nfac.f ;        // normalize
    u32 += 2 ;
    u16 ++ ;
  }
  if(i == n-1) {
    z.u = u16[0] & 0xFFFF0000u ;     // upper 16 bits of u16
    sign = z.u & 0x80000000u ;       // get sign
    z.u &= 0x7FFFFFFFu ;             // abs value
    z.u >>= s8 ;                     // shift back to proper position
    z.u |= sign ;                    // re apply sign
    u32[0].f = z.f * nfac.f ;        // normalize
  }
}

// restore the upper 24 bits of 4 32 bit words from 4 24 bit tokens in 3 32 bit words
// normally used to restore floats
void bf24_fp32(void *f32, uint32_t *u24, int32_t n){
  int i ;
  uint32_t *u32 = f32 ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32a, v24a, v32b, v24b, vm ;
  vx = _mm_loadu_si128((__m128i*) _32_to_upper24) ;
  for(i=0 ; i<n-8 ; i+=8){
    v24a = _mm_loadu_si128((__m128i*)  u24) ;
    v24b = _mm_loadu_si128((__m128i*) (u24+3)) ;
    v32a = _mm_shuffle_epi8(v24a, vx) ;
    v32b = _mm_shuffle_epi8(v24b, vx) ;
    _mm_storeu_si128((__m128i*)  u32,    v32a) ;
    _mm_storeu_si128((__m128i*) (u32+4), v32b) ;
    u32 += 8 ;
    u24 += 6 ;
  }
  while(i < n-3){   // 4 words into 3 words
    v24a = _mm_loadu_si128((__m128i*) u24) ;
    v32a = _mm_shuffle_epi8(v24a, vx) ;
    _mm_storeu_si128((__m128i*) u32, v32a) ;
    u32 += 4 ;
    u24 += 3 ;
    i += 4 ;
  }
#else
  for(i=0 ; i<n-3 ; i+=4){
    u32[0]  = u24[0] & (~0xFF) ;            // upper 24 bits of u24[0]
    u32[1]  = (u24[0] & 0xFF) << 24 ;       // upper 8 bits from lower 8 bits of u24[0]
    u32[1] |= ((u24[1] >> 8) & 0xFFFF00 ) ; // next 16 bits from upper 16 bits of u24[1]
    u32[2]  = u24[1] << 16 ;                // upper 16 bits from lower 16 bits of u24[1]
    u32[2] |= (u24[2] >> 16) & 0xFF00 ;     // next 8 bits from upper 8 bits of u24[2]
    u32[3]  = u24[2] << 8 ;                 // lower 24 bits of u24[2]
    u32 += 4 ;
    u24 += 3 ;
  }
#endif
  // process leftovers
// printf("fp32_bf24 bf24_fp32 = %d %8.8x %8.8x %8.8x\n", n3, u24[0], u24[1], u24[2]);
  if(n3 > 0) u32[0] = u24[0] & (~0xFF) ;    // upper 24 bits of u24[0]
  if(n3 > 1) {
    u32[1]  = (u24[0] & 0xFF) << 24 ;       // upper 8 bits from lower 8 bits of u24[0]
    u32[1] |= ((u24[1] >> 8) & 0xFFFF00 ) ; // next 16 bits from upper 16 bits of u24[1]
  }
  if(n3 > 2) {
    u32[2]  = u24[1] << 16 ;                // upper 16 bits from lower 16 bits of u24[1]
    u32[2] |= (u24[2] >> 16) & 0xFF00 ;     // next 8 bits from upper 8 bits of u24[2]
  }
}

// restore the lower 24 bits of 4 32 bit unsigned integers from 4 24 bit tokens in 3 32 bit words
void u24u32(uint32_t *u32, uint32_t *u24, int32_t n){
  int i ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32a, v24a, v32b, v24b, vm ;
  vx = _mm_loadu_si128((__m128i*) _32_to_lower24) ;
  for(i=0 ; i<n-7 ; i+=8){
    v24a = _mm_loadu_si128((__m128i*)  u24) ;
    v24b = _mm_loadu_si128((__m128i*) (u24+3)) ;
    v32a = _mm_shuffle_epi8(v24a, vx) ;
    v32b = _mm_shuffle_epi8(v24b, vx) ;
    _mm_storeu_si128((__m128i*)  u32,    v32a) ;
    _mm_storeu_si128((__m128i*) (u32+4), v32b) ;
    u32 += 8 ;
    u24 += 6 ;
  }
  if(i < n-3){
    v24a = _mm_loadu_si128((__m128i*) u24) ;
    v32a = _mm_shuffle_epi8(v24a, vx) ;
    _mm_storeu_si128((__m128i*) u32, v32a) ;
    u32 += 4 ;
    u24 += 3 ;
  }
#else
  for(i=0 ; i<n-3 ; i+=4){
    u32[0]  = u24[0] >> 8 ;             // upper 24 bits of u24[0]
    u32[1]  = (u24[0] & 0xFF) << 16 ;   // upper 8 bits from lower 8 bits of u24[0]
    u32[1] |= (u24[1] >> 16) ;          // next 16 bits from upper 16 bits of u24[1]
    u32[2]  = (u24[1] & 0xFFFF) << 8 ;  // upper 16 bits from lower 16 bits of u24[1]
    u32[2] |= (u24[2] >> 24) ;          // next 8 bits from upper 8 bits of u24[2]
    u32[3]  = u24[2] & 0xFFFFFF ;       // lower 24 bits of u24[2]
    u32 += 4 ;
    u24 += 3 ;
  }
#endif
  // process leftovers
  if(n3 > 0) u32[0] = u24[0] >> 8 ;     // upper 24 bits of u24[0]
  if(n3 > 1) {
    u32[1]  = (u24[0] & 0xFF) << 16 ;   // upper 8 bits from lower 8 bits of u24[0]
    u32[1] |= (u24[1] >> 16) ;          // next 16 bits from upper 16 bits of u24[1]
  }
  if(n3 > 2) {
    u32[2]  = (u24[1] & 0xFFFF) << 8 ;  // upper 16 bits from lower 16 bits of u24[1]
    u32[2] |= (u24[2] >> 24) ;          // next 8 bits from upper 8 bits of u24[2]
  }
}

// restore the lower 24 bits of 4 32 bit signed integers from 4 24 bit tokens in 3 32 bit words
// the most significant of the 24 bits is treated as a sign bit
// same as bf24_fp32 but with an 8 bit arithmetic right shift to propagate sign
void i24i32(int32_t *i32, uint32_t *u24, int32_t n){
  int i ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32a, v24a, v32b, v24b, vm ;
  vx = _mm_loadu_si128((__m128i*) _32_to_upper24) ;
  for(i=0 ; i<n-7 ; i+=8){
    v24a = _mm_loadu_si128((__m128i*)  u24) ;
    v24b = _mm_loadu_si128((__m128i*) (u24+3)) ;
    v32a = _mm_shuffle_epi8(v24a, vx) ;
    v32b = _mm_shuffle_epi8(v24b, vx) ;
    v32a = _mm_srai_epi32(v32a, 8) ;    // right shift by 8 positions
    v32b = _mm_srai_epi32(v32b, 8) ;    // right shift by 8 positions
    _mm_storeu_si128((__m128i*)  i32,    v32a) ;
    _mm_storeu_si128((__m128i*) (i32+4), v32b) ;
    i32 += 8 ;
    u24 += 6 ;
  }
  if(i < n-3) {
    v24a = _mm_loadu_si128((__m128i*) u24) ;
    v32a = _mm_shuffle_epi8(v24a, vx) ;
    v32a = _mm_srai_epi32(v32a, 8) ;    // right shift by 8 positions
    _mm_storeu_si128((__m128i*) i32, v32a) ;
    i32 += 4 ;
    u24 += 3 ;
  }
#else
  for(i=0 ; i<n-3 ; i+=4){
    i32[0]  = u24[0] ;                      // upper 24 bits of u24[0]
    i32[0] >>= 8 ;
    i32[1]  = (u24[0] & 0xFF) << 24 ;       // upper 8 bits from lower 8 bits of u24[0]
    i32[1] |= ((u24[1] >> 8) & 0xFFFF00 ) ; // next 16 bits from upper 16 bits of u24[1]
    i32[1] >>= 8 ;
    i32[2]  = u24[1] << 16 ;                // upper 16 bits from lower 16 bits of u24[1]
    i32[2] |= (u24[2] >> 16) & 0xFF00 ;     // next 8 bits from upper 8 bits of u24[2]
    i32[2] >>= 8 ;
    i32[3]  = (u24[2] << 8) >> 8 ;          // lower 24 bits of u24[2]
    i32 += 4 ;
    u24 += 3 ;
  }
#endif
  // process leftovers
  if(n3 > 0) {
    i32[0] = u24[0] ;                       // upper 24 bits of u24[0]
    i32[0] >>= 8 ;
  }
  if(n3 > 1) {
    i32[1]  = (u24[0] & 0xFF) << 24 ;       // upper 8 bits from lower 8 bits of u24[0]
    i32[1] |= ((u24[1] >> 8) & 0xFFFF00 ) ; // next 16 bits from upper 16 bits of u24[1]
    i32[1] >>= 8 ;
  }
  if(n3 > 2) {
    i32[2]  = u24[1] << 16 ;                // upper 16 bits from lower 16 bits of u24[1]
    i32[2] |= (u24[2] >> 16) & 0xFF00 ;     // next 8 bits from upper 8 bits of u24[2]
    i32[2] >>= 8 ;
  }
}
#if defined(SELF_TEST)
// int verify_upper(uint32_t *a, uint32_t *b, int n);
// int verify_lower(uint32_t *a, uint32_t *b, int n);
int verify_masked(void *a, void *b, uint32_t mask, int n);
void apply_round(void *a, void *b, uint32_t round, int n);

#define NPTS 1031
#define NTIMES 10000

#include <misc_timers.h>

int main(int argc, char **argv){
  int i, j, nt, errors,  n = 7 ;
  static uint32_t fp32[8] = { 0x00810203, 0x04850607, 0x08890A0B, 0x0C8D0E0F, \
                              0x10911213, 0x14951617, 0x18991A1B, 0x1C9D1E1F  } ;
  static uint32_t uf24[8]  ;
  static uint32_t fp32b [8] ;
  uint32_t lfp32[NPTS+1], luf24[NPTS+1], lfp32b[NPTS+1] ;
  float float32[NPTS], float32b[NPTS], float32m[NPTS] ;
  uint32_t *ifloat32 = (uint32_t *) float32 ;
  uint32_t *ifloat32b = (uint32_t *) float32b ;
  uint32_t *ifloat32m = (uint32_t *) float32m ;
  uint64_t t0, t1, tmin, tmax, tavg, freq ;
  double nano ;
  uint64_t t[NTIMES] ;
  double bias, avg ;

  freq = cycles_counter_freq() ;
  nano = 1000000000 ;
  nano /= freq ;
  for(i=0 ; i<NPTS ; i++) { lfp32[i] = i ; luf24[i] = 0 ; lfp32b[i] = 0 ; }
//   for(i=0 ; i<NPTS ; i++) { float32[i] = NPTS + (i - .111)/NPTS ; float32b[i] = -1 ; }
  for(i=0 ; i<NPTS ; i++) { float32[i] = (i + .111)/NPTS ; float32b[i] = -1 ; }

  memset(uf24, 0xFF, sizeof(uf24));
  fp32_bf24((void *) fp32, uf24, n) ;
  printf("fp32      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;

  printf("uf24      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;

  memset(fp32b, 0xFF, sizeof(fp32b));
  bf24_fp32(fp32b, uf24, n) ;
//   errors = verify_upper((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  errors = verify_masked((uint32_t *) fp32, (uint32_t *) fp32b, 0xFFFFFF00u, n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;

  printf(">bf24     : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    fp32_bf24((void *) float32, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf(">fp32     : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    bf24_fp32(float32b, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_upper((uint32_t *) float32, (uint32_t *) float32b, NPTS) ;
  apply_round(float32, float32m, 0x80, NPTS);
  errors = verify_masked((uint32_t *) float32m, (uint32_t *) float32b, 0xFFFFFF00u,  NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;
  bias = 0.0 ; avg = 0.0 ;
  for(i=0 ; i<NPTS ; i++) {bias += (float32[i] - float32b[i]) ; avg += float32[i] ; }
  bias /= NPTS ; avg /= NPTS ;
  printf("sample float : %f %f, relative error = %8.3g, avg = %f, bias = %8.3g\n", 
         float32[NPTS/2], float32b[NPTS/2], (float32[NPTS/2]-float32b[NPTS/2])/float32[NPTS/2], avg, bias) ;
  printf("===============================================================\n") ;

  memset(uf24, 0xFF, sizeof(uf24));
  fp32_bf16((void *) fp32, uf24, n) ;
  printf("fp32      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;

  printf("uf16      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;

  memset(fp32b, 0xFF, sizeof(fp32b));
  bf16_fp32(fp32b, uf24, n) ;
//   errors = verify_upper((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  errors = verify_masked((uint32_t *) fp32, (uint32_t *) fp32b, 0xFFFF0000u, n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;

  printf(">bf16     : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    fp32_bf16((void *) float32, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  apply_round(float32, float32m, 0x8000, NPTS);
  for(i=0 ; i<16 ; i++) printf("%8.8x ",ifloat32m[i]) ; printf("\n") ;
  printf("luf16     : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x          ",luf24[i])     ; printf("\n") ;
//   fp32_nf16((void *) float32, luf24, NPTS, 8) ;
//   printf("luf16b    : ");
//   for(i=0 ; i<8 ; i++) printf("%8.8x          ",luf24[i])     ; printf("\n") ;
  printf("fp32      : ");
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf(">fp32     : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    bf16_fp32(float32b, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_upper((uint32_t *) float32, (uint32_t *) float32b, NPTS) ;
//   for(i=0 ; i<16 ; i++) printf("%8.8x ",ifloat32b[i]) ; printf("\n") ;
  apply_round(float32, float32m, 0x8000, NPTS);
  errors = verify_masked((uint32_t *) float32m, (uint32_t *) float32b, 0xFFFF0000u,  NPTS) ;
//   printf(">fp32     : ");
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;
  bias = 0.0 ; avg = 0.0 ;
  for(i=0 ; i<NPTS ; i++) {bias += (float32[i] - float32b[i]) ; avg += float32[i] ; }
  bias /= NPTS ; avg /= NPTS ;
  printf("sample float : %f %f, relative error = %8.3e, avg = %f, bias = %8.3g\n", 
         float32[NPTS/2], float32b[NPTS/2], (float32[NPTS/2]-float32b[NPTS/2])/float32[NPTS/2], avg, bias) ;
//   printf("===============================================================\n") ;

  printf(">nf16-8   : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    fp32_nf16((void *) float32, luf24, NPTS, 8) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  apply_round(float32, float32m, 0x8000, NPTS);
  for(i=0 ; i<16 ; i++) printf("%8.8x ",ifloat32m[i]) ; printf("\n") ;
  printf("luf16     : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x          ",luf24[i])     ; printf("\n") ;
  fp32_nf16((void *) float32, luf24, NPTS, 8) ;
//   printf("luf16b    : ");
//   for(i=0 ; i<8 ; i++) printf("%8.8x          ",luf24[i])     ; printf("\n") ;
  printf("fp32      : ");
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf(">fp32     : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    nf16_fp32(float32b, luf24, NPTS, 8) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_upper((uint32_t *) float32, (uint32_t *) float32b, NPTS) ;
  apply_round(float32, float32m, 0x8000, NPTS);
  errors = verify_masked((uint32_t *) float32m, (uint32_t *) float32b, 0xFFFF0000u,  NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;
  bias = 0.0 ; avg = 0.0 ;
  for(i=0 ; i<NPTS ; i++) {bias += (float32[i] - float32b[i]) ; avg += float32[i] ; }
  bias /= NPTS ; avg /= NPTS ;
  printf("sample float : %f %f, relative error = %8.3e, avg = %f, bias = %8.3g\n", 
         float32[NPTS/2], float32b[NPTS/2], (float32[NPTS/2]-float32b[NPTS/2])/float32[NPTS/2], avg, bias) ;
  printf("===============================================================\n") ;
return 0 ;
  memset(uf24, 0xFF, sizeof(uf24));
  u32_u24(fp32, uf24, n) ;
  printf("u32       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    u32_u24(lfp32, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("u24       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;
  memset(fp32b, 0xFF, sizeof(fp32b));
  u24u32(fp32b, uf24, n) ;
//   errors = verify_lower((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  errors = verify_masked((uint32_t *) fp32, (uint32_t *) fp32b, 0xFFFFFF, n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    u24u32(lfp32b, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_lower((uint32_t *) lfp32, (uint32_t *) lfp32b, NPTS) ;
  errors = verify_masked((uint32_t *) lfp32, (uint32_t *) lfp32b, 0xFFFFFF, NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;

  printf("===============================================================\n") ;

  memset(uf24, 0xFF, sizeof(uf24));
  u32_u24(fp32, uf24, n) ;
  printf("i32       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    u32_u24(lfp32, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("i24       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;
  memset(fp32b, 0xFF, sizeof(fp32b));
  i24i32((int32_t *) fp32b, uf24, n) ;
//   errors = verify_lower((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  errors = verify_masked((uint32_t *) fp32, (uint32_t *) fp32b, 0xFFFFFF,  n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    i24i32((int32_t *) lfp32b, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 3*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_lower((uint32_t *) lfp32, (uint32_t *) lfp32b, NPTS) ;
  errors = verify_masked((uint32_t *) lfp32, (uint32_t *) lfp32b, 0xFFFFFF,  NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;

  printf("\n\n\n\n\n");
}
#endif
