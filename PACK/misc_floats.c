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
#include <rmn/misc_types.h>
#include <misc_floats.h>

#if ! defined(STATIC)
#define STATIC static
#endif

#include <with_simd.h>
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)

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

STATIC inline uint32_t f32_max_exp(void *f32, int32_t n){
  int i, n15 ;
  uint32_t sign, maxabs, tmp ;
  uint32_t *u32 = f32 ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vsign, vfabs, vtemp, vmask ;
  __m128i vsign0, vsign1, vfmax0, vfmax1 ;
#endif
  sign = 0 ;
  maxabs = 0 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  n15 = (n & 15) ;
  if(n15 == 0) n15 = 16 ;
  vmask = _mm256_set1_epi32(0x7FFFFFFFu) ;        // mask out leading bit

  vsign = _mm256_loadu_si256((__m256i *)  u32) ;
  vfabs = _mm256_and_si256(vmask, vsign) ;
  if(n15 > 7){
    vtemp = _mm256_loadu_si256((__m256i *) (u32+8)) ;
    vsign = _mm256_or_si256(vsign, vtemp) ;
    vtemp = _mm256_and_si256(vmask, vtemp) ;
    vfabs = _mm256_max_epu32(vfabs, vtemp) ;
  }
  u32 += n15 ;
  for(i = n15; i < n-15 ; i += 16){
    vtemp = _mm256_loadu_si256((__m256i *)  u32   ) ;    // fetch 8 floats
    vsign = _mm256_or_si256(vsign, vtemp) ;              // OR with sign
    vtemp = _mm256_and_si256(vmask, vtemp) ;             // absolute value
    vfabs = _mm256_max_epu32(vfabs, vtemp) ;             // MAX absolute value

    vtemp = _mm256_loadu_si256((__m256i *) (u32+8)) ;    // fetch next 8 floats
    vsign = _mm256_or_si256(vsign, vtemp) ;              // OR with sign
    vtemp = _mm256_and_si256(vmask, vtemp) ;             // absolute value
    vfabs = _mm256_max_epu32(vfabs, vtemp) ;             // MAX absolute value
    u32 += 16 ;
  }
  vsign0 = _mm256_extracti128_si256(vsign, 0) ;
  vsign1 = _mm256_extracti128_si256(vsign, 1) ;
  vfmax0 = _mm256_extracti128_si256(vfabs, 0) ;
  vfmax1 = _mm256_extracti128_si256(vfabs, 1) ;
  vsign0 = _mm_or_si128(vsign0, vsign1) ;
  vfmax0 = _mm_max_epu32(vfmax0, vfmax1) ;
  vsign1 = _mm_shuffle_epi32(vsign0, 0xF) ;       // bring words 2 and 3 into positions 0 and 1
  vfmax1 = _mm_shuffle_epi32(vfmax0, 0xF) ;       // bring words 2 and 3 into positions 0 and 1
  vsign0 = _mm_or_si128(vsign0, vsign1) ;
  vfmax0 = _mm_max_epu32(vfmax0, vfmax1) ;
  vsign1 = _mm_shuffle_epi32(vsign0, 1) ;       // bring words 1 into position 0
  vfmax1 = _mm_shuffle_epi32(vfmax0, 1) ;       // bring words 1 into position 0
  vsign0 = _mm_or_si128(vsign0, vsign1) ;
  vfmax0 = _mm_max_epu32(vfmax0, vfmax1) ;
  _mm_storeu_si32( (&sign),   vsign0) ;
  _mm_storeu_si32( (&maxabs), vfmax0) ;
#else
  for(i = 0; i < n ; i ++){
    tmp   = u32[i] ;
    sign |= tmp ;
    tmp  &= 0x7FFFFFFFu ;
    maxabs = (tmp > maxabs) ? tmp : maxabs ;
  }
#endif
  sign &= 0x80000000u ;
  return sign | maxabs ;
}

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
//   int tmp = f32_max_exp(f32, n) ;
//   if(tmp == 0) i = 1 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  for(    ; i < n-15 ; i += 16){
    v32  = _mm256_loadu_si256((__m256i*)  u32) ;     // get 8 floats
    v32  = _mm256_add_epi32(v32, vr) ;               // add rounding term
    v32  = _mm256_srli_epi32(v32,16) ;               // to lower 16 bits
    v16a = _mm256_extracti128_si256(v32, 0) ;        // get lower 128 bits
    v16b = _mm256_extracti128_si256(v32, 1) ;        // get upper 128 bits
    v16a = _mm_shuffle_epi8(v16a, vx) ;              // shuffle to lower 64 bits
    v16b = _mm_shuffle_epi8(v16b, vx) ;              // shuffle to lower 64 bits
    v16c = _mm_unpacklo_epi64(v16a, v16b) ;          // merge to 128 bits
    _mm_storeu_si128((__m128i *) u16  , v16c) ;      // store 4 words
    v32  = _mm256_loadu_si256((__m256i*)  (u32+8)) ; // get next 8 floats
    v32  = _mm256_add_epi32(v32, vr) ;
    v32  = _mm256_srli_epi32(v32,16) ;
    v16a = _mm256_extracti128_si256(v32, 0) ;
    v16b = _mm256_extracti128_si256(v32, 1) ;
    v16a = _mm_shuffle_epi8(v16a, vx) ;
    v16b = _mm_shuffle_epi8(v16b, vx) ;
    v16c = _mm_unpacklo_epi64(v16a, v16b) ;
    _mm_storeu_si128((__m128i *) (u16+4), v16c) ;      // store next 4 words
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

// return maximum absolute value that can be represented by a 16 bit float
// using nexp (<=8) bits for the exponent field
float nf16_floatmax(int nexp){
  FloatUint top ;
  uint32_t s8 = (8 - nexp) ;
  uint32_t fexp = (0x7F >> s8) ;
  if(nexp == 8) fexp-- ;
  top.u = ((128+fexp) << 23) | 0x7FFFFF ;
  top.u >>= (9 + nexp) ;
  top.u <<= (9 + nexp) ;
  return top.f ;
}

// pack upper 16 bits from 4 32 bit words into 2 32 bit unsigned integers
// normally used to pack floats 
// 1 bit : sign / nexp bits : exponent / (16 - nexp - 1) bits : mantissa
// when nexp == 8, "brain float 16" is produced (same as fp32_bf16)
// NOTE: add overflow to infinity argument, set top = 0x7FFF0000u in that case
void fp32_nf16(void *f32, uint32_t *u16, int32_t n, int32_t nexp){
  int32_t i ;
  FloatUint *u32 = f32 ;
  uint32_t round, fexp, sign, top ;
  uint32_t s8 ;
  FloatUint nfac, z ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vi, v8, vs, vt, vr, vn, v1, v32, vm ;
  __m256  vf ;
  __m128i v16a, v16b, v16c, vx ;
#endif
  i=0 ;
//   int tmp = f32_max_exp(f32, n) ;
//   if(tmp == 0) i = 1 ;
  s8 = (8 - nexp) ;                           // shift count for exp / mantissa
  round = (0x8000 >> s8) ;                    // rounding term to be added before clipping mantissa
  fexp = (0x7F >> s8) ;                       // new offset for exponent
  nfac.u = (fexp  << 23) ;                    // normalisation factor
  //  largest source absolute value that can be represented = ((128+fexp) << 23) | 0x7FFFFF
  //   top = ((fexp+fexp+1) << 23) | 0x7FFFFF ;    // maximum allowable value after normalisation
  top = 0x7FFE0000u ;                         // maximum value once coded
  if(nexp == 8) top = 0x7FFF0000u ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  vr = _mm256_set1_epi32(round) ;
  vs = _mm256_set1_epi32(s8) ;
//   vt = _mm256_set1_epi32(top) ;
//   vt = _mm256_set1_epi32(0x7FFF0000u) ;
  vt = _mm256_set1_epi32(top) ;
  vn = _mm256_set1_epi32(nfac.u) ;
  vm = _mm256_set1_epi32(0x80000000u) ;    // sign mask
  vx = _mm_loadu_si128((__m128i*) _32_from_lower16) ;
  for(    ; i < n-15 ; i += 16){
    vf = _mm256_loadu_ps((float *) u32 ) ;          // get 8 floats
    vi = (__m256i) _mm256_mul_ps(vf, (__m256) vn) ; // normalize
    v1 = _mm256_and_si256(vm, vi) ;                 // get sign
    vi = _mm256_andnot_si256(vm, vi) ;              // absolute value
    vi = _mm256_add_epi32(vr, vi) ;                 // add rounding factor
//     vi = _mm256_min_epi32(vt, vi) ;                 // "infinity" clip
    vi = _mm256_sllv_epi32(vi, vs) ;                // shift
    vi = _mm256_min_epi32(vt, vi) ;                 // "infinity" clip
    vi = _mm256_or_si256(v1, vi) ;                  // restore sign
    vi  = _mm256_srli_epi32(vi,16) ;                // shift right by 16 bits
    v16a = _mm256_extracti128_si256(vi, 0) ;        // lower 128 bits
    v16b = _mm256_extracti128_si256(vi, 1) ;        // upper 128 bits
    v16a = _mm_shuffle_epi8(v16a, vx) ;             // shuffle to lower 64 bits
    v16b = _mm_shuffle_epi8(v16b, vx) ;             // shuffle to lower 64 bits
    v16c = _mm_unpacklo_epi64(v16a, v16b) ;         // merge to 128 bits
    _mm_storeu_si128((__m128i *) u16, v16c) ;       // store 4 words

    vf = _mm256_loadu_ps((float *) (u32+8) ) ;      // get next 8 floats
    vi = (__m256i) _mm256_mul_ps(vf, (__m256) vn) ; // normalize
    v1 = _mm256_and_si256(vm, vi) ;                 // get sign
    vi = _mm256_andnot_si256(vm, vi) ;              // absolute value
    vi = _mm256_add_epi32(vr, vi) ;                 // add rounding factor
//     vi = _mm256_min_epi32(vt, vi) ;                 // "infinity" clip
    vi = _mm256_sllv_epi32(vi, vs) ;                // shift
    vi = _mm256_min_epi32(vt, vi) ;                 // "infinity" clip
    vi = _mm256_or_si256(v1, vi) ;                  // restore sign
    vi  = _mm256_srli_epi32(vi,16) ;                // shift right by 16 bits
    v16a = _mm256_extracti128_si256(vi, 0) ;
    v16b = _mm256_extracti128_si256(vi, 1) ;
    v16a = _mm_shuffle_epi8(v16a, vx) ;
    v16b = _mm_shuffle_epi8(v16b, vx) ;
    v16c = _mm_unpacklo_epi64(v16a, v16b) ;
    _mm_storeu_si128((__m128i *) (u16+4), v16c) ;   // store next 4 words
    u32 += 16 ;
    u16 +=  8 ;
  }
#endif
  // process leftovers
// printf("fp32_nf16, fac = %g, round = %8.8x, s8 = %d", nfac.f, round, s8) ;
// int clip = 0;
  for(    ; i < n-1 ; i += 2){
    z.f = u32[0].f ;
    sign = z.u & 0x80000000u ;           // keep sign
    z.u &= 0x7FFFFFFFu ;                 // absolute value
    z.f *= nfac.f ;                      // normalize
    z.u += round ;                       // add rounding
//     z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
    z.u <<= s8 ;                         // shift
    z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
    z.u |= sign ;                        // restore sign
    u16[0]  =  z.u & top ;               // upper 16 bits to upper 16 bits of u16

    z.f = u32[1].f ;
    sign = z.u & 0x80000000u ;           // keep sign
    z.u &= 0x7FFFFFFFu ;                 // absolute value
    z.f *= nfac.f ;                      // normalize
    z.u += round ;                       // add rounding
//     z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
    z.u <<= s8 ;                         // shift
    z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
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
//     z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
    z.u <<= s8 ;                         // shift
    z.u = (z.u > top) ? top : z.u ;      // "infinity" clip
    z.u |= sign ;                        // restore sign
    u16[0]  =  z.u & 0xFFFF0000u ;       // upper 16 bits to upper 16 bits of u16
  }
// printf("clip = %d\n            ", clip) ;
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
  __m256i v32, va, vb, vl, vh, vm ;
  __m128i v16 ;
#endif
  i=0 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  vm = _mm256_set1_epi32(0xFFFF0000) ;
  for(    ; i < n-15 ; i += 16){
    v32 = _mm256_loadu_si256((__m256i *) u16) ;      // 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b c, d,e ,f
    vl  = _mm256_and_si256(v32, vm) ;                // 0, 2, 4, 6, 8, a, c, e
    vh  = _mm256_slli_epi32(v32, 16) ;               // 1, 3, 5, 7, 9, b, d, f
    va  = _mm256_unpacklo_epi32(vl, vh) ;            // 0, 1, 2, 3, 8, 9, a, b
    vb  = _mm256_unpackhi_epi32(vl, vh) ;            // 4, 5, 6, 7, c, d, e, f
    vl  = _mm256_permute2x128_si256(va, vb, 32) ;    // 0, 1, 2, 3, 4, 5, 6, 7
    vh  = _mm256_permute2x128_si256(va, vb, 49) ;    // 8, 9, a, b, c, d, e, f
    _mm256_storeu_si256((__m256i *) (u32),   vl) ;
    _mm256_storeu_si256((__m256i *) (u32+8), vh) ;
//     v16 = _mm_loadu_si128((__m128i *) u16) ;
//     v32 = _mm256_cvtepu16_epi32(v16) ;
//     v32 = _mm256_slli_epi32(v32, 16) ;
//     v32 = _mm256_shuffle_epi32(v32, 0xB1) ;
//     _mm256_storeu_si256((__m256i *) u32, v32) ;
//     v16 = _mm_loadu_si128((__m128i *) (u16+4)) ;
//     v32 = _mm256_cvtepu16_epi32(v16) ;
//     v32 = _mm256_slli_epi32(v32, 16) ;
//     v32 = _mm256_shuffle_epi32(v32, 0xB1) ;
//     _mm256_storeu_si256((__m256i *) (u32+8), v32) ;
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
// NOTE: may want to add infinity option, and decode 7FFF as infinity
void nf16_fp32(void *f32, uint32_t *u16, int32_t n, int nexp){
  int32_t i ;
  FloatUint *u32 = f32 ;
  uint32_t top, fexp, sign ;
  uint32_t s8 ;
  FloatUint nfac, z ;
//   bf16_fp32(f32, u16, n) ;
//   return ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vi, v8, vs, vt, vn, v1, v32, vm, vu, v1l, v1h, vl, vh, va,vb ;
  __m256  vf, vfl, vfh ;
  __m128i v16 ;
#endif
  i = 0 ;
  s8 = (8 - nexp) ;                           // shift count for exp / mantissa
  fexp = (0x7F >> s8) ;                       // new offset for exponent
  fexp = 254 - fexp ;                         // new exponent for normalisation factor
  nfac.u = (fexp  << 23) ;                    // normalisation factor
//   top = nfac.u | 0x7FFFFF ;                   // maximum allowable value after normalisation
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  vs = _mm256_set1_epi32(s8) ;
  vu = _mm256_set1_epi32(0xFFFF0000) ;
  vt = _mm256_set1_epi32(top) ;
  vn = _mm256_set1_epi32(nfac.u) ;
  vm = _mm256_set1_epi32(0x80000000u) ;    // sign mask
  for(    ; i < n-15 ; i += 16){
    v32 = _mm256_loadu_si256((__m256i *) u16) ;      // 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b c, d,e ,f
    vl  = _mm256_and_si256(v32, vu) ;                // 0, 2, 4, 6, 8, a, c, e  (lower portion)
    vh  = _mm256_slli_epi32(v32, 16) ;               // 1, 3, 5, 7, 9, b, d, f  (upper portion)
    va  = _mm256_unpacklo_epi32(vl, vh) ;            // 0, 1, 2, 3, 8, 9, a, b
    vb  = _mm256_unpackhi_epi32(vl, vh) ;            // 4, 5, 6, 7, c, d, e, f
    vl  = _mm256_permute2x128_si256(va, vb, 32) ;    // 0, 1, 2, 3, 4, 5, 6, 7
    vh  = _mm256_permute2x128_si256(va, vb, 49) ;    // 8, 9, a, b, c, d, e, f
    v1l = _mm256_and_si256(vm, vl) ;                 // get sign of vl
    v1h = _mm256_and_si256(vm, vh) ;                 // get sign of vh
    vl  = _mm256_andnot_si256(vm, vl) ;              // absolute value of vl
    vh  = _mm256_andnot_si256(vm, vh) ;              // absolute value of vh
    vl  = _mm256_srlv_epi32(vl, vs) ;                // shift vl to proper position
    vh  = _mm256_srlv_epi32(vh, vs) ;                // shift vh to proper position
    vl  = _mm256_or_si256(v1l, vl) ;                 // restore sign of vl
    vh  = _mm256_or_si256(v1h, vh) ;                 // restore sign of vh
    vfl = _mm256_mul_ps((__m256) vl, (__m256) vn) ;  // apply normalization factor to vl
    vfh = _mm256_mul_ps((__m256) vh, (__m256) vn) ;  // apply normalization factor to vh
    _mm256_storeu_ps((float *) (u32  ), vfl) ;       // store 8 floats from vfl
    _mm256_storeu_ps((float *) (u32+8), vfh) ;       // store 8 floats from vfl

    u32 +=  16 ;
    u16 +=   8 ;
  }
//   for(    ; i < n-7 ; i += 8){
//     v16 = _mm_loadu_si128((__m128i *) u16) ;         // get 4 words of packed floats
//     v32 = _mm256_cvtepu16_epi32(v16) ;               // expand to 8 words
//     v32 = _mm256_slli_epi32(v32, 16) ;               // shift left 16 positions
//     v1  = _mm256_and_si256(vm, v32) ;                // get sign
//     v32 = _mm256_andnot_si256(vm, v32) ;             // absolute value
//     v32 = _mm256_shuffle_epi32(v32, 0xB1) ;          // shuffle to proper position (undo cvtepu16 LE swap)
//     v32 = _mm256_srlv_epi32(v32, vs) ;               // shift to proper position
//     v32 = _mm256_or_si256(v1, v32) ;                 // restore sign
//     vf  = _mm256_mul_ps((__m256) v32, (__m256) vn) ; // apply normalization factor
//     _mm256_storeu_ps((float *) u32, vf) ;            // store 8 floats
// 
//     u32 +=  8 ;
//     u16 +=  4 ;
//   }
#endif
  // process leftovers
// printf("nf16_fp32, fac = %g", nfac.f) ;
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
// if(i==0) printf(", u16 = %8.8x, f32 = %g %g", u16[0], u32[0].f, u32[1].f);
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
// printf("\n            ") ;
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
void u24i32(int32_t *i32, uint32_t *u24, int32_t n){
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
    v32a = _mm_srai_epi32(v32a, 8) ;    // arithmetic right shift by 8 positions
    v32b = _mm_srai_epi32(v32b, 8) ;    // arithmetic right shift by 8 positions
    _mm_storeu_si128((__m128i*)  i32,    v32a) ;
    _mm_storeu_si128((__m128i*) (i32+4), v32b) ;
    i32 += 8 ;
    u24 += 6 ;
  }
  if(i < n-3) {
    v24a = _mm_loadu_si128((__m128i*) u24) ;
    v32a = _mm_shuffle_epi8(v24a, vx) ;
    v32a = _mm_srai_epi32(v32a, 8) ;    // arithmetic right shift by 8 positions
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
    i32[3]  = (u24[2] << 8)      ;          // lower 24 bits of u24[2]
    i32[3] >>= 8 ;
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
