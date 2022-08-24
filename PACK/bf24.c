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

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
#include <emmintrin.h>
#include <immintrin.h>

// shuffle table to pack the upper 24 bits from 4 32 bit words into 3 32 bit words
static uint8_t upper_32_to_24[] = { 7, 1, 2, 3, 10, 11, 5, 6, 13, 14, 15, 9, 128, 128, 128, 128 } ;
// shuffle table to pack the lower 24 bits from 4 32 bit words into 3 32 bit words
static uint8_t lower_32_to_24[] = { 6, 0, 1, 2,9, 10, 4, 5, 12, 13, 14, 8, 128, 128, 128, 128 } ;

// shuffle table to extract 4 24 bit items from 3 32 bit words into the upper part of a 32 bit word
static uint8_t upper_24_to_32[] = {128, 1, 2, 3, 128, 6, 7, 0, 128, 11, 4, 5, 128, 8, 9, 10 } ;
// shuffle table to extract 4 24 bit items from 3 32 bit words into the lower part of a 32 bit word
static uint8_t lower_24_to_32[] = { 1, 2, 3, 128, 6, 7, 0, 128, 11, 4, 5, 128, 8, 9, 10, 128 } ;

// mask for masked store of first 3 elements in vector
static uint32_t mask24[4] = { 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0 } ;

#endif
/*
void c_fp32_bf24(void *f32, void *f24, uint32_t n){
  int i, j ;
  uint8_t *c32 = (uint8_t *) f32 ;
  uint8_t *c24 = (uint8_t *) f24 ;
  for(i=0 ; i<n-3 ; i++){
    for(j=0 ; j<12 ; j++){
      c24[j] = c32[f32_to_24[j]] ;
    }
    c32 += 16 ;
    c24 += 12 ;
  }
}*/
/*
void i_fp32_bf24(void *f32, void *f24, uint32_t n){
  int i ;
  uint32_t *fp32 = f32 ;
  uint32_t *bf24 = f24 ;
  uint32_t mask = ~0xFF ;
  for(i=0 ; i<n-3 ; i++){
    bf24[0] =  (fp32[0] & mask)        | (fp32[1] >> 24) ;
    bf24[1] = ((fp32[1] & mask) <<  8) | (fp32[2] >> 16) ;
    bf24[2] = ((fp32[2] & mask) << 16) | (fp32[3] >> 8) ;
    fp32 += 4 ;
    bf24 += 3 ;
  }
}*/

// pack upper 24 bits from 4 32 bit words into 3 32 bit unsigned integers
// normally used to pack floats
void fp32_bf24(void *f32, uint32_t *u24, uint32_t n){
  int i ;
  uint32_t *u32 = f32 ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32, v24, vm ;
  vx = _mm_loadu_si128((__m128i*) upper_32_to_24) ;
  vm = _mm_loadu_si128((__m128i*) mask24) ;
  for(i=0 ; i<n-7 ; i+=8){   // 8 words into 6 words
    v32 = _mm_loadu_si128((__m128i*) u32) ;
    v24 = _mm_shuffle_epi8(v32, vx) ;
    _mm_storeu_si128((__m128i*) u24, v24) ;      // plenty of extra room
    u32 += 4 ;
    u24 += 3 ;
    v32 = _mm_loadu_si128((__m128i*) u32) ;
    v24 = _mm_shuffle_epi8(v32, vx) ;
    _mm_storeu_si128((__m128i*) u24, v24) ;      // plenty of extra room
    u32 += 4 ;
    u24 += 3 ;
  }
  if(i < n-3){   // 4 words into 3 words
    v32 = _mm_loadu_si128((__m128i*) u32) ;
    v24 = _mm_shuffle_epi8(v32, vx) ;
    _mm_maskstore_epi32((int *) u24, vm, v24) ;  // mask used to store only 3 words
    u32 += 4 ;
    u24 += 3 ;
  }
#else
  for(i=0 ; i<n-3 ; i+=4){
    u24[0] = u32[0] & (~0xFF) ;             // upper 24 bits of u32[0]
    u24[0] |= u32[1] >> 24 ;                // upper 8 bits of u32[1]
    u24[1]  = (u32[1] & (~0xFF)) << 8 ;     // lower 16 of upper 24 bits of u32[1]
    u24[1] |= (u32[2] >> 16) ;              // upper 16 of upper 24 bits of u32[2]
    u24[2]  = (u32[2] & (~0xFF)) << 16 ;    // lower 8 of upper 24 bits of u32[2]
    u24[2] |= (u32[3] >> 8) ;               // upper 24 bits of u32[3]
    u32 += 4 ;
    u24 += 3 ;
  }
#endif
  // process leftovers
  if(n3 > 0) u24[0] = u32[0] & (~0xFF) ;    // upper 24 bits of u32[0]
  if(n3 > 1) {
    u24[0] |= u32[1] >> 24 ;                // upper 8 bits of u32[1]
    u24[1]  = (u32[1] & (~0xFF)) << 8 ;     // lower 16 of upper 24 bits of u32[1]
  }
  if(n3 > 2) {
    u24[1] |= (u32[2] >> 16) ;              // upper 16 of upper 24 bits of u32[2]
    u24[2]  = (u32[2] & (~0xFF)) << 16 ;    // lower 8 of upper 24 bits of u32[2]
  }
}

// pack lower 24 bits from 4 32 bit unsigned integers into 3 32 bit unsigned integers
void u32_u24(uint32_t *u32, uint32_t *u24, uint32_t n){
  int i ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32, v24, vm ;
  vx = _mm_loadu_si128((__m128i*) lower_32_to_24) ;
  vm = _mm_loadu_si128((__m128i*) mask24) ;
//   for(i=0 ; i<n-3 ; i+=4){
//     v32 = _mm_loadu_si128((__m128i*) u32) ;
//     v24 = _mm_shuffle_epi8(v32, vx) ;
//     if(i < n-4) {
//       _mm_storeu_si128((__m128i*) u24, v24) ;      // plenty of extra room
//     }else{
//       _mm_maskstore_epi32((int *) u24, vm, v24) ;  // mask used to store only 3 words
//     }
//     u32 += 4 ;
//     u24 += 3 ;
//   }
  for(i=0 ; i<n-7 ; i+=8){
    v32 = _mm_loadu_si128((__m128i*) u32) ;
    v24 = _mm_shuffle_epi8(v32, vx) ;
    _mm_storeu_si128((__m128i*) u24, v24) ;      // plenty of extra room
    u32 += 4 ;
    u24 += 3 ;
    v32 = _mm_loadu_si128((__m128i*) u32) ;
    v24 = _mm_shuffle_epi8(v32, vx) ;
    _mm_storeu_si128((__m128i*) u24, v24) ;      // plenty of extra room
    u32 += 4 ;
    u24 += 3 ;
  }
  if(i < n-3){   // 4 words into 3 words
    v32 = _mm_loadu_si128((__m128i*) u32) ;
    v24 = _mm_shuffle_epi8(v32, vx) ;
    _mm_maskstore_epi32((int *) u24, vm, v24) ;  // mask used to store only 3 words
    u32 += 4 ;
    u24 += 3 ;
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
void i32_u24(int32_t *i32, uint32_t *u24, uint32_t n){
  u32_u24((uint32_t *) i32, u24, n) ;
}

// restore the upper 24 bits of 4 32 bit words from 4 24 bit tokens in 3 32 bit words
// normally used to restore floats
void bf24_fp32(void *f32, uint32_t *u24, uint32_t n){
  int i ;
  uint32_t *u32 = f32 ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32, v24 ;
  vx = _mm_loadu_si128((__m128i*) upper_24_to_32) ;
  for(i=0 ; i<n-7 ; i+=8){
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    _mm_storeu_si128((__m128i*) u32, v32) ;
    u32 += 4 ;
    u24 += 3 ;
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    _mm_storeu_si128((__m128i*) u32, v32) ;
    u32 += 4 ;
    u24 += 3 ;
  }
  if(i < n-3){
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    _mm_storeu_si128((__m128i*) u32, v32) ;
    u32 += 4 ;
    u24 += 3 ;
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
void u24u32(uint32_t *u32, uint32_t *u24, uint32_t n){
  int i ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32, v24 ;
  vx = _mm_loadu_si128((__m128i*) lower_24_to_32) ;
  for(i=0 ; i<n-7 ; i+=8){
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    _mm_storeu_si128((__m128i*) u32, v32) ;
    u32 += 4 ;
    u24 += 3 ;
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    _mm_storeu_si128((__m128i*) u32, v32) ;
    u32 += 4 ;
    u24 += 3 ;
  }
  if(i < n-3){
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    _mm_storeu_si128((__m128i*) u32, v32) ;
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
void i24i32(int32_t *i32, uint32_t *u24, uint32_t n){
  int i ;
  uint32_t n3 = (n & 3) ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vx, v32, v24 ;
  vx = _mm_loadu_si128((__m128i*) upper_24_to_32) ;
  for(i=0 ; i<n-7 ; i+=8){
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    v32 = _mm_srai_epi32(v32, 8) ;    // right shift by 8 positions
    _mm_storeu_si128((__m128i*) i32, v32) ;
    i32 += 4 ;
    u24 += 3 ;
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    v32 = _mm_srai_epi32(v32, 8) ;    // right shift by 8 positions
    _mm_storeu_si128((__m128i*) i32, v32) ;
    i32 += 4 ;
    u24 += 3 ;
  }
  if(i < n-3) {
    v24 = _mm_loadu_si128((__m128i*) u24) ;
    v32 = _mm_shuffle_epi8(v24, vx) ;
    v32 = _mm_srai_epi32(v32, 8) ;    // right shift by 8 positions
    _mm_storeu_si128((__m128i*) i32, v32) ;
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
int verify_upper(uint32_t *a, uint32_t *b, int n);
int verify_lower(uint32_t *a, uint32_t *b, int n);
int verify_upper(uint32_t *a, uint32_t *b, int n){
  int i, errors ;
  uint32_t mask = 0xFFFFFF00 ;   // upper 24 bits
  errors = 0 ;
  for(i=0 ; i<n ; i++) {
    if( (a[i] & mask) != (b[i] & mask) ) errors++ ;
  }
  return errors ;
}

int verify_lower(uint32_t *a, uint32_t *b, int n){
  int i, errors ;
  uint32_t mask = 0x00FFFFFF ;   // upper 24 bits
  errors = 0 ;
  for(i=0 ; i<n ; i++) {
    if( (a[i] & mask) != (b[i] & mask) ) errors++ ;
  }
  return errors ;
}

#define NPTS 4096
#define NTIMES 1000000

#include <misc_timers.h>

int main(int argc, char **argv){
  int i, j, nt, errors,  n = 7 ;
  static uint32_t fp32[8] = { 0x00810203, 0x04850607, 0x08890A0B, 0x0C8D0E0F, \
                              0x10911213, 0x14951617, 0x18991A1B, 0x1C9D1E1F  } ;
  static uint32_t uf24[8]  ;
  static uint32_t fp32b [8] ;
  uint32_t lfp32[NPTS], luf24[NPTS], lfp32b[NPTS] ;
  float float32[NPTS], float32b[NPTS] ;
  uint64_t t0, t1, tmin, tmax, tavg, freq ;
  double nano ;
  uint64_t t[NTIMES] ;

  freq = cycles_counter_freq() ;
  nano = 1000000000 ;
  nano /= freq ;
  for(i=0 ; i<NPTS ; i++) { lfp32[i] = i ; luf24[i] = 0 ; lfp32b[i] = 0 ; }
  for(i=0 ; i<NPTS ; i++) { float32[i] = -i - .111 ; float32b[i] = -1 ; }

  memset(uf24, 0xFF, sizeof(uf24));
  fp32_bf24(fp32, uf24, n) ;
  printf("fp32      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    fp32_bf24(float32, luf24, NPTS) ;
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
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %5.2f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("uf24      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;
  memset(fp32b, 0xFF, sizeof(fp32b));
  bf24_fp32(fp32b, uf24, n) ;
  errors = verify_upper((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;
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
  errors = verify_upper((uint32_t *) float32, (uint32_t *) float32b, NPTS) ;
  printf("floats : %f %f, relative error = %8.3g\n", 
         float32[NPTS/2], float32b[NPTS/2], (float32[NPTS/2]-float32b[NPTS/2])/float32[NPTS/2]) ;
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %5.2f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;
  printf("\n") ;

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
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %5.2f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("u24       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;
  memset(fp32b, 0xFF, sizeof(fp32b));
  u24u32(fp32b, uf24, n) ;
  errors = verify_lower((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
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
  errors = verify_lower((uint32_t *) lfp32, (uint32_t *) lfp32b, NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %5.2f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;

  printf("\n") ;

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
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %5.2f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("i24       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;
  memset(fp32b, 0xFF, sizeof(fp32b));
  i24i32((int32_t *) fp32b, uf24, n) ;
  errors = verify_lower((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
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
  errors = verify_lower((uint32_t *) lfp32, (uint32_t *) lfp32b, NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%6.1f [avg = %6.1f, %5.2f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;

  printf("\n\n\n\n\n");
}
#endif
