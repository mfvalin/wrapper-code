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

#define STATIC extern
#include <rmn/misc_types.h>
#include <rmn/misc_operators.h>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// ieee_prop struct to block header
extern inline uint16_t encode_ieee_header(ieee_prop prop){
  uint16_t header ;   // MSB -> LSB emax:8, full:1, same_exp:1, npti:3, nptj:3          (short block)       (full = 0)
                      //            emax:8, full:1, same_exp:1, allp:1, allm:1, spare:4 (full 8x8 block)    (same_exp = 0)
                      //            emax:8, full:1, same_exp:1, sign:1, cbit:3, nbit:2                      (full = 1, same_exp = 1)

  header = prop.emax << 8 ;             // common part
  header |= (prop.mima << 6) ;
  if(prop.npti == 8 && prop.nptj == 8){ // full block (npti * nptj == 64)
    header |= (1 << 7) ;                // full block flag
    if(prop.mima){
      header |= (prop.allm << 5) ;      // sign is 1 if all negative
      uint32_t ncbits = prop.emin >> 4 ;
      uint32_t cbits  = prop.emin & 0x7 ;
      header |= (cbits << 2) ;          // cbit field
      header |= ncbits ;                // nbit field
// printf("encode_ieee_header : ncbits = %d, cbits = %1x, header = %4.4x(%4.4x)\n", ncbits, cbits, header, header & 0x1F) ;
    }else{
      header |= (prop.allp << 5) ;      // not same exponent, keep allp and allm flags
      header |= (prop.allm << 4) ;
    }
  }else{                                // short block (npti * nptj < 64)
    header |= ((prop.npti - 1) << 3) ;
    header |= (prop.nptj - 1) ;
  }
  return header ;
}

// block header to ieee_prop struct
extern inline ieee_prop decode_ieee_header(uint16_t header){
  ieee_prop prop ;
  int full ;

  prop.emax = header >> 8 ;            // common part
  prop.emin = 0 ;                      // unknown
  prop.mima = (1 & (header >> 6)) ;
  prop.errf = 0 ;
  prop.xtra = 0 ;
  prop.zero = 0 ;                      // don't care
  full = (1 & (header >> 7)) ;
  if(full){                            // full block (npti * nptj == 64)
    prop.npti = 8 ;
    prop.nptj = 8 ;
    if(prop.mima){
      prop.allm = (1 & (header >> 5)) ;   // sign field
      prop.allp = 1 - prop.allm ;         // opposite of allm
      uint32_t ncbits =  header & 3 ;
      uint32_t cbits  =  (header >> 2) & 7 ;
      prop.emin = (ncbits << 4) | cbits ;
// printf("decode_ieee_header : ncbits = %d, cbits = %1x, header = %4.4x(%4.4x)\n", ncbits, cbits, header, header & 0x1F) ;
    }else{
      prop.allp = (1 & (header >> 5)) ;
      prop.allm = (1 & (header >> 4)) ;
    }
  }else{                                // short block (npti * nptj < 64)
    prop.npti = ((header >> 3) & 7) + 1 ;
    prop.nptj = (header & 7) + 1 ;
    prop.allp = 0 ;                     // always false for short blocks
    prop.allm = 0 ;
  }
  return prop ;
}

// encode properties into ieee_prop sruct
// allp  : logical OR reduction of all floats    (sign bit 0 only if all numbers are non negative)
// allm  : logical AND  reduction of all floats  (gign bit 1 only if all numbers are negative)
// emin  : exponent of the smallest non zero absolute value
// emax  : exponent of the largest absolute value
// zero  : non zero if no zero value was present
static inline ieee_prop make_ieee_prop(uint32_t allp, uint32_t allm, uint32_t emin, uint32_t emax, uint32_t zero){
  ieee_prop prop;
  uint32_t cbits, ncbits ;

  cbits = (allp ^ allm) << 9 ;                   // look for most significant common bits in mantissa
  ncbits = lzcnt_32(cbits) ;                     // number of identical most significant  bits in mantissa
  ncbits = (ncbits > 3) ? 3 : ncbits ;           // max allowed is 3
  cbits = (emax << 8) >> (32 - ncbits) ;         // extract common bits at top of mantissa
  allp >>= 31 ;  allp = allp ? 0 : 1 ;           // 1 if all numbers are non negative
  allm >>= 31 ;                                  // 1 if all numbers are negative
  prop.emax = emax >> 24 ;                       // largest exponent WITH IEEE bias (127)
  prop.emin = emin >> 24 ;                       // smallest exponent WITH IEEE bias (127)
  prop.allp = allp ;                             // 1 if all numbers are non negative
  prop.allm = allm ;                             // 1 if all numbers are negative
  prop.zero = zero ? 0 : 1 ;                     // 1 if zero value detected
  prop.mima = ((prop.emax == prop.emin) && (allp | allm) && (zero)) ? 1 : 0 ;    // same exponent and sign throughout
  if(prop.mima == 0){ ncbits = 0 ; cbits = 0 ; } // common mantissa bits used only if same exponent throughout and full block (8x8)
  else { prop.emin = (ncbits << 4) | cbits ;} ;  // store cbits and ncbits in emin field if same exponent throughout
  prop.xtra = 0 ;                                // reserved for future use, MUST be zero for now
  prop.errf = 0 ;                                // no error so far
  prop.npti = 0 ;                                // unknown at this point, will be inserted later
  prop.nptj = 0 ;                                // unknown at this point, will be inserted later
// printf("make_ieee_prop : emax = %d, emin = %d, mima = %d, allp = %d, allm = %d, zero = %d\n", 
//        prop.emax, prop.emin, prop.mima, prop.allp, prop.allm, prop.zero) ;
// if(ncbits > 0) printf("ncbits = %d, (%4.4x)\n", ncbits, cbits) ;
  return prop ;
}
// get the IEEE exponent of the largest float (absolute value)
// and the smallest non zero float (absolute value)  in float array f
// if all values of f are 0.0, 255 will be returned for the minimum exponent
// get sign properties (all >-0 , all <0)
ieee_prop ieee_properties(float *f, int n){
  int i ;
  uint32_t *uf = (uint32_t *) f ;
  uint32_t t, emax = 0 , emin = 0xFFu << 24, zero = 0xFFu << 24 ;
  uint32_t allp = 0, allm = 0xFFFFFFFFu ;
  ieee_prop prop;

  if(n > 64 || n < 0){
    printf("ieee_properties ERROR : invalid number of points : expected 1 -> 64, got %d\n", n) ;
    prop.errf = 1 ;
    return prop ;
  }
  for(i=0 ; i<n ; i++){
    allp |= uf[i] ;                   // upper bit will remain 0 if >= 0 floats only
    allm &= uf[i] ;                   // upper bit will remain 1 if < 0 floats only
    t = uf[i] << 1 ;                  // get rid of sign, exponent in upper 8 bits
    zero = (t < zero) ? t : zero ;    // smallest value, including zero
    emax = (t > emax) ? t : emax ;
    t = (t == 0) ? (0xFFu << 24) : t ;         // ignore values of 0
    emin = (t < emin) ? t : emin ;
  }
  prop = make_ieee_prop(allp, allm, emin, emax, zero) ;
  if(n != 64) {                           // short block, only flag kept is mima (same exponent)
    if( prop.allp != 1) prop.mima = 0 ;   // and only if all numbers are non negative
    prop.allp = prop.allm = 0 ;           // suppress allm and allp flags
  }else{
    prop.npti = prop.nptj = 8 ;          // full 8x8 block
  }
  return prop ;
}
ieee_prop ieee_properties_64(float *f){  // special case (frequent occurrance) used for 8x8 blocks
  int i ;
  uint32_t *uf = (uint32_t *) f ;
  uint32_t t, emax = 0 , emin = 0xFFu << 24, zero = 0xFFu << 24 ;
  uint32_t allp = 0, allm = 0xFFFFFFFFu ;
  ieee_prop prop;

  for(i=0 ; i<64 ; i++){
    allp |= uf[i] ;                   // upper bit will remain 0 if >= 0 floats only
    allm &= uf[i] ;                   // upper bit will remain 1 if < 0 floats only
    t = uf[i] << 1 ;                  // get rid of sign, exponent in upper 8 bits
    zero = (t < zero) ? t : zero ;    // smallest absolute value, including zero
    emax = (t > emax) ? t : emax ;    // largest absolute value
    t = (t == 0) ? (0xFFu << 24) : t ;         // ignore values of 0 for emin
    emin = (t < emin) ? t : emin ;    // smallest absolute value (larger than 0)
  }
  prop = make_ieee_prop(allp, allm, emin, emax, zero) ;
  prop.npti = prop.nptj = 8 ;         // full 8x8 block
// printf("ieee_properties_64 mima = %d, allp = %d, allm = %d\n", prop.mima, prop.allp, prop.allm) ;
  return prop ;
}

// encode IEEE float block as a sequence of 16 bit tokens
ieee_prop ieee_encode_block_16(float xf[64], int ni, int nj, uint16_t *restrict stream){
  int i, n = ni*nj ;
  uint32_t *xi = (uint32_t *) xf ;
  uint16_t header  ;
  FloatInt factor ;
  int16_t *t = (int16_t *) stream ;
  int32_t tmp ;
  uint32_t utmp ;
  ieee_prop prop ;

  if(ni < 1 || ni > 8 || nj < 1 || nj > 8){
    printf("ieee_encode_block_16 ERROR : invalid number of points : expected (1-8) x (1-8), got %d x %d\n", ni, nj) ;
    prop.errf = 1 ;
    return prop ;
  }
  prop = (ni == 8 && nj == 8) ? ieee_properties_64(xf) : ieee_properties(xf, ni*nj) ;
  prop.npti = ni ;   // add block dimensions
  prop.nptj = nj ;
  stream[0] = encode_ieee_header(prop) ;
// printf("ieee_encode_block_16 : emax = %d, emin = %d, mima = %d, allp = %d, allm = %d, zero = %d\n", 
//        prop.emax, prop.emin, prop.mima, prop.allp, prop.allm, prop.zero) ;

//   allp and allm CANNOT be set if ni*nj != 64
  if(prop.mima) {                                                   // same max exponent, same sign, no zero
    uint32_t ncbits = (n == 64) ? (prop.emin >> 4) : 0 ;
// printf("encode : mima, emin = %2.2x, ncbits = %d ", prop.emin, ncbits) ;
    for(i=0 ; i<n ; i++){
      utmp = xi[i] << ncbits ;                                      // push common bits left
      utmp = (utmp & 0x7FFFFF) + 0x40 ;                             // mantissa after rounding
      utmp >>= 7 ;                                                  // keep upper 16 bits of mantissa (rounded)
      utmp = (utmp > 0xFFFF) ? 0xFFFF : utmp ;                      // clipped at 0xFFFF
      stream[i+1] = utmp ;
    }
  }else if(prop.allp){                                              // all numbers are non negative
    factor.i = (127 + 142 - prop.emax) << 23 ;                      // largest number will be 2**16 - 1
// printf("encode : allp factor = %f ", factor.f) ;
    for(i=0 ; i<n ; i++){
      tmp = factor.f * xf[i] + 0.5f ;
      tmp = (tmp > 0xFFFF) ? 0xFFFF : tmp ;
      stream[i+1] = tmp ;
    }
  }else if(prop.allm){                                              // all numbers negative
    factor.i = (127 + 142 - prop.emax) << 23 ;                      // largest number will be 2**16 - 1
    factor.f = -factor.f ;
// printf("encode : allm factor = %f ", factor.f) ;
    for(i=0 ; i<n ; i++){
      tmp = factor.f * xf[i] + 0.5f ;
      tmp = (tmp > 0xFFFF) ? 0xFFFF : tmp ;
      stream[i+1] = tmp ;
    }
  }else{                                                            // range of exponents, mixed signs
    factor.i = (127 + 141 - prop.emax) << 23 ;                      // largest number will be 2**15 - 1
// printf("encode : factor = %f ", factor.f) ;
    for(i=0 ; i<n ; i++){
      tmp = factor.f * xf[i] + ((xf[i] < 0) ? -0.5f : 0.5f) ;
      tmp = (tmp > 0x7FFF) ? 0x7FFF : tmp ;
//       tmp = (tmp < -0x7FFF) ? -0x7FFF : tmp ;
      t[i+1] = tmp ;
    }
  }
// printf(" ni = %d(%d), nj = %d(%d)\n", prop.npti, ni, prop.nptj, nj) ;
  return prop ;
}

// encode a float array
void ieee_encode_array_16(float *restrict xf, int ni, int lni, int nj, uint16_t *restrict stream){
  int i, j, indx, nil, njl ;
  ieee_prop prop ;
  float x64[64] ;
  for(j=0 ; j<nj ; j+= 8){
    njl = (nj-j < 8) ? (nj - j) : 8 ;
    for(i=0 ; i<ni ; i+=8){
      indx = i + j * lni ;                                  // index of lower left corner of block to encode
      nil = (ni-i < 8) ? (ni - i) : 8 ;
      prop = ieee_get_block(&xf[indx], x64, 8, nil, 8) ;    // get nil x njl block
      prop = ieee_encode_block_16(x64, nil, njl, stream) ;  // encode block
      stream += 1 + nil * njl ;                             // space taken by previous block
    }
  }
}

// decode IEEE float block as a sequence of 16 bit tokens
// xf 
ieee_prop ieee_decode_block_16(float xf[64], int ni, int nj, uint16_t *restrict stream){
  int i, n = ni*nj ;
  uint32_t *xi = (uint32_t *) xf ;
  FloatInt factor, fi ;
  int16_t *t = (int16_t *) stream ;
  int32_t tmp ;
  uint32_t sign = 0 ;
  ieee_prop prop ;

  prop = decode_ieee_header(stream[0]) ;
  if(prop.npti  != ni || prop.nptj != nj){
    for(i=0 ; i<n ; i++) xi[i] = 0xFF8 << 19 ;         // NaN
    printf("ieee_decode_block_16 : ERROR, inconsistent dimensions : expected %d x %d , got %d x %d\n", prop.npti, prop.nptj, ni, nj) ;
    prop.errf = 1 ;   // error flag
    return prop ;
  }
// printf("decode : ni = %d(%d), nj = %d(%d), allm = %d, allp = %d, mima = %d\n", 
//        prop.npti, ni, prop.nptj, nj, prop.allm, prop.allp, prop0.mima) ;
  if(prop.mima) {                                                   // same max exponent, same sign
    if(prop.allm) sign = 1 << 31 ;
    uint32_t ncbits = (n == 64) ? (prop.emin >> 4) : 0 ;
    uint32_t cbits  = (n == 64) ? (prop.emin &  7) : 0 ;
// printf("decode mima sign = %8.8x, emin = %2.2x, ncbits = %d, cbits = %4.4x\n", sign, prop.emin, ncbits, cbits) ;
    for(i=0 ; i<n ; i++){
      tmp = stream[i+1] | (cbits << 16) ;                           // insert common bits
      tmp <<= (7 - ncbits) ;                                        // align mantissa to proper bit position
      fi.i = (prop.emax << 23) | tmp | sign ;                       // restore upper 16 bits of mantissa
      xf[i] = fi.f ;
    }
  }else if(prop.allp){                                              // all numbers positive
    factor.i = (127 - 142 + prop.emax) << 23 ;                      // largest number will be 2**16 - 1
//  printf("decode allp factor = %f\n", factor.f) ;
    for(i=0 ; i<n ; i++){
      xf[i] = factor.f * stream[i+1] ;
    }
  }else if(prop.allm){                                              // all numbers negative
    factor.i = (127 - 142 + prop.emax) << 23 ;                      // largest number will be 2**16 - 1
    factor.f = -factor.f ;
//  printf("decode allm factor = %f\n", factor.f) ;
    for(i=0 ; i<n ; i++){
      xf[i] = factor.f * stream[i+1] ;
    }
  }else{                                                            // range of exponents, mixed signs
    factor.i = (127 - 141 + prop.emax) << 23 ;                      // largest number will be 2**15 - 1
// printf("decode factor = %f\n", factor.f) ;
    for(i=0 ; i<n ; i++){
      xf[i] = factor.f * t[i+1] ;
    }
  }
  return prop ;
}

// copy a block of (1 <= ni <= 8) x (1 <= nj <=8) 32 bit words from array blk back into f
// f     : destination array
// ni    : row size
// nj    : number of rows
// lni   : storage length of rows
// blk   : 32 bit words array to provide copied block
void put_w32_block(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict d = (uint32_t *) f ;
  uint32_t *restrict s = (uint32_t *) blk ;
  int i, j ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vm = _mm256_memmask_si256(ni) ;  // mask for load and store operations
  if(nj == 8){
    if(ni == 8){
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += lni ; s += ni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += lni ; s += ni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += lni ; s += ni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += lni ; s += ni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += lni ; s += ni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += lni ; s += ni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += lni ; s += ni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += lni ; s += ni ;
    }else{                // ni < 8
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += lni ; s += ni ;
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += lni ; s += ni ;
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += lni ; s += ni ;
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += lni ; s += ni ;
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += lni ; s += ni ;
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += lni ; s += ni ;
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += lni ; s += ni ;
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += lni ; s += ni ;
    }
    return ;
  }else{                   // nj < 8
    for(j=0 ; j<nj ; j++){
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ;
      d += lni ; s += ni ;
    }
  }
#else

  if(nj == 8){
    if(ni == 8){
      for(j=0 ; j<8  ; j++){ for(i=0 ; i<8 ; i++) d[i] = s[i] ; d += lni ; s += 8 ; }
    }else{
      ni &= 7 ;
      for(j=0 ; j<8  ; j++){ for(i=0 ; i<ni; i++) d[i] = s[i] ; d += lni ; s += ni ; }
    }
  }else{
    if(ni == 8){
      for(j=0 ; j<nj ; j++){ for(i=0 ; i<8 ; i++) d[i] = s[i] ; d += lni ; s += 8 ; }
    }else{
      ni &= 7 ;
      for(j=0 ; j<nj ; j++){ for(i=0 ; i<ni; i++) d[i] = s[i] ; d += lni ; s += ni ; }
    }
  }
#endif
}

// copy a block of (1 <= ni <= 8) x (1 <= nj <=8) 32 bit words into array blk
// return ieee properties (analysis on the fly)
ieee_prop get_ieee32_block(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i, j ;
  uint32_t emax, emin, allm, allp ;
  ieee_prop prop ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vo0, vma, vmi, va0, vor, vand ;
  __m128i t ;
  vo0  = _mm256_loadu_si256((__m256i *)s) ; s += lni ;     // get row
  _mm256_storeu_si256((__m256i *)d, vo0) ; d += 8 ;        // store row
  vor  = vo0 ;                                             // wired OR accumulator, MSB 0 only if all numbers non negative
  vand = vo0 ;                                             // wired AND accumulator, MSB 1 only if all numbers negative
  va0  = _mm256_slli_si256(vo0, 1) ;                       // get rid of sign, exponent in upper 8 bits
  vma  = va0 ;                                             // largest exponent
  vmi  = va0 ;                                             // smallest exponent
  for(j=1 ; j<8 ; j++){
    vo0  = _mm256_loadu_si256((__m256i *)s) ; s += lni ;
    vand = _mm256_and_si256(vand, vo0) ;                   // wired AND, MSB 1 if all numbers negative
    vor  = _mm256_or_si256(vor, vo0) ;                     // wired OR, MSB 0 if all numbers non negative
    va0  = _mm256_slli_si256(vo0, 1) ;                     // get rid of sign, exponent in upper 8 bits
    vma  = _mm256_max_epu32(vma, va0) ;                    // largest exponent
    vmi  = _mm256_min_epu32(vmi, va0) ;                    // smallest exponent
    _mm256_storeu_si256((__m256i *)d, vo0) ; d += 8 ;
  }
  t = _mm_max_epu32(_mm_lower128(vma) , _mm_upper128(vma) ) ;  // 8 -> 4 max( 0,4 1,5 2,6 3,7 )
  t = _mm_max_epu32(t, _mm_shuffle_epi32(t, 0b11101110) ) ;    // 4 -> 2 max( 0,4,2,6 1,5,3,7 xxxxx xxxxx)
  t = _mm_max_epu32(t, _mm_shuffle_epi32(t, 0b01010101) ) ;    // 2 -> 1 max( 0,4,2,6,1,5,3,7 xxxxx xxxxx xxxxx)
  _mm_storeu_si32(&emax, t) ;                                  // store reduced value
  t = _mm_min_epu32(_mm_lower128(vmi) , _mm_upper128(vmi) ) ;
  t = _mm_min_epu32(t, _mm_shuffle_epi32(t, 0b10111011) ) ;
  t = _mm_min_epu32(t, _mm_shuffle_epi32(t, 0b01010101) ) ;
  _mm_storeu_si32(&emin, t) ;                                  // store reduced value
  prop.allp = ( _mm256_movemask_ps((__m256) vor) == 0 ) ;      // all upper bits MUST be 0
  prop.allm = ( _mm256_movemask_ps((__m256) vand) == 0xFF ) ;  // all upper bits MUST be 1
#else
#endif
  prop.emin = emin ;
  prop.mima = (emax == emin) ;
  prop.emax = emax ;
  return prop ;
}

// copy a block of (1 <= ni <= 8) x (1 <= nj <=8) 32 bit words into array blk
// f     : source array
// ni    : row size
// nj    : number of rows
// lni   : storage length of rows
// blk   : 32 bit words array to receive copied block
// N.B. this code is not entirely safe as it may "overread" from array f by up to 3 locations
void get_w32_block(void *restrict f, void *restrict blk, int ni, int lni, int nj){
  uint32_t *restrict s = (uint32_t *) f ;
  uint32_t *restrict d = (uint32_t *) blk ;
  int i, j ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vm = _mm256_memmask_si256(ni) ;  // mask for load and store operations
  if(nj == 8){
    if(ni == 8){
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += ni ; s += lni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += ni ; s += lni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += ni ; s += lni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += ni ; s += lni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += ni ; s += lni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += ni ; s += lni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += ni ; s += lni ;
      _mm256_storeu_si256((__m256i *)d, _mm256_loadu_si256((__m256i *)s)) ; d += ni ; s += lni ;
    }else{   // ni < 8, use a masked load and store
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += ni ; s += lni ;
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += ni ; s += lni ;
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += ni ; s += lni ;
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += ni ; s += lni ;
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += ni ; s += lni ;
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += ni ; s += lni ;
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ; d += ni ; s += lni ;
    _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ;
    }
  }else{    // nj < 8
    for(j=0 ; j<nj ; j++){
      _mm256_maskstore_epi32 ((int *)d, vm, _mm256_maskload_epi32 ((int const *) s, vm) ) ;
      d += ni ; s += lni ;
    }
  }
  return ;
#else

  if(ni > 4){            // possible "overread" by up to 3 elements
    if(nj & 8) {
      for(j=0 ; j<8 ; j++){ for(i=0 ; i<8; i++) d[i] = s[i] ; d += ni ; s += lni ; }
      return ;
    }
    if(nj & 4) for(j=0 ; j<4 ; j++){ for(i=0 ; i<8; i++) d[i] = s[i] ; d += ni ; s += lni ; }
    if(nj & 2) for(j=0 ; j<2 ; j++){ for(i=0 ; i<8; i++) d[i] = s[i] ; d += ni ; s += lni ; }
    if(nj & 1) for(i=0 ; i<8; i++) d[i] = s[i] ;
    return ;
  }else{                 // possible "overread" by up to 3 elements
    if(nj & 8) for(j=0 ; j<8 ; j++){ for(i=0 ; i<4; i++) d[i] = s[i] ; d += ni ; s += lni ; }
    if(nj & 4) for(j=0 ; j<4 ; j++){ for(i=0 ; i<4; i++) d[i] = s[i] ; d += ni ; s += lni ; }
    if(nj & 2) for(j=0 ; j<2 ; j++){ for(i=0 ; i<4; i++) d[i] = s[i] ; d += ni ; s += lni ; }
    if(nj & 1)                       for(i=0 ; i<4; i++) d[i] = s[i] ;
    return ;
  }
#endif
}
// copy a block of (1 <= ni <= 8) x (1 <= nj <=8) floats into array blk
// f     : float array
// ni    : row size
// nj    : number of rows
// lni   : storage length of rows
// blk   : float array to receive copied block
// return ieee properties for block 
// N.B.  : Fortran storage order is assumed
ieee_prop ieee_get_block(float *restrict f, float blk[64], int ni, int lni, int nj){
  int i, j ;
  float *blk0 = blk ;
  ieee_prop prop ;

  if(ni < 1 || nj < 1 || ni > 8 || nj > 8){
    printf("ieee_get_block ERROR : invalid number of points : expected 1 -> 64, got %d\n", ni*nj) ;
    prop.errf = 1 ;
    return prop ;
  }
  get_w32_block(f, blk, ni, lni, nj) ;
  if(ni == 8 && nj == 8){      // full 8x8 block
    prop = ieee_properties_64(blk0) ;
  }else{                       // short block
    prop = ieee_properties(blk0, ni*nj) ;
  }
  prop.npti = ni ;
  prop.nptj = nj ;
// printf("ieee_get_block mima = %d, allp = %d, allm = %d\n", prop.mima, prop.allp, prop.allm) ;
  return prop ;
}

// get the IEEE exponent of the largest float (absolute value)
// and the smallest non zero float (absolute value)  in float array f
// if all values of f are 0.0, 255 will be returned for the minimum exponent
// uint32_t ieee_minmax_exponent(float *f, int n){
//   int i ;
//   uint32_t *uf = (uint32_t *) f ;
//   uint32_t t, emax = 0 , emin = 0xFF << 24 ;
// 
//   for(i=0 ; i<n ; i++){
//     t = uf[i] << 1 ;                  // get rid of sign, exponent in upper 8 bits
//     emax = (t > emax) ? t : emax ;
//     t = (t == 0) ? emin : t ;         // ignore values of 0
//     emin = (t < emin) ? t : emin ;
//   }
//   emax >>= 24 ;                       // exponent WITH IEEE bias (127)
//   emin >>= 24 ;                       // exponent WITH IEEE bias (127)
//   return emax | (emin << 8) ;         // emax and emin
// }

// get the IEEE exponent of the largest float (absolute value) in float array f
// n : number of values in array f
uint32_t ieee_max_exponent(float *f, int n){
  int i ;
  uint32_t *uf = (uint32_t *) f ;
  uint32_t t, emax = 0 ;

  for(i=0 ; i<n ; i++){
    t = uf[i] << 1 ;                  // get rid of sign, exponent in upper 8 bits
    emax = (t > emax) ? t : emax ;
  }
  emax >>= 24 ;                       // exponent WITH IEEE bias (127)
  return emax ;
}

// get the IEEE exponent of the smallest non zero float (absolute value) in float array f
// if all values of f are 0.0, 255 will be returned
// n : number of values in array f
uint32_t ieee_min_exponent(float *f, int n){
  int i ;
  uint32_t *uf = (uint32_t *) f ;
  uint32_t t, emin = 0xFF << 24 ;     // largest possible exponent

  for(i=0 ; i<n ; i++){
    t = uf[i] << 1 ;                  // get rid of sign, exponent in upper 8 bits
    t = (t == 0) ? emin : t ;         // ignore values of 0
    emin = (t < emin) ? t : emin ;
  }
  emin >>= 24 ;                       // exponent WITH IEEE bias (127)
  return emin ;
}

// vector version of BitsNeeded_u32
// what      : array of unsigned 32 bit integers
// nn        : ABS(nn) = number of values in array what
//             if nn < 0, do not update table bits with population counts
// bits      : cumulative bits needed count
int32_t vBitsNeeded_u32(uint32_t * restrict what, int32_t * restrict bits, int nn){
  int i, needed=0 ;
  int n = (nn < 0) ? -nn : nn ;
  for(i=0 ; i<n ; i++) {
    int32_t nbits = BitsNeeded_u32(what[i]) ;
    if(nn >= 0) bits[nbits]++ ;
    needed = (nbits > needed) ? nbits : needed ;
  }
  return needed ;
}

// vector version of BitsNeeded_32
// what      : array of signed 32 bit integers
// nn        : ABS(nn) = number of values in array what
//             if nn < 0, do not update table bits with population counts
// bits      : cumulative bits needed count
int32_t vBitsNeeded_32(int32_t * restrict what, int32_t * restrict bits, int nn){
  int i, needed=0 ;
  int n = (nn < 0) ? -nn : nn ;
  for(i=0 ; i<n ; i++) {
    int32_t nbits = BitsNeeded_u32(what[i]) ;
    if(nn >= 0) bits[nbits]++ ;
    needed = (nbits > needed) ? nbits : needed ;
  }
  return needed ;
}

// vector version of BitsNeeded_u64
// what      : array of unsigned 64 bit integers
// nn        : ABS(nn) = number of values in array what
//             if nn < 0, do not update table bits with population counts
// bits      : cumulative bits needed count
int32_t vBitsNeeded_u64(uint64_t * restrict what, int32_t * restrict bits, int nn){
  int i, needed=0 ;
  int n = (nn < 0) ? -nn : nn ;
  for(i=0 ; i<n ; i++) {
    int32_t nbits = BitsNeeded_u64(what[i]) ;
    if(nn >= 0) bits[nbits]++ ;
    needed = (nbits > needed) ? nbits : needed ;
  }
  return needed ;
}

// vector version of BitsNeeded_64
// what      : array of signed 64 bit integers
// nn        : ABS(nn) = number of values in array what
//             if nn < 0, do not update table bits with population counts
// bits      : cumulative bits needed count
int32_t vBitsNeeded_64(int64_t * restrict what, int32_t * restrict bits, int nn){
  int i, needed=0 ;
  int n = (nn < 0) ? -nn : nn ;
  for(i=0 ; i<n ; i++) {
    int32_t nbits = BitsNeeded_64(what[i]) ;
    if(nn >= 0) bits[nbits]++ ;
    needed = (nbits > needed) ? nbits : needed ;
  }
  return needed ;
}

void BitEntropy4(float entropy[4], uint32_t *bitstream, int npts, int nbits, int rshift) {
  uint32_t bins[4][16] ;
  uint32_t mask = 0xF ;
  int i, j ;
  uint32_t t ;
  float e[4] ;
  float prob, scale, ovlog2 ;

  scale = 1.0f / (float)(npts) ;
  ovlog2 = 1.0f / logf(2.0f) ;
  for(j=0 ; j<4 ; j++) for(i=0 ; i<16 ; i++) bins[j][i] = 0 ;
  for(j=0 ; j<4 ; j++) e[j] = 0.0f ;
  for (i = 0; i < npts; i++) {
    t = bitstream[i] ;
    bins[0][t & mask]++ ; t >>= 4 ;
    bins[1][t & mask]++ ; t >>= 4 ;
    bins[2][t & mask]++ ; t >>= 4 ;
    bins[3][t & mask]++ ;
  }
  for(j=0 ; j<4 ; j++){
    for(i=0 ; i<16 ; i++){
      if (bins[j][i] != 0) {
        prob = (float)bins[j][i] * scale ;
        e[j] += (prob * logf(prob) * ovlog2) ;
      }
    }
  }
  for(j=0 ; j<4 ; j++) entropy[j] = (-e[j]) ;
}

float BitEntropy(int32_t *bitstream, int npts, int nbits, int rshift) {
    int lbin, sizebins;
//     int imin, imax, range, nbits_local ;
    float prob, entropie, scale, ovlog2 ;
    uint32_t *bins;
    uint32_t mask ;
    uint32_t samples = 0 ;

    entropie = 0.0f;
    mask = RMASK32(nbits) ;
    // printf("imin : %d imax : %d range: %d, nbits : %d\n", imin, imax, imax-imin, nbits_local);

    sizebins = 1 << (nbits-rshift) ;
    bins = (uint32_t *) calloc(sizebins,sizeof(uint32_t));

    for (int i = 0; i < npts; i++) {
        lbin = (bitstream[i] & mask) >> rshift ;
        bins[lbin]++;
    }

    scale = 1.0f / (float)(npts) ;
    ovlog2 = 1.0f / logf(2.0f) ;
    for (int i = 0; i < sizebins; i++) {
        if (bins[i] != 0) {
            prob = (float)(bins[i]) * scale ;
            // printf("i: %d count: %d prob: %f contrib : %f \n", i, bins[i], prob, (prob * log(prob)/log(2.0)));
            entropie += prob * logf(prob) * ovlog2 ;
            samples++ ;
        }
    }

    free(bins) ;
// printf("sizebins = %d, mask = %8.8x, rshift = %d, samples = %d, entropie = %g\n", sizebins, mask, rshift, samples, entropie);
    return -entropie ;
}

