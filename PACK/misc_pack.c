/* Hopefully useful routines for C and FORTRAN
 * Copyright (C) 2020  Recherche en Prevision Numerique
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#include <stdio.h>
#include <stdint.h>

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
#include <immintrin.h>
#endif

#include <misc_types.h>
#include <misc_operators.h>
#include <misc_pack.h>

#if ! defined(__INTEL_COMPILER_UPDATE)
#pragma GCC optimize "tree-vectorize"
#endif

// adjust quantum to first power of 2 <= quantum
float quantum_adjust(float quantum){
  FloatInt rq ;

  rq.f = quantum ;
  rq.i &= 0x7F800000 ;   // drop mantissa bits
//   printf("quantum = %g, adjusted quantum = %g\n", quantum, rq.f) ;
  return rq.f ;
}

// pre-quantizer, prepare the packing header that will be used by the next stages
// p       : packing header to be initialized                         [OUT]
// nbits   : quantize using at most nbits bits (assumed to be <= 24)  [IN]
// maxval  : highest value in array  to be quantized (float)          [IN]
// minval  : lowest value in array  to be quantized (float)           [IN]
// quantum : if nonzero and positive,                                 [IN]
//           use inverse of quantum as the quantization factor
//           nbits will then be recomputed internally as a function of quantum
//           quantum will be adjusted to the first power of 2 <= quantum
void float_quantize_prep(int nbits, QuantizeHeader *p, float maxval, float minval, float quantum) {
  FloatInt   m1, m2, m3;  // access as float or 32 bit int
  DoubleLong m0 ;         // access as double or 64 bit long
  int exp1, exp2, exp3 ;
  float fac32, range;
  double fac64 ;
  int mask_trunc ;
  int offset ;
  int64_t irange ;

  range = maxval - minval;                   // range of floating point numbers
  if(quantum > 0.0) {
    int neededbits ;
    m1.f = quantum_adjust(quantum) ;         // first power of 2 <= quantum
    irange = (range / m1.f) ;
    NEEDBITS( irange, neededbits) ;
    nbits = neededbits ;
  }

  if(nbits <= 0 ) return;
  if(nbits > 24 ) nbits = 24 ;               // no more than 24 bits will be kept
  p->nbits = nbits ;

  mask_trunc = ( -1 << (24 - nbits) ) ;      // truncation mask
  m1.f = maxval ;
  m2.f = minval ;
  m3.f = range ;
  exp1 = 0xFF & (m1.i >> 23) ;               // exponents for max, min, range
  exp2 = 0xFF & (m2.i >> 23) ;
  exp3 = 0xFF & (m3.i >> 23) ;
  exp1 = (exp2 > exp1) ? exp2 : exp1 ;       // largest exponent ( max, min ) (biased)
  exp1 = (exp3 > exp1) ? exp3 : exp1 ;       // largest exponent ( max, min, range ) (biased)
  p->e = exp1 - 127 ;                        // largest exponent (with IEEE bias removed)

  if(exp1 < 255 && exp1 > 0){                // O.K. for IEEE float 32 (most of the time)
    m1.i = (127 + (23 - (exp1 - 127))) ;     // factor to bring largest exponent to 23
    m1.i <<= 23 ;
    fac32 = m1.f;                            // normalizing multiplier
    offset = minval * fac32 ;                // quantized minimum value
//     printf(" fac32 = %10.4g(%d) ", fac32, p->e);
  }else{                                     // must use IEEE float 64 (rare cases)
    m0.l = (1023 + (23 - (exp1 - 127)));     // factor to bring largest exponent to 23
    m0.l = m0.l << 52 ;
    fac64 = m0.d;                            // normalizing multiplier
    offset = minval * fac64 ;                // quantized minimum value
//     printf(" fac64 = %10.4g(%d) ", fac64, p->e);
  }

  m3.i = (exp3 - nbits + 1) << 23 ;          // subtract nbits - 1 from range exponent (quick divide)
  if( (exp3 - nbits + 1) < 0) {              // range is too small, compute quantum the hard way
    m3.i = exp3 << 23 ;                      // exp3 is IEEE exponent from range
    m3.f = m3.f / (1 << (nbits -1)) ;        // divide by 2**(nbits - 1)
  }
  p->q = m3.f ;                              // effective quantum
  offset = offset & mask_trunc ;             // drop lower (24 - nbits) bits
  p->o = offset ;                            // truncated minimum value
}

// "linear" quantizer for 32 bit floats producing a stream of unsigned 32 bit integers
// nbits is assumed to be <= 24 and taken from the packing header
// iz      : quantized stream (32 bit elements)  [OUT]
// z       : float array to be quantized         [IN]
// ni      : number of useful points in rows     [IN]
// lni     : storage size of rows in array z     [IN]
// lniz    : storage size of rows in array iz    [IN]
// nj      : number of rows in z and iz          [IN]
// p       : packing header as initialized by float_quantize_prep  [IN]
void float_quantize(void *iz, float *z, int ni, int lni, int lniz, int nj, QuantizeHeader *p ) {
  int i, j, it;
  int offset = p->o ;   // offset reflecting minimum value
  int expmax = p->e ;   // largest exponent in original float values (bias removed)
  int round = 0 ;
  int roundp, roundm ;
  int *izw = (int   *) iz;
  int nbits = p->nbits ;

  if(nbits <= 0 ) return ;
  if(nbits > 24 ) nbits = 24 ;           // no more than 24 bits will be kept
  if(nbits < 24)                         // no rounding if nbits = 24
    round = 1 << (23 - nbits);           // rounding for quantized value
  roundp = round - offset ;              // rounding for positive values combined with offset
  roundm = (nbits < 24) ? roundp - 1 : roundp ;   // rounding for negative values combined with offset
  if(ni == lni && ni == lniz) {          // can we fuse loops ?
    ni = ni * nj ;
    nj = 1 ;
  }

  if(expmax > -127 && expmax < 127){       // "civilized" exponent for largest value
    FloatInt m1;                           // a float can be used to perform the quantification
    float fac32 ;
    m1.i = (127 + (23 - expmax)) << 23 ;   // factor to bring largest exponent to 23
    fac32 = m1.f ;
    for(j=0 ; j<nj ; j++){
      for (i=0 ; i<ni ; i++){
        // vz = _mm256_loadu_ps(z)             // load z[i:i+7]
        // t  = _mm256_mul_ps(vz , V(fac32))   // * fac32
        // it = _mm356_cvtps_epi32(t)          // convert to integer
        it = z[i] * fac32;                 // "normalize" to largest exponent
        // t  = _mm256_cmpgt_epi21(vz, V(0))   // where vz > 0 use roundp, else use roundm
        // r  = _mm256_blendv_epi8(roundm, roundp, t)
        round = (z[i] < 0) ? roundm : roundp ;   // different adjustment for positive and negative values
        // it = _mm256_add_epi32(it, r)        // add adjustment term and shift
        // it = _mm256_srlv_epi32(it, V(24-nbits))
        it = (it + round) >> (24-nbits);   // remove offset, add rounding term
        // it = _mm256_max_epi32(it, V(0))     // make sure result is >= 0
        // _mm256_storeu_si256(izw, it)        // store izw[i:i+7]
        izw[i] = (it < 0) ? 0 : it ;       // quantized result MUST be >= 0
      }
      z += lni ;
      izw += lniz ;
    }
  }else{                                   // need to use a double to perform the quantification correctly
    DoubleLong m1;
    double fac64 ;
    m1.l = (1023 + (23 - expmax)) ;        // factor to bring largest exponent to 23
    m1.l = m1.l << 52 ;
    fac64 = m1.d ;
    for(j=0 ; j<nj ; j++){
      for (i=0 ; i<ni ; i++){
        it = z[i] * fac64;                 // "normalize" to largest exponent
        round = (z[i] < 0) ? roundm : roundp ;   // different adjustment for positive and negative values
        it = (it + round) >> (24-nbits);   // remove offset, add rounding term
        izw[i] = (it < 0) ? 0 : it ;       // quantized result MUST be >= 0
      }
      z += lni ;
      izw += lniz ;
    }
  }
}

// quantize an array of floats (transform floats into integers using a quantization unit)
// z       : array of floats to quantize                 [IN]
// q       : integer array, quantized values             [OUT]
// ni      : number of values                            [IN]
// quantum : quantization unit                           [IN]
// t       : lowest, highest quantized values            [OUT]
// return  : number of bits needed to represent the lowest->highest quantized values
//           quantum will be internally adjusted to the first power of 2 <= quantum
// the vectorisers seem to have no problem with this function
uint32_t float_quantize_simple_1D(float * restrict z, int32_t * restrict q, int ni, float quantum, IntPair *t){
  return float_quantize_simple(z, q, ni, ni, ni, 1, quantum, t) ;
}

// 2D version of float_quantize_simple_1D
// lniz  : storage length of z rows      [IN]
// lniq  : storage length of q rows      [IN]
// nj    : number of rows                [IN]
uint32_t float_quantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, IntPair *t){
  int i, i0 ;
  int32_t min = 0x7FFFFFFF, max = -min ;
  float ovq ;
  FloatInt rq ;
  uint32_t needed1, needed2 ;

  rq.f = quantum ; rq.i &= 0x7F800000 ; quantum = rq.f ;   // adjust quantum to power of 2 <= quantum
  ovq = 1.0f / quantum ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 vz, vf ;
  __m256i vq, vmi, vma ;
  __m128i vma0, vma1, vmi0, vmi1 ;
  int n7 = (ni & 7) ? (ni & 7) : 8 ;
  uint32_t round_mode ;

  if(ni < 8) goto less_than_8 ;
  round_mode = _MM_GET_ROUNDING_MODE() ;
  _MM_SET_ROUNDING_MODE(_MM_ROUND_NEAREST) ;
  vf  = _mm256_set1_ps(ovq) ;
  vmi = _mm256_set1_epi32(min) ;
  vma = _mm256_set1_epi32(max) ;
  while(nj-- > 0) {
    // N.B. the first and second pass will process some values twice (overlap)
    for(i0 = 0 ; i0 < ni-7 ; i0+=n7 , n7=8){          // batches of 8 values
      vz = _mm256_loadu_ps(z+i0) ;                    // fetch z
      vz = _mm256_mul_ps(vz, vf) ;                    // * ovq
      vq = _mm256_cvtps_epi32(vz) ;                   // convert to integer (use current rounding mode)
      vma = _mm256_max_epi32(vma, vq) ;               // max(q, max)
      vmi = _mm256_min_epi32(vmi, vq) ;               // min(q, min)
      _mm256_storeu_si256( (__m256i *)(q+i0), vq) ;   // store q
    }
    z += lniz ;
    q += lniq ;
  }
  _MM_SET_ROUNDING_MODE(round_mode) ;
  vma0 = _mm256_extracti128_si256(vma, 0) ;     // fold max to 1 value
  vma1 = _mm256_extracti128_si256(vma, 1) ;
  vma0 = _mm_max_epi32(vma0, vma1) ;            // max( 0|4 , 1|5 , 2|6 , 1|7 )
  vma1 = _mm_shuffle_epi32(vma0, 0b10110001) ;  // vma0[2], vma0[3], vma0[0], vma0[1]
  vma0 = _mm_max_epi32(vma0, vma1) ;            // max( 0|2|4|6 , 1|3|5|7 , dont' t care )
  vma1 = _mm_shuffle_epi32(vma0, 0b01001011) ;  // vma0[1], vma0[0], vma0[2], vma0[3]
  vma0 = _mm_max_epi32(vma0, vma1) ;            // max( 0|1|2|3|4|5|6|7 , dont' t care )
  _mm_storeu_si32( &max, vma0) ;                // store max

  vmi0 = _mm256_extracti128_si256(vmi, 0) ;     // fold min to 1 value
  vmi1 = _mm256_extracti128_si256(vmi, 1) ;
  vmi0 = _mm_min_epi32(vmi0, vmi1) ;
  vmi1 = _mm_shuffle_epi32(vmi0, 0b10110001) ;
  vmi0 = _mm_min_epi32(vmi0, vmi1) ;
  vmi1 = _mm_shuffle_epi32(vmi0, 0b01001011) ;
  vmi0 = _mm_min_epi32(vmi0, vmi1) ;
  _mm_storeu_si32( &min, vmi0) ;                // store min
  goto end ;
less_than_8:
#endif
  while(nj-- > 0) {
    for(i = 0 ; i < ni ; i++){    // a multiple of 8 values is left
      q[i] = (z[i] * ovq) + ((z[i] < 0) ? -.5f : .5f) ;
      min = (q[i] < min) ? q[i] : min ;
      max = (q[i] > max) ? q[i] : max ;
    }
    z += lniz ;
    q += lniq ;
  }
end:
  t->t[0] = min ; needed1 = NeedBits(min) ;
  t->t[1] = max ; needed2 = NeedBits(max) ;
  return  (needed1 > needed2) ? needed1 : needed2 ;
}

// inverse quantizer for 32 bit floats producing a stream of unsigned 32 bit integers
// nbits is taken from the packing header
// iz      : quantized stream (32 bit elements)  [IN]
// z       : float array to be quantized         [OUT]
// ni      : number of useful points in rows     [IN]
// lni     : storage size of rows in array z     [IN]
// lniz    : storage size of rows in array iz    [IN]
// nj      : number of rows in z and iz          [IN]
// p       : packing header as initialized by float_quantize_prep     [IN]
void float_unquantize(void *iz, float *z, int ni, int lni, int lniz, int nj, QuantizeHeader *p) {
  int i, j, t;
  unsigned int   *izw = (unsigned int   *) iz;
  int offset = p->o;
  int exp = p->e;
  int nbits = p->nbits ;

  if(nbits <= 0 ) return ;
  if(nbits > 24 ) nbits = 24 ;           // no more than 24 bits

  if(ni == lni && ni == lniz) {          // can we fuse loops
    ni = ni * nj ;
    nj = 1 ;
  }
  if(exp > -127 && exp < 127){           // "civilized" exponent for largest value
    FloatInt m1;
    float fac32;
    m1.i = (127 + exp - 23) ;  // inverse of factor to bring largest exponent to 23
    m1.i <<= 23 ;
    fac32 = m1.f;
    for(j=0 ; j<nj ; j++){
      for ( i=0 ; i<ni ; i++){
        t = izw[i] << (24-nbits) ;
        z[i] = (t + offset) * fac32;
      }
      z += lni ;
      izw += lniz ;
    }
  }else{
    DoubleLong m1;
    double fac64;
    exp = (exp > 127) ? 127 : exp ;
    m1.l = (exp + 1023 - 23) ;  // inverse of factor to bring largest exponent to 23
    m1.l = m1.l << 52 ;
    fac64 = m1.d;
    for(j=0 ; j<nj ; j++){
      for ( i=0 ; i<ni ; i++){
        t = izw[i] << (24-nbits) ;
        z[i] = (t + offset) * fac64;
      }
      z += lni ;
      izw += lniz ;
    }
  }
}
// restore float values from integer quantized values
// z       : array of restored floats                   [OUT]
// q       : integer array, quantized values            [IN]
// ni      : number of values                           [IN]
// quantum : quantization unit                          [IN]
// t       : lowest, highest restored values            [OUT]
//           quantum will be internally adjusted to the first power of 2 <= quantum
// there is a SIMD version because the min/max in the loop seem to be a problem for some vectorisers
void float_unquantize_simple_1D(float * restrict z, int32_t * restrict q, int ni, float quantum, FloatPair *t){
  float_unquantize_simple(z, q, ni, ni, ni, 1, quantum, t) ;
}

// 2D version of float_unquantize_simple_1D
// lniz  : storage length of z rows      [IN]
// lniq  : storage length of q rows      [IN]
// nj    : number of rows                [IN]
void float_unquantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, FloatPair *t){
  int i, i0 ;
  float min = 1.0E+38f, max = -min ;
  float vmin[8], vmax[8] ;
  int n7 = (ni & 7) ? (ni & 7) : 8 ;
  FloatInt rq ;

  rq.f = quantum ; rq.i &= 0x7F800000 ; quantum = rq.f ;   // adjust quantum to power of 2 <= quantum
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 vz, vf = _mm256_set1_ps(quantum), vmi = _mm256_set1_ps(+1.0E38), vma = _mm256_set1_ps(-1.0E38) ;
  __m256i vq ;

  if(ni < 8) goto less_than_8 ;                     // less than 8 values, use non SIMD code

  while(nj-- > 0) {
    // N.B. the first and second pass will process some values twice (overlap)
    for(i0 = 0 ; i0 < ni-7 ; i0+=n7 , n7=8){          // batches of 8 values
      vq = _mm256_loadu_si256( (__m256i *)(q+i0) ) ;  // fetch q
      vz = _mm256_cvtepi32_ps(vq) ;                   // convert to float
      vz = _mm256_mul_ps(vz, vf) ;                    // * quantum
      vma = _mm256_max_ps(vma, vz) ;                  // max(z, max)
      vmi = _mm256_min_ps(vmi, vz) ;                  // min(z, min)
      _mm256_storeu_ps(z+i0, vz) ;                    // store z
    }
    z += lniz ;
    q += lniq ;
  }
  _mm256_storeu_ps(vmax, vma) ;                     // fold max, min to 1 value each
  _mm256_storeu_ps(vmin, vmi) ;
  for(i=1 ; i<8 ; i++){                             // fold into vmin[0], vmax[0]
    vmin[0] = (vmin[i] < vmin[0]) ? vmin[i] : vmin[0] ;
    vmax[0] = (vmax[i] > vmax[0]) ? vmax[i] : vmax[0] ;
  }
  min = vmin[0] ;
  max = vmax[0] ;
  goto end ;

less_than_8 :
#endif
  // non SIMD code, also used by SIMD version if less than 8 values
  while(nj-- > 0) {
    for(i=0 ; i<ni ; i++){
      z[i] = quantum * q[i] ;                 // unquantize
      min = (z[i] < min) ? z[i] : min ;       // min
      max = (z[i] > max) ? z[i] : max ;       // max
    }
    z += lniz ;
    q += lniq ;
  }
end:
  t->t[0] = min ;      // build return value
  t->t[1] = max ;
}

