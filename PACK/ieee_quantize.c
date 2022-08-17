/* Hopefully useful routines for C and FORTRAN
 * Copyright (C) 2021  Recherche en Prevision Numerique
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
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
#include <limits.h>
#include <float.h>

#include <misc_operators.h>
#include <misc_types.h>

#include <ieee_quantize.h>

// largest exponent as a function of exponent bit field width
static int32_t e[] = {0, 1, 3, 7, 15, 31, 63, 127, 255} ;

#define VL 8

// vectorised by VL (POWER of 2) version
void quantize_setup(float *z,            // array to be quantized (IEEE 754 32 bit float) (INPUT)
                        int n,           // number of data elements
                        qhead *h)        // quantization control information (OUTPUT)
{
  int i, j, j0 ;
  float mi[VL], ma[VL], mia[VL], za[VL], z0[VL] ;
  FloatUint rng ;
  float maxabs, minabs ;

  for(i=0 ; i<VL && i<n ; i++) z0[i] = z[i] ;
  for(    ; i<VL ; i++) z0[i] = z[0] ;              // if less than VL values, pad with z0[0]

  for(i=0 ; i < VL ; i++){
    mi[i]  = z0[i] ;                                // minimum value
    ma[i]  = z0[i] ;                                // maximum value
    mia[i] = (z0[i] > 0) ? z0[i] : -z0[i] ;         // smallest absolute value (larger than 0)
  }

  j0 = n & (VL-1) ;                               // skip first mod(n,VL) values
  if(j0 == 0) j0 = VL ;                           // n was a multiple of VL
  for(j=j0 ; j<(n-VL+1) ; j+=VL){
    for(i=0; i<VL ; i++){                         // loops not fused to help gcc optimizer
      ma[i] = MAX(ma[i] , z[j+i]);                // maximum signed value
      mi[i] = MIN(mi[i] , z[j+i]);                // minimum signed value
    }
    for(i=0; i<VL ; i++){                         // loops not fused to help gcc optimizer
      za[i]  = z[j+i]  > 0 ? z[j+i] : -z[j+i] ;   // absolute value of z[j+i]
      za[i]  = z[j+i] != 0 ? za[i]  : mia[i] ;    // if 0, replace with current absolute value minimum
      mia[i] = MIN(mia[i], za[i]) ;               // minimum non zero absolute value 
    }
  }
  for(i=1; i<VL ; i++) {               // final folding pass
    ma[0]  = MAX(ma[0] , ma[i]);
    mi[0]  = MIN(mi[0] , mi[i]);
    mia[0] = MIN(mia[0] , mia[i]);
  }
  h->fmax = ma[0];                             // highest value
  h->fmin = mi[0];                             // lowest value
  h->amin = mia[0] ;                           // smallest non zero absolute value
  h->nbits = 0 ;                               // to be set later
  h->sbit = (h->fmax * h->fmin < 0) ? 1 : 0 ;  // a sign bit is needed, there are positive and negative numbers
  h->negative = ((h->fmax * h->fmin >= 0) && (h->fmin < 0)) ? 1 : 0 ;  // all values are negative
  rng.f = h->fmax - h->fmin ;                     // signed range
  rng.i = ((rng.i >> 23) + 1) << 23 ;             // next power of 2 > rng.f
  h->rng = rng.f ;                                // signed range
  maxabs = h->fmax >= 0 ? h->fmax : -h->fmax ;    // |maxval|
  minabs = h->fmin >= 0 ? h->fmin : -h->fmin ;    // |minval|
  maxabs = maxabs > minabs ? maxabs : minabs ;    // max( |maxval| , |minval| )
  h->fmaxa = maxabs ;
  rng.f = maxabs - h->amin ;
  rng.i = ((rng.i >> 23) + 1) << 23 ;             // next power of 2 > rng.f
  h->rnga = rng.f ;                               // range of absolute values
}

// only keep nbits in the mantissa of IEEE 754 floating point numbers
void ieee_clip(void *f, int n, int nbits){
  int i, j ;
  uint32_t *fi = (uint32_t *) f ;
  uint32_t mask = 0xFFFFFFFF ;

  mask >>= (32 -(23-nbits))  ; // lower 23 -nbits bits
  mask = ~mask ;
  for(i=0 ; i<4 ; i++) fi[i] = fi[i] & mask ;
  for(j=(n&3) ; j<n ; j+=4){
    for(i=0; i<4 ; i++) fi[i+j] = fi[i+j] & mask ;
  }
}

// transform a float into an integer (ieee style)
// f     : float to be quantized
// e0    : true IEEE exponent of largest absolute value (no 127 bias)
// round : mantissa rounding (integer, single bit in appropriate potition)
// t0    : float scaling factor
// nm    : number of effective mantissa bits ( 1 - 23 )
// limit : maximum permitted absolute value once quantized (normally right mask of nbits bits)
// sbit  : (0/1) mask for sign bit. if 0, sign will be ignored
//
// return signed quantized (ieee style) integer value for f
// the absolute value of the input float is converted to a new floating format
// with a reduced size exponent and a reduced size mantissa (nm bits)
// this new "reduced size float" is treated as an integer when restoring sign if needed
static inline int32_t ieee_f_to_q(float f, int32_t e0, int32_t round, float t0, int nm, int32_t limit, int32_t sbit){
  FloatUint z, x1, y ;
  int sign, ex ;
  int q ;                       // final result
  z.f = f ;                     // float will be mostly manipulated as an unsigned integer
  sign = z.i >> 31 ;            // get sign bit (most significant bit)
  sign = sign & sbit ;          // possibly ignore sign bit
  z.i &= 0x7FFFFFFF ;           // suppress FP sign bit
  z.i += round ;                // apply mantissa rounding (may cause an exponent increase by 1)
  ex = (z.i >> 23) ;            // get exponent (including IEEE bias of 127)
  x1.i = ((ex - e0) << 23) |    // alter exponent (largest value becomes 127), keep mantissa intact
         (z.i & 0x7FFFFF) ;     // (equivalent to dividing by 2**e0)
  y.f = x1.f * t0 ;             // apply scaling (may produce a denormalized float)
  q = y.i >> (23 - nm) ;        // get rid of unused rightmost mantissa bits
  q = q > limit ? limit : q ;   // limit is expected to have the nbits rightmost bits set
  q = sign ? -q : q ;           // restore sign by negating integer quantized value if necessary
  return q ;
}

// restore the float value from the transformed value
// q    : quantized value representing a float
// t1   : first scaling factor
// t2   : second scaling factor
// nm   : number of effective mantissa bits ( 1 - 23 )
// sbit : (0/1) mask for sign bit. if 0, sign will be ignored
// neg  : (0/1) restored float value will be negative if 1 (usually sbit is 0 if neg is 1)
//
// return restored float value
// t1 * t2 might generate an overflow, which is why they MUST be applied separately using 2 multiplies
static inline float ieee_q_to_f_2(int32_t q, float t1, float t2, int nm, int32_t sbit, int32_t neg){
  int sign ;
  FloatUint q1 ;
  float f ;                     // final result
  sign = (q < 0) ? 1 : 0 ;      // get sign
  sign = sign & sbit ;          // possibly unsigned quantized value
  q1.i = sign ? -q : q ;        // get absolute value if signed
  q1.i = q1.i << (23 - nm) ;    // shift left into proper place
  f = (q1.f * t1) * t2 ;        // apply the 2 scaling factors
  sign = sign | neg ;
  f = sign ? -f : f ;           // restore sign if negative
  return f ;
}
// single factor version if t1 * t2 are known not to create an overflow
static inline float ieee_q_to_f_1(int32_t q, float t1t2, int nm, int32_t sbit, int32_t neg){
  int sign ;
  FloatUint q1 ;
  float f ;                     // final result
  sign = (q < 0) ? 1 : 0 ;      // get sign
  sign = sign & sbit ;          // possibly unsigned quantized value
  q1.i = sign ? -q : q ;        // get absolute value if signed
  q1.i = q1.i << (23 - nm) ;    // shift left into proper place
  f = q1.f * t1t2 ;             // apply the combined scaling factor
  sign = sign | neg ;
  f = sign ? -f : f ;           // restore sign if negative
  return f ;
}

// transform IEEE 754 floating point numbers into signed integers
// a positive float value will be transformed into a nbits integer number
// the upper nexp bits contain a reduced binary exponent
// the lower (nbits -nexp) bits contain a mantissa
// the "hidden one" and "denormalization" features of IEEE 754 are used
// example : nbits = 16, nexp = 4
//           the largest value  0 eeeeeeee mmmmmmmmmmmmxxxxxxxxxxx becomes
//                                    1111 mmmmmmmmmmmm
//           denormalized value 0 00000000 ddddddddddddxxxxxxxxxxx becomes
//                                    0000 dddddddddddd
// numbers 0.0 -> *fmaxa get converted to range
//         0.0 -> 2 ** (e[nexp] + 1 -127)
//         taking advantage of denormalized numbers to extend the dynamic range
int32_t ieee_quantize(float *f,        // array to quantize (IEEE 754 32 bit float) (INPUT)
                      int32_t *q,      // quantized data (OUTPUT)
                      int n,           // number of data elements (INPUT)
                      int nexp,        // number of bits for the exponent part of quantized data (INPUT)
                      int nbits,       // number of bits in quantized data (INPUT)
                      qhead *h)        // quantization control information (INPUT+OUTPUT)
{
  int i, e0, nm ;
  FloatUint z0, t0 ;
  int32_t round, min, max ;
  int32_t limit, sbit ;
  float fmaxa ;    // largest absolute value in array

  e0 = -128 ;                                // invalid true exponent
  if(h == NULL) return e0 ;
  fmaxa = h->fmaxa ;
  if(nexp < 1 || nexp > 8) return e0 ;       // nexp too small or too large
  sbit = h->sbit ;
  nbits = nbits - sbit ;                     // need for a sign bit ? (if so reduce nbits)
  limit = ~((-1) << nbits) ;
  nm = nbits - nexp ;                        // number of effective mantissa bits
  if(nm < 1 || nm >23) return e0 ;           // too few or too many mantissa bits

  z0.f = fmaxa ;
  e0 = (z0.i >> 23) - 127 ;                  // true exponent of largest absolute value
  t0.i = e[nexp] << 23 ;                     // final scaling factor
  round = (1 << (23 + nexp - nbits -1)) ;    // rounding term

  min = INT_MAX ; max = INT_MIN ;
  for(i = 0 ; i < n ; i++) {
    q[i] = ieee_f_to_q(f[i], e0, round, t0.f, nm, limit, sbit) ;
    min = (q[i] < min) ? q[i] : min ;
    max = (q[i] > max) ? q[i] : max ;
  }
  if(h != NULL){
    h->e0 = e0 ;             // true exponent of largest absolute value float
    h->nbits = nbits ;       // number of bits per token
    h->nexp = nexp ;         // number of exponent bits
//     h->min = min ;           // lowest quantized signed value
//     h->max = max ;           // largest quantized signed value
    h->max = ieee_f_to_q(h->fmax, e0, round, t0.f, nm, limit, sbit) ;
    h->min = ieee_f_to_q(h->fmin, e0, round, t0.f, nm, limit, sbit) ;
    h->limit = limit ;       // keep limit mask
  }
  return e0 ;
}

// vector version of above
int32_t ieee_quantize_v4(float *f,        // array to quantize (IEEE 754 32 bit float) (INPUT)
                      int32_t *q,      // quantized data (OUTPUT)
                      int n,           // number of data elements (INPUT)
                      int nexp,        // number of bits for the exponent part of quantized data (INPUT)
                      int nbits,       // number of bits in quantized data (INPUT)
                      qhead *h)        // quantization control information (INPUT+OUTPUT)
{
  int i, j, e0, nm ;
  FloatUint z0, t0 ;
  int32_t round, min, max ;
  int32_t limit, sbit ;
  int32_t vl ;
  float fmaxa ;    // largest absolute value in array

  e0 = -128 ;                                // invalid true exponent
  if(h == NULL) return e0 ;
  fmaxa = h->fmaxa ;
  if(nexp < 1 || nexp > 8) return e0 ;       // nexp too small or too large
  sbit = h->sbit ;
  nbits = nbits - sbit ;                     // need for a sign bit ?
  limit = ~((-1) << nbits) ;
  nm = nbits - nexp ;                        // number of effective mantissa bits
  if(nm < 1 || nm >23) return e0 ;           // too few or too many mantissa bits

  z0.f = fmaxa ;
  e0 = (z0.i >> 23) - 127 ;                  // true exponent of largest absolute value
  t0.i = e[nexp] << 23 ;                     // final scaling factor
  round = (1 << (23 + nexp - nbits -1)) ;

  vl = (n & 3) ; vl = (vl == 0) ? 4 : vl ;
  for(j = 0 ; j < n-3 ;){          // n is ASSUMED TO BE > 3
    for(i = 0 ; i < 4 ; i++){
      q[j+i] = ieee_f_to_q(f[j+i], e0, round, t0.f, nm, limit, sbit) ;
    }
    j += vl ;
    vl = 4 ;
  }
  if(h != NULL){
    h->e0 = e0 ;
    h->nbits = nbits ;
    h->nexp = nexp ;
    h->max = ieee_f_to_q(h->fmax, e0, round, t0.f, nm, limit, sbit) ;
    h->min = ieee_f_to_q(h->fmin, e0, round, t0.f, nm, limit, sbit) ;
    h->limit = limit ;       // keep limit mask
  }
  return e0 ;
}

// restore float values from quantized (ieee style) values
int32_t ieee_unquantize(float *f,      // restored array (IEEE 754 32 bit float) (OUTPUT)
                        int32_t *q,    // quantized array (INPUT)
                        int n,         // number of data elements (INPUT)
                        qhead *h)      // quantization control information (INPUT)
{
  int i ;
  FloatUint t1, t2 ;
  int nm ;
  int nexp, e0, nbits ;
  int32_t sbit, neg ;

  if(h == NULL) return 0 ;
  nbits = h->nbits ;                         // number of bits in quantized data
  nexp = h->nexp ;                           // number of bits for the exponent part of quantized data
  if(nexp < 1 || nexp > 8) return 0 ;        // nexp too small or too large
  e0 = h->e0 ;                               // reference exponent (from ieee_quantize)
  if(e0 > 127 || e0 < -127 ) return 0 ;      // invalid reference exponent

  nm = nbits - nexp ;                        // number of effective mantissa bits
  if(nm < 1 || nm >23) return 0 ;

  sbit = h->sbit ;
  neg  = h->negative ;
  if(e0 > e[nexp]) {                         // must use 2 factors if e0 > e[nexp]
fprintf(stdout,"BEEP\n");
    t1.i = ((254 - e[nexp]) << 23) ;
    t2.i = ((127 + e0)      << 23) ;         // t1.f * t2.f would be too large ( > 2**128 )
    for(i = 0 ; i < n ; i++) {
      f[i] = ieee_q_to_f_2(q[i], t1.f, t2.f, nm, sbit, neg) ;
    }
  }else{                                     // can use 1 factor if e0 <= e[nexp]
fprintf(stdout,"BOP\n");
    t1.i = ((254 - e[nexp] + e0) << 23) ;
    for(i = 0 ; i < n ; i++) {
      f[i] = ieee_q_to_f_1(q[i], t1.f, nm, sbit, neg) ;
    }
  }
  return 1 ;
}

// IEEE 32 bit floating point to half precision (16 bit) IEEE floating point
// any number >= 65520 will be coded as infinity in FP16
static inline uint16_t ieee_fp32_to_fp16(float f){
  FloatUint z, y ;
  uint32_t sign ;
  uint32_t round = 0x1000 ;
  uint32_t limit = ((127+16) << 23) | 0x7FFFFF ; // largest representable FP16
  z.f = f ;                     // float will be mostly manipulated as an integer
  sign = (z.i >> 16) & 0x8000 ; // position of FP16 sign bit
  z.i &= 0x7FFFFFFF ;           // suppress FP sign bit
  z.i += round ;                // apply mantissa rounding
  z.i = (z.i > limit) ? limit : z.i ;
  y.i = 15 << 23 ;              // scale by 2 ** (15 -127)
  y.f *= z.f ;
  y.i = (y.i >> 13) & 0xFFFF ;  // suppress lower 16 bits of mantissa & reduce exp + mantissa to 15 bits
  y.i |= sign ;                 // apply sign
  return y.i ;
}

// IEEE 32 bit floating point to half precision (16 bit) IEEE floating point
void fp32_to_fp16(float *f, uint16_t *q, int n){
  int i ;
  for(i = 0 ; i < n ; i++) q[i] = ieee_fp32_to_fp16(f[i]) ;
}

// scaled IEEE 32 bit floating point to half precision (16 bit) IEEE floating point
void fp32_to_fp16_scaled(float *f, uint16_t *q, int n, float scale){
  int i ;
  for(i = 0 ; i < n ; i++) q[i] = ieee_fp32_to_fp16(f[i]*scale) ;
}

// IEEE 32 bit floating point to brain float (16 bit)
// current code will not behave correctly if upper 16 bits after sign bit in float are 1 (NaN)
// removing comments makes code safe
void fp32_to_bf16(float *f, uint16_t *q, int n){
  FloatUint z ;
  uint32_t round = 1 << 15 ;
  uint32_t sign ;
  int i ;
  for(i = 0 ; i < n ; i++) {
    z.f = f[i] ;
//     sign = z.i & (~0x7FFFFFFF) ;
//     z.i &= 0x7FFFFFFF ;
    z.i += round ;
//     z.i &= 0x7FFFFFFF ;
//     z.i |= sign ;
    q[i] = z.i >> 16 ;
  }
}

// half precision (16 bit) IEEE floating point to IEEE 32 bit floating point
// the infinity argment value is used as a result if the FP16 exponent is 31
// (the sign of the FP16 value will be preserved)
static inline float ieee_fp16_to_fp32(uint16_t q, uint32_t infinity){
  FloatUint z, y ;
  int sign ;
  int e0;
  sign = q & 0x8000 ;    // position of FP16 sign bit
  sign <<= 16 ;          // position of FP32 sign bit
  z.i = q & 0x7FFF ;     // suppress FP16 sign bit
  e0 = z.i >> 10 ;       // FP16 exponent (with bias 15)
  z.i <<= 13 ;
  y.i = (254 - 15) << 23 ;       // scale by 2 ** (127 - 15)
  z.f *= y.f ;
  z.i = (e0 == 31) ? infinity : z.i ;  // infinity if biased FP16 exponent == 31
  z.i |= sign ;
  return z.f ;
}

// half precision (16 bit) IEEE floating point to IEEE 32 bit floating point
// if inf is not NULL, *inf is used instead of IEEE Infinity (sign of FP value is preserved)
// 0X7F800000 is infinity, 0X7F800001 -> 0X7FFFFFFF is a NaN (not a number)
// f  : 32 bit IEEE floating point array
// q  : 16 bit IEEE floating point array
void fp16_to_fp32(float *f, void *f16, int n, void *inf){
  int i ;
  uint32_t Inf = 0X7F800000 ; //  IEEE Infinity
  uint16_t *q = f16 ;

  Inf = inf ? *(uint32_t *)inf : Inf ;
  for(i = 0 ; i < n ; i++) f[i] = ieee_fp16_to_fp32(q[i], Inf) ;
}

#if defined(__clang__) || defined(__ICC) || defined(__PGIC__) || defined(VANILLA)
#else
// NOTE: strange method to defeat some overzealous optimizers at O2 and above
// static uint32_t Fetch_32(uint32_t *what) {  // return dereferenced what or IEEE Infinity if what is NULL
//   return what ? *what : 0X7F800000 ; //  IEEE Infinity
// }
// uint32_t (*FudgedWordFetch)(uint32_t *what) = &Fetch_32 ;
// NOTE: strange method to defeat some overzealous optimizers at O2 and above
static float Reciprocal(float what){  // return reciprocal of a float
  return 1.0f / what ;
}
float (*FudgeFloatReciprocal)(float what) = &Reciprocal ;
#endif

// scaled half precision (16 bit) IEEE floating point to IEEE 32 bit floating point
// if inf is not NULL, *inf is used instead of IEEE Infinity (sign of FP value is preserved)
// 0X7F800000 is infinity, 0X7F800001 -> 0X7FFFFFFF is a NaN (not a number)
void fp16_to_fp32_scaled(float *f, void *f16, int n, void *inf, float scale){
  int i ;
  uint32_t Inf ;
  uint16_t *q = f16 ;
  float rscale ;

#if defined(__clang__) || defined(__ICC) || defined(__PGIC__) || defined(VANILLA)
  rscale = 1.0f / scale ;
  Inf = inf ? *((uint32_t *)inf) : 0X7F800000 ;    // IEEE infinity
#else
  rscale = (*FudgeFloatReciprocal)(scale) ;  // NOTE: defeat some overzealous at O2 and above
//   Inf = (*FudgedWordFetch)(inf) ;  // NOTE: defeat some overzealous at O2 and above
  Inf = inf ? *((uint32_t *)inf) : 0X7F800000 ;    // IEEE infinity
#endif
  for(i = 0 ; i < n ; i++) f[i] = rscale * ieee_fp16_to_fp32(q[i], Inf) ;
}

// brain float (16 bit) to IEEE 32 bit floating point
void bf16_to_fp32(float *f, uint16_t *q, int n){
  FloatUint z ;
  int i ;
  for(i = 0 ; i < n ; i++){
    z.i = q[i] ;
    z.i <<= 16 ;   // shift to upper 16 bits
    q[i] = z.f ;
  }
}

