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
#include <misc_types.h>
#include <misc_operators.h>
#include <misc_pack.h>
#if ! defined(__INTEL_COMPILER_UPDATE)
#pragma GCC optimize "tree-vectorize"
#endif

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
    m1.f = m2.f = quantum ;
    m1.i &= (~0x7FFFFF) ;                    // first power of 2 <= quantum
    irange = (range / m1.f) ;
    NEEDBITS( irange, neededbits) ;
// printf("new nbits = %d, quantum = %g (%8.8x) -> %g (%8.8x), range = %g\n", neededbits, quantum, m2.i, m1.f, m1.i, range) ;
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

