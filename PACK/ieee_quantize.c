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
typedef union{
  uint32_t i ;
  float    f ;
} intflt ;

typedef struct{
  int32_t e0 ;
  int32_t nbits ;
  int32_t nexp ;
  int32_t min ;
  int32_t max ;
} qhead ;

static int32_t e[] = {0, 1, 3, 7, 15, 31, 63, 127, 255} ;

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
// numbers 0.0 -> *fmaxa get converted tot range
//         0.0 -> 2 ** (e[nexp] + 1 -127)
//         taking advantage of denormalized numbers to extend the dynamic range
int32_t ieee_quantize(float *f,        // array to quantize (IEEE 754 32 bit float) (INPUT)
                      float *fmaxa,    // largest absolute value inarray (INPUT)
                      int32_t *q,      // quantized data (OUTPUT)
                      int n,           // number of data elements
                      int nexp,        // number of bits for the exponent part of quantized data (INPUT)
                      int nbits,       // number of bits in quantized data (INPUT)
                      qhead *h)        // quantization control information (OUTPUT)
{
  int i, ex, e0, nm, sign ;
  intflt z0, z, x1, t0, y ;
  int32_t round, min, max ;

  e0 = -128 ;                                // invalid true exponent
  if(h == NULL) return e0 ;
  if(nexp < 1 || nexp > 8) return e0 ;       // nexp too small or too large
  nm = nbits - nexp ;                        // number of effective mantissa bits
  if(nm < 1 || nm >23) return e0 ;           // too few or too many mantissa bits
  z0.f = *fmaxa ;
  e0 = (z0.i >> 23) - 127 ;                  // true exponent of largest absolute value
  t0.i = e[nexp] << 23 ;                     // final scaling factor
  round = (1 << (23 + nexp - nbits -1)) ;
  min = 0x7FFFFFFF ; max = -min ;
  for(i = 0 ; i < n ; i++) {
    z.f = f[i] ;                                          // floats manipulated as integers
    sign = z.i >> 31 ;                                    // get sign bit
    z.i &= 0x7FFFFFFF ;                                   // get rid of sign
    z.i += round ;                                        // apply rounding
    ex = (z.i >> 23) - 127 ;                              // exponent
    x1.i = ((127 + ex - e0) << 23) | (z.i & 0x7FFFFF) ;   // keep mantissa, alter exponent (largest value gets 127)
    y.f = x1.f * t0.f ;                                   // apply final scale factor, result may be denormalized
    q[i] = y.i >> (23 - nm) ;                             // shift to reduce to nbits
    q[i] = sign ? -q[i] : q[i] ;                          // restore sign
    min = (q[i] < min) ? q[i] : min ;
    max = (q[i] > max) ? q[i] : max ;
  }
  if(h != NULL){
    h->e0 = e0 ;
    h->nbits = nbits ;
    h->nexp = nexp ;
    h->min = min ;
    h->max = max ;
  }
  return e0 ;
}

int32_t ieee_unquantize(float *f,      // restored array (IEEE 754 32 bit float) (OUTPUT)
                        int32_t *q,    // quantized array (INPUT)
                        int n,         // number of data elements (INPUT)
                        qhead *h)      // quantization control information (INPUT)
{
  int i ;
  intflt t1, t2, q1 ;
  int nm, sign ;
  int nexp, e0, nbits ;

  if(h == NULL) return 0 ;
  nbits = h->nbits ;                         // number of bits in quantized data
  nexp = h->nexp ;                           // number of bits for the exponent part of quantized data
  if(nexp < 1 || nexp > 8) return 0 ;        // nexp too small or too large
  e0 = h->e0 ;                               // reference exponent (from ieee_quantize)
  if(e0 > 127 || e0 < -127 ) return 0 ;      // invalid reference exponent

  nm = nbits - nexp ;                        // number of effective mantissa bits
  if(nm < 1 || nm >23) return 0 ;

  if(e0 > e[nexp]) {                         // must use 2 factors if e0 > e[nexp]
    t1.i = ((254 - e[nexp]) << 23) ;
    t2.i = ((127 + e0)      << 23) ;         // t1.f * t2.f would be too large ( > 2**128 )
    for(i = 0 ; i < n ; i++) {
      sign = (q[i] < 0) ? 1 : 0 ;            // extract sign
      q1.i = sign ? -q[i] : q[i] ;           // positive value
      q1.i = q1.i << (23 - nm) ;
      f[i] = q1.f * t1.f * t2.f ;
//       f[i] = (q1.i == 0) ? 0 : q1.f * t1.f * t2.f ;
      f[i] = sign ? -f[i] : f[i] ;           // restore sign
    }
  }else{                                     // can use 1 factor if e0 <= e[nexp]
    t1.i = ((254 - e[nexp] + e0) << 23) ;
    for(i = 0 ; i < n ; i++) {
      sign = (q[i] < 0) ? 1 : 0 ;            // extract sign
      q1.i = sign ? -q[i] : q[i] ;           // positive value
      q1.i = q1.i << (23 - nm) ;
      f[i] = q1.f * t1.f ;
//       f[i] = (q1.i == 0) ? 0 : q1.f * t1.f ;
      f[i] = sign ? -f[i] : f[i] ;           // restore sign
    }
  }
  return 1 ;
}

#if defined(SELF_TEST)

#define NPTS 35
#define N    35
#define NEXP 4

int main(){
  intflt x1, x2, x3, y, z0, z, t0, t1, t2, fi0, fo0 ;
  uint32_t e, e0, m ;
  int i ;
  float r ;
  int denorm , first_denorm = 1;
  float zmax = 65432.898 ;
  float fi[NPTS], fo[NPTS] ;
  int32_t q[NPTS] ;
  qhead h ;

//   z0.f = zmax ;
//   z.i &= 0xFFFFFF80 ;
//   e0 = (z0.i >> 23) - 127 ; 
//   fprintf(stdout,"z = %8.8x %f\n", z0.i, z0.f) ;

  fi[0] = zmax ;
  for(i=1 ; i<NPTS ; i++) fi[i] = fi[i-1] * .499 ;
  fi[3] = -fi[3] ;
  for(i=1 ; i<NPTS ; i++) q[i] = 0 ;
  e0 = ieee_quantize( fi, &zmax,  q, N, NEXP, 16, &h) ;
  fprintf(stdout,"nexp = %d, nbits = %d, e0 = %d %d, min = %d, max = %d, span = %d\n",h.nexp, h.nbits, h.e0, e0, h.min, h.max, h.max-h.min) ;
//   ieee_unquantize( fo, q, N, NEXP, e0, 16, &h) ;
  ieee_unquantize( fo, q, N, &h) ;
  for(i = 0 ; i < N ; i++) {
    fi0.f = fi[i] ; fo0.f = fo[i] ;
    r = fo[i]/fi[i] ; r = r - 1.0 ; r = (r>0) ? r : -r ;
    denorm = (q[i] >> 11) == 0 ;
    if(denorm && first_denorm) {
      fprintf(stdout,"\n") ;
      first_denorm = 0 ;
    }
    fprintf(stdout,"z0 = %8.8x %12g, coded = %8.8x, z1 = %8.8x %12g, r = %12f %10.0f %12.0f, d = %9f\n",
                   fi0.i, fi[i],     q[i],          fo0.i, fo[i],    r, 1.0/r, zmax / (fi[i] - fo[i]), fo[i] - fi[i]);
    if(q[i] == 0) break ;
  }
return 0 ;
#if 0
fprintf(stdout,"====================================\n");
  t0.i = 15 << 23 ;
  t1.i = (254 - 15) << 23 ;
first_denorm = 1;
  for(i = 0 ; i < 35 ; i++) {
    e0 = ieee_quantize( &z0.f, &zmax,  &y.i, 1, NEXP, 16, h) ;
    denorm = (y.i >> 11) == 0 ;
    if(denorm && first_denorm) {
      fprintf(stdout,"\n") ;
      first_denorm = 0 ;
    }
    ieee_unquantize( &x3.f, &y.i, 1, NEXP, e0, 16, h) ;
    r = x3.f / z0.f ;
    fprintf(stdout,"z0 = %8.8x %12g, x1 = %8.8x %12g, coded = %8.8x, z1 = %8.8x %12g, r = %f\n", 
                    z0.i, z0.f,      x1.i, x1.f,      y.i,       x3.i, x3.f,      r) ;
    if(y.i == 0) break ;
    z0.f *= .33 ;
  }

  return 0 ;
#endif
}
#endif
