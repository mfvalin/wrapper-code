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

static int32_t e[] = {0, 1, 3, 7, 15, 31, 63, 127, 255} ;

int32_t ieee_quantize(float *f,        // array to quantize (INPUT)
                      float *fmaxa,    // largest absolute value inarray (INPUT)
                      int32_t *q,      // quantized data (OUTPUT)
                      int n,           // number of data elements
                      int nexp,        // number of bits for the exponent part of quantized data
                      int nbits)       // number of bits in quantized data
{
  int i, ex, e0, nm, sign ;
  intflt z0, z, x1, t0, y ;
  int32_t round ;

  if(nexp < 1 || nexp > 8) return 0 ;        // nexp too small or too large
  nm = nbits - nexp ;                        // number of effective mantissa bits
  if(nm < 1 || nm >23) return 0 ;

  z0.f = *fmaxa ;
  e0 = (z0.i >> 23) - 127 ;
  t0.i = e[nexp] << 23 ;
  round = (1 << (23 + nexp - nbits -1)) ;
  for(i = 0 ; i < n ; i++) {
    z.f = f[i] ;
    sign = z.i >> 31 ;
    z.i += round ;
    z.i &= 0x7FFFFFFF ;                                   // get rid of sign
    ex = (z.i >> 23) - 127 ;                              // exponent
    x1.i = ((127 + ex - e0) << 23) | (z.i & 0x7FFFFF) ;   // keep mantissa, alter exponent
    y.f = x1.f * t0.f ;
    q[i] = y.i >> (23 - nm) ;
    q[i] = sign ? -q[i] : q[i] ;
  }
  return e0 ;
}

int32_t ieee_unquantize(float *f,      // restored array (OUTPUT)
                        int32_t *q,    // quantized array (INPUT)
                        int n,         // number of data elements
                        int nexp,      // number of bits for the exponent part of quantized data
                        int e0,        // reference exponent (from ieee_quantize)
                        int nbits)     // number of bits in quantized data
{
  int i ;
  intflt t1, q1 ;
  int nm, sign ;

  if(nexp < 1 || nexp > 8) return 0 ;        // nexp too small or too large

  nm = nbits - nexp ;                        // number of effective mantissa bits
  if(nm < 1 || nm >23) return 0 ;

  t1.i = ((254 - e[nexp]) + e0) << 23 ;
  for(i = 0 ; i < n ; i++) {
    sign = (q[i] < 0) ? 1 : 0 ;
    q1.i = sign ? -q[i] : q[i] ;
    q1.i = q1.i << (23 - nm) ;
    f[i] = (q1.i == 0) ? 0 : q1.f * t1.f ;
    f[i] = sign ? -f[i] : f[i] ;
  }
  return 1 ;
}

#if defined(SELF_TEST)
#define NPTS 35
int main(){
  intflt x1, x2, x3, y, z0, z, t0, t1, t2, fi0, fo0 ;
  uint32_t e, e0, m ;
  int i ;
  float r ;
  int denorm , first_denorm = 1;
  float zmax = 35432.198 ;
  float fi[NPTS], fo[NPTS] ;
  int32_t q[NPTS] ;

  z0.f = zmax ;
//   z.i &= 0xFFFFFF80 ;
  e0 = (z0.i >> 23) - 127 ; 
  fprintf(stdout,"z = %8.8x %f\n", z0.i, z0.f) ;

  fi[0] = zmax ;
  for(i=1 ; i<NPTS ; i++) fi[i] = fi[i-1] * .33 ;
  fi[3] = -fi[3] ;
  e0 = ieee_quantize( fi, &zmax,  q, NPTS, 4, 16) ;
  ieee_unquantize( fo, q, NPTS, 4, e0, 16) ;
  for(i = 0 ; i < 35 ; i++) {
    fi0.f = fi[i] ; fo0.f = fo[i] ; r = fo[i]/fi[i] ;
    denorm = (q[i] >> 11) == 0 ;
    if(denorm && first_denorm) {
      fprintf(stdout,"\n") ;
      first_denorm = 0 ;
    }
    fprintf(stdout,"z0 = %8.8x %12g, coded = %8.8x, z1 = %8.8x %12g, r = %12f %12.2f, d = %9f\n",
                   fi0.i, fi[i],     q[i],          fo0.i, fo[i],    r, fi[i] / (fi[i] - fo[i]), fo[i] - fi[i]);
    if(q[i] == 0) break ;
  }
fprintf(stdout,"====================================\n");
  t0.i = 15 << 23 ;
  t1.i = (254 - 15) << 23 ;
first_denorm = 1;
  for(i = 0 ; i < 35 ; i++) {
    e0 = ieee_quantize( &z0.f, &zmax,  &y.i, 1, 4, 16) ;
    denorm = (y.i >> 11) == 0 ;
    if(denorm && first_denorm) {
      fprintf(stdout,"\n") ;
      first_denorm = 0 ;
    }
    ieee_unquantize( &x3.f, &y.i, 1, 4, e0, 16) ;
    r = x3.f / z0.f ;
    fprintf(stdout,"z0 = %8.8x %12g, x1 = %8.8x %12g, coded = %8.8x, z1 = %8.8x %12g, r = %f\n", 
                    z0.i, z0.f,      x1.i, x1.f,      y.i,       x3.i, x3.f,      r) ;
    if(y.i == 0) break ;
    z0.f *= .33 ;
  }

  return 0 ;
}
#endif
