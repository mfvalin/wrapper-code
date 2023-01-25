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

#include <rmn/misc_operators.h>
#include <rmn/misc_types.h>

#include <ieee_quantize.h>
#include <ieee_quantize.h>

#define NPT  8
#define NPTS 38
#define N    35
#define NEXP 4

int main(){
  FloatUint x1, x2, x3, y, z0, z, t0, t1, t2, fi0, fo0 ;
  uint32_t e, e0, m ;
  int i, j, nbits, nbits0 ;
  float r ;
  int denorm , first_denorm = 1;
  float    Big = FLT_MAX ;
//   float    Big = 66000.1 ;
  uint32_t Inf = 0X7F800000 ;
  uint32_t NaN = 0X7F800001 ;
  float zmax = 65432.898 ;
  float scale = 1.0f ;
//   float zmax = 8190.99 ;
//   float zmax = 7.0E+4 ;
  float fi[NPTS], fo[NPTS] ;
  int32_t q[NPTS] ;
  qhead h ;
  float fz0[NPT] ;
//   float fp32 ;
  qhead he ;
//   uint16_t fp16 ;
  uint16_t vfp16[NPTS] ;
  uint32_t limit16 = ((127+14) << 23) | 0x7FFFFF ; // largest representable FP16

  fprintf(stdout, "limit16 = %8.8x, %8d, %8.8x\n", limit16, limit16 >> 23, limit16 & 0x7FFFFF);
  fi[0] = 1.0 ;
  for(i=1 ; i < NPTS ; i++) { fi[i] = -fi[i-1] * 2.0f ; }
//   fi[17] = 0.000060975552 ;
  fi[16] = 65519.9980468749999999f ;  // largest value that does not decode to infinity
  fi[16] = 65519.998046875f ;      // generates infinity after decoding
  fi[18] = 0.00006104 ;
  fi[19] = 65536 - 56 ;
  for(i = 20 ; i < NPTS ; i++) { fi[i] = fi[i-1] + 1.0f ; }
  fp32_to_fp16_scaled(fi, vfp16, NPTS, scale) ;
  fp16_to_fp32_scaled(fo, vfp16, NPTS, NULL, scale) ;
//   fp16_to_fp32_scaled(fo, vfp16, NPTS, (void *) &Big, scale) ;
//   fp16_to_fp32(fo, vfp16, NPTS, (void *) &Big) ;
//   fp16_to_fp32_scaled(fo, vfp16, NPTS, &NaN, scale) ;
//   fp16_to_fp32_scaled(fo, vfp16, NPTS, &Inf, scale) ;
  for(i=0 ; i < NPTS ; i++){
    x1.f = fi[i] ;
//     fp16 = vfp16[i] ;
//     fp32 = fo[i] ;
    fprintf(stdout, "fp32 = %12g (%12g) (%8.8x), fp16 = %8.8x (%2d,%4.4x)\n", 
            fi[i], fo[i], x1.u, vfp16[i], vfp16[i]>>10, vfp16[i] & 0x3FF) ;
  }
// return 0 ;
//   for(i=0 ; i<NPT ; i++) fz0[i] = i - (NPT-1)/2.0f ;
//   fprintf(stdout,"fz0 :");
//   for(i=0 ; i<NPT ; i++) fprintf(stdout," %f",fz0[i]) ; fprintf(stdout,"\n\n");
//   for(j=1 ; j<=NPT ; j++){
//     quantize_setup(fz0, j, &he);
//     fprintf(stdout,"[%2d] max = %10f, min = %10f, mina = %10f, rng = %10f, rnga = %10f\n", 
//             j, he.fmax, he.fmin, he.amin, he.rng, he.rnga);
//   }
// return 0 ;
  ieee_clip(fi, 0, 16) ;
  fi[0] = zmax ;
  for(i=1 ; i<NPTS ; i++) fi[i] = fi[i-1] * .499 ;
//   fi[5] = -fi[5] ;
//   fi[3] = -fi[0] ;
  for(i=0 ; i<NPTS ; i++) fi[i] = -fi[i] ;
  for(i=1 ; i<NPTS ; i++) q[i] = 0 ;
  quantize_setup(fi, N, &h);
  fprintf(stdout,"fi  : max = %f, min = %f, mina = %f, rng = %f, rnga = %f\n", 
          h.fmax, h.fmin, h.amin, h.rng, h.rnga);
  nbits0 = 16 ;
  nbits = nbits0 ;
//   if(h.fmin * h.fmax < 0) nbits-- ;  // positive and negative numbers, need to reserve a bit for the sign
//   e0 = ieee_quantize( fi, &zmax,  q, N, NEXP, nbits, &h) ;
  e0 = ieee_quantize_v4( fi, q, N, NEXP, nbits, &h) ;
  fprintf(stdout,"nexp = %d, nbits = %d, e0 = %d %d, min = %d, max = %d, span = %d, limit = %8.8x, sbit = %d, neg = %d\n",
          h.nexp, h.nbits, h.e0, e0, h.min, h.max, h.max-h.min, h.limit, h.sbit, h.negative) ;
//   ieee_unquantize( fo, q, N, NEXP, e0, 16, &h) ;
  ieee_unquantize( fo, q, N, &h) ;
  for(i = 0 ; i < N ; i++) {
    fi0.f = fi[i] ; fo0.f = fo[i] ;
    r = fo[i]/fi[i] ; r = r - 1.0 ; r = (r>0) ? r : -r ;
    denorm = (ABS(q[i]) >> (nbits0-NEXP)) == 0 ;
    if(denorm && first_denorm) {
      fprintf(stdout,"\n") ;
      first_denorm = 0 ;
    }
    fprintf(stdout,"z0 = %8.8x %12g, coded = %8.8x, z1 = %8.8x %12g, r = %12f %10.0f %12.0f, d = %9f\n",
                   fi0.u, fi[i],     q[i],          fo0.u, fo[i],    r, 1.0/r, zmax / (fi[i] - fo[i]), fo[i] - fi[i]);
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
