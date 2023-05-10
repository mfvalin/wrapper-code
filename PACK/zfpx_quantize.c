//
// Hopefully useful code
// Copyright (C) 2022  Recherche en Prevision Numerique
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
// 
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// 
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// transforms and helper functions for eventual use in floating point and integer compression
//
// some code has been borrowed from / inspired by / derived from
// https://github.com/LLNL/zfp (Copyright (c) 2014-2022, Lawrence Livermore National Security, LLC)

#include <stdint.h>
#include <stdio.h>

#include <zfpx_quantize.h>

#if ! defined(__INTEL_COMPILER_UPDATE) && ! defined(__PGIC__)
#pragma GCC optimize "tree-vectorize"
#endif
#include <with_simd.h>

#if 0
// borrowed from the source code of zfp
#define ptrdiff_t int
/* reversible forward lifting transform of 4-vector */
static void
rev_fwd_lift_Int(int32_t* p, int s)
{
  int32_t x, y, z, w;
  x = *p; p += s;
  y = *p; p += s;
  z = *p; p += s;
  w = *p; p += s;

  /*                                       (P4 Pascal matrix)
  ** high-order Lorenzo transform       ** high-order inverse Lorenzo transform
  ** ( 1  0  0  0) (x)                  ** ( 1  0  0  0) (x)
  ** (-1  1  0  0) (y)                  ** ( 1  1  0  0) (y)
  ** ( 1 -2  1  0) (z)                  ** ( 1  2  1  0) (z)
  ** (-1  3 -3  1) (w)                  ** ( 1  3  3  1) (w)
  */
  w -= z; z -= y; y -= x;               w += z;
  w -= z; z -= y;                       z += y; w += z;
  w -= z;                               y += x; z += y; w += z;

  p -= s; *p = w;
  p -= s; *p = z;
  p -= s; *p = y;
  p -= s; *p = x;
}
// borrowed from the source code of zfp
/* reversible inverse lifting transform of 4-vector */
static void
rev_inv_lift_Int(int32_t* p, int s)
{
  int32_t x, y, z, w;
  x = *p; p += s;
  y = *p; p += s;
  z = *p; p += s;
  w = *p; p += s;

  /*
  ** high-order Lorenzo transform (P4 Pascal matrix)
  ** ( 1  0  0  0) (x)
  ** ( 1  1  0  0) (y)
  ** ( 1  2  1  0) (z)
  ** ( 1  3  3  1) (w)
  */
  w += z;
  z += y; w += z;
  y += x; z += y; w += z;

  p -= s; *p = w;
  p -= s; *p = z;
  p -= s; *p = y;
  p -= s; *p = x;
}
// borrowed from the source code of zfp
/* forward lifting transform of 4-vector */
static void
fwd_lift_Int(int32_t* p, ptrdiff_t s)
{
  int32_t x, y, z, w;
  x = *p; p += s;
  y = *p; p += s;
  z = *p; p += s;
  w = *p; p += s;
  /*
  ** non-orthogonal transform   ** non-orthogonal inverse transform
  **        ( 4  4  4  4) (x)   **       ( 4  6 -4 -1) (x)
  ** 1/16 * ( 5  1 -1 -5) (y)   ** 1/4 * ( 4  2  4  5) (y)
  **        (-4  4  4 -4) (z)   **       ( 4 -2  4 -5) (z)
  **        (-2  6 -6  2) (w)   **       ( 4 -6 -4  1) (w)
  */
  x += w; x >>= 1; w -= x;
  z += y; z >>= 1; y -= z;
  x += z; x >>= 1; z -= x;
  w += y; w >>= 1; y -= w;
  w += y >> 1; y -= w >> 1;

  p -= s; *p = w;
  p -= s; *p = z;
  p -= s; *p = y;
  p -= s; *p = x;
}
// borrowed from the source code of zfp
/* inverse lifting transform of 4-vector */
static void
inv_lift_Int(int32_t* p, ptrdiff_t s)
{
  int32_t x, y, z, w;
  x = *p; p += s;
  y = *p; p += s;
  z = *p; p += s;
  w = *p; p += s;

  /*
  ** non-orthogonal transform
  **       ( 4  6 -4 -1) (x)
  ** 1/4 * ( 4  2  4  5) (y)
  **       ( 4 -2  4 -5) (z)
  **       ( 4 -6 -4  1) (w)
  */
  y += w >> 1; w -= y >> 1;
  y += w; w <<= 1; w -= y;
  z += x; x <<= 1; x -= z;
  y += z; z <<= 1; z -= y;
  w += x; x <<= 1; x -= w;

  p -= s; *p = w;
  p -= s; *p = z;
  p -= s; *p = y;
  p -= s; *p = x;
}
#endif

// borrowed from the source code of zfp
// reorder table for 2 dimensional transform
#define cache_align_(x) x __attribute__((aligned(64)))
#define index(i, j) ((i) + 4 * (j))
/* order coefficients (i, j) by i + j, then i^2 + j^2 */
// cache_align_(static const uint8_t perm_2[16]) = {
cache_align_(const uint8_t zfpx_perm_2[16]) = {
  index(0, 0), /*  0 : 0 */

  index(1, 0), /*  1 : 1 */
  index(0, 1), /*  2 : 1 */

  index(1, 1), /*  3 : 2 */

  index(2, 0), /*  4 : 2 */
  index(0, 2), /*  5 : 2 */

  index(2, 1), /*  6 : 3 */
  index(1, 2), /*  7 : 3 */

  index(3, 0), /*  8 : 3 */
  index(0, 3), /*  9 : 3 */

  index(2, 2), /* 10 : 4 */

  index(3, 1), /* 11 : 4 */
  index(1, 3), /* 12 : 4 */

  index(3, 2), /* 13 : 5 */
  index(2, 3), /* 14 : 5 */

  index(3, 3), /* 15 : 6 */
};
// reorder after 2 dimensional forward transform, explicit "no gather" version
void zfpx_shuffle_2d(int32_t *src, int32_t *dst){
  dst[ 0] = src[index(0, 0)] ;
  dst[ 1] = src[index(1, 0)] ;
  dst[ 2] = src[index(0, 1)] ;
  dst[ 3] = src[index(1, 1)] ;
  dst[ 4] = src[index(2, 0)] ;
  dst[ 5] = src[index(0, 2)] ;
  dst[ 6] = src[index(2, 1)] ;
  dst[ 7] = src[index(1, 2)] ;
  dst[ 8] = src[index(3, 0)] ;
  dst[ 9] = src[index(0, 3)] ;
  dst[10] = src[index(2, 2)] ;
  dst[11] = src[index(3, 1)] ;
  dst[12] = src[index(1, 3)] ;
  dst[13] = src[index(3, 2)] ;
  dst[14] = src[index(2, 3)] ;
  dst[15] = src[index(3, 3)] ;
}
// reorder before 2 dimensional inverse transform, explicit "no scatter" version
void zfpx_unshuffle_2d(int32_t *src, int32_t *dst){
  dst[index(0, 0)] = src[ 0] ;
  dst[index(1, 0)] = src[ 1] ;
  dst[index(0, 1)] = src[ 2] ;
  dst[index(1, 1)] = src[ 3] ;
  dst[index(2, 0)] = src[ 4] ;
  dst[index(0, 2)] = src[ 5] ;
  dst[index(2, 1)] = src[ 6] ;
  dst[index(1, 2)] = src[ 7] ;
  dst[index(3, 0)] = src[ 8] ;
  dst[index(0, 3)] = src[ 9] ;
  dst[index(2, 2)] = src[10] ;
  dst[index(3, 1)] = src[11] ;
  dst[index(1, 3)] = src[12] ;
  dst[index(3, 2)] = src[13] ;
  dst[index(2, 3)] = src[14] ;
  dst[index(3, 3)] = src[15] ;
}
#undef index

// borrowed from the source code of zfp
// reorder table for 3 dimensional transform
#define index(i, j, k) ((i) + 4 * ((j) + 4 * (k)))
/* order coefficients (i, j, k) by i + j + k, then i^2 + j^2 + k^2 */
// cache_align_(static const uint8_t perm_3[64]) = {
cache_align_(const uint8_t zfpx_perm_3[64]) = {
  index(0, 0, 0), /*  0 : 0 */

  index(1, 0, 0), /*  1 : 1 */
  index(0, 1, 0), /*  2 : 1 */
  index(0, 0, 1), /*  3 : 1 */

  index(0, 1, 1), /*  4 : 2 */
  index(1, 0, 1), /*  5 : 2 */
  index(1, 1, 0), /*  6 : 2 */

  index(2, 0, 0), /*  7 : 2 */
  index(0, 2, 0), /*  8 : 2 */
  index(0, 0, 2), /*  9 : 2 */

  index(1, 1, 1), /* 10 : 3 */

  index(2, 1, 0), /* 11 : 3 */
  index(2, 0, 1), /* 12 : 3 */
  index(0, 2, 1), /* 13 : 3 */
  index(1, 2, 0), /* 14 : 3 */
  index(1, 0, 2), /* 15 : 3 */
  index(0, 1, 2), /* 16 : 3 */

  index(3, 0, 0), /* 17 : 3 */
  index(0, 3, 0), /* 18 : 3 */
  index(0, 0, 3), /* 19 : 3 */

  index(2, 1, 1), /* 20 : 4 */
  index(1, 2, 1), /* 21 : 4 */
  index(1, 1, 2), /* 22 : 4 */

  index(0, 2, 2), /* 23 : 4 */
  index(2, 0, 2), /* 24 : 4 */
  index(2, 2, 0), /* 25 : 4 */

  index(3, 1, 0), /* 26 : 4 */
  index(3, 0, 1), /* 27 : 4 */
  index(0, 3, 1), /* 28 : 4 */
  index(1, 3, 0), /* 29 : 4 */
  index(1, 0, 3), /* 30 : 4 */
  index(0, 1, 3), /* 31 : 4 */

  index(1, 2, 2), /* 32 : 5 */
  index(2, 1, 2), /* 33 : 5 */
  index(2, 2, 1), /* 34 : 5 */

  index(3, 1, 1), /* 35 : 5 */
  index(1, 3, 1), /* 36 : 5 */
  index(1, 1, 3), /* 37 : 5 */

  index(3, 2, 0), /* 38 : 5 */
  index(3, 0, 2), /* 39 : 5 */
  index(0, 3, 2), /* 40 : 5 */
  index(2, 3, 0), /* 41 : 5 */
  index(2, 0, 3), /* 42 : 5 */
  index(0, 2, 3), /* 43 : 5 */

  index(2, 2, 2), /* 44 : 6 */

  index(3, 2, 1), /* 45 : 6 */
  index(3, 1, 2), /* 46 : 6 */
  index(1, 3, 2), /* 47 : 6 */
  index(2, 3, 1), /* 48 : 6 */
  index(2, 1, 3), /* 49 : 6 */
  index(1, 2, 3), /* 50 : 6 */

  index(0, 3, 3), /* 51 : 6 */
  index(3, 0, 3), /* 52 : 6 */
  index(3, 3, 0), /* 53 : 6 */

  index(3, 2, 2), /* 54 : 7 */
  index(2, 3, 2), /* 55 : 7 */
  index(2, 2, 3), /* 56 : 7 */

  index(1, 3, 3), /* 57 : 7 */
  index(3, 1, 3), /* 58 : 7 */
  index(3, 3, 1), /* 59 : 7 */

  index(2, 3, 3), /* 60 : 8 */
  index(3, 2, 3), /* 61 : 8 */
  index(3, 3, 2), /* 62 : 8 */

  index(3, 3, 3), /* 63 : 9 */
};
// reorder after 3 dimensional forward transform, explicit "no gather" version
void zfpx_shuffle_3d(int32_t *src, int32_t *dst){
  dst[ 0] = src[index(0, 0, 0)] ;
  dst[ 1] = src[index(1, 0, 0)] ;
  dst[ 2] = src[index(0, 1, 0)] ;
  dst[ 3] = src[index(0, 0, 1)] ;
  dst[ 4] = src[index(0, 1, 1)] ;
  dst[ 5] = src[index(1, 0, 1)] ;
  dst[ 6] = src[index(1, 1, 0)] ;
  dst[ 7] = src[index(2, 0, 0)] ;
  dst[ 8] = src[index(0, 2, 0)] ;
  dst[ 9] = src[index(0, 0, 2)] ;
  dst[10] = src[index(1, 1, 1)] ;
  dst[11] = src[index(2, 1, 0)] ;
  dst[12] = src[index(2, 0, 1)] ;
  dst[13] = src[index(0, 2, 1)] ;
  dst[14] = src[index(1, 2, 0)] ;
  dst[15] = src[index(1, 0, 2)] ;
  dst[16] = src[index(0, 1, 2)] ;
  dst[17] = src[index(3, 0, 0)] ;
  dst[18] = src[index(0, 3, 0)] ;
  dst[19] = src[index(0, 0, 3)] ;
  dst[20] = src[index(2, 1, 1)] ;
  dst[21] = src[index(1, 2, 1)] ;
  dst[22] = src[index(1, 1, 2)] ;
  dst[23] = src[index(0, 2, 2)] ;
  dst[24] = src[index(2, 0, 2)] ;
  dst[25] = src[index(2, 2, 0)] ;
  dst[26] = src[index(3, 1, 0)] ;
  dst[27] = src[index(3, 0, 1)] ;
  dst[28] = src[index(0, 3, 1)] ;
  dst[29] = src[index(1, 3, 0)] ;
  dst[30] = src[index(1, 0, 3)] ;
  dst[31] = src[index(0, 1, 3)] ;
  dst[32] = src[index(1, 2, 2)] ;
  dst[33] = src[index(2, 1, 2)] ;
  dst[34] = src[index(2, 2, 1)] ;
  dst[35] = src[index(3, 1, 1)] ;
  dst[36] = src[index(1, 3, 1)] ;
  dst[37] = src[index(1, 1, 3)] ;
  dst[38] = src[index(3, 2, 0)] ;
  dst[39] = src[index(3, 0, 2)] ;
  dst[40] = src[index(0, 3, 2)] ;
  dst[41] = src[index(2, 3, 0)] ;
  dst[42] = src[index(2, 0, 3)] ;
  dst[43] = src[index(0, 2, 3)] ;
  dst[44] = src[index(2, 2, 2)] ;
  dst[45] = src[index(3, 2, 1)] ;
  dst[46] = src[index(3, 1, 2)] ;
  dst[47] = src[index(1, 3, 2)] ;
  dst[48] = src[index(2, 3, 1)] ;
  dst[49] = src[index(2, 1, 3)] ;
  dst[50] = src[index(1, 2, 3)] ;
  dst[51] = src[index(0, 3, 3)] ;
  dst[52] = src[index(3, 0, 3)] ;
  dst[53] = src[index(3, 3, 0)] ;
  dst[54] = src[index(3, 2, 2)] ;
  dst[55] = src[index(2, 3, 2)] ;
  dst[56] = src[index(2, 2, 3)] ;
  dst[57] = src[index(1, 3, 3)] ;
  dst[58] = src[index(3, 1, 3)] ;
  dst[59] = src[index(3, 3, 1)] ;
  dst[60] = src[index(2, 3, 3)] ;
  dst[61] = src[index(3, 2, 3)] ;
  dst[62] = src[index(3, 3, 2)] ;
  dst[63] = src[index(3, 3, 3)] ;
}
// reorder before 3 dimensional inverse transform, explicit "no scatter" version
void zfpx_unshuffle_3d(int32_t *src, int32_t *dst){
  dst[index(0, 0, 0)] = src[ 0] ;
  dst[index(1, 0, 0)] = src[ 1] ;
  dst[index(0, 1, 0)] = src[ 2] ;
  dst[index(0, 0, 1)] = src[ 3] ;
  dst[index(0, 1, 1)] = src[ 4] ;
  dst[index(1, 0, 1)] = src[ 5] ;
  dst[index(1, 1, 0)] = src[ 6] ;
  dst[index(2, 0, 0)] = src[ 7] ;
  dst[index(0, 2, 0)] = src[ 8] ;
  dst[index(0, 0, 2)] = src[ 9] ;
  dst[index(1, 1, 1)] = src[10] ;
  dst[index(2, 1, 0)] = src[11] ;
  dst[index(2, 0, 1)] = src[12] ;
  dst[index(0, 2, 1)] = src[13] ;
  dst[index(1, 2, 0)] = src[14] ;
  dst[index(1, 0, 2)] = src[15] ;
  dst[index(0, 1, 2)] = src[16] ;
  dst[index(3, 0, 0)] = src[17] ;
  dst[index(0, 3, 0)] = src[18] ;
  dst[index(0, 0, 3)] = src[19] ;
  dst[index(2, 1, 1)] = src[20] ;
  dst[index(1, 2, 1)] = src[21] ;
  dst[index(1, 1, 2)] = src[22] ;
  dst[index(0, 2, 2)] = src[23] ;
  dst[index(2, 0, 2)] = src[24] ;
  dst[index(2, 2, 0)] = src[25] ;
  dst[index(3, 1, 0)] = src[26] ;
  dst[index(3, 0, 1)] = src[27] ;
  dst[index(0, 3, 1)] = src[28] ;
  dst[index(1, 3, 0)] = src[29] ;
  dst[index(1, 0, 3)] = src[30] ;
  dst[index(0, 1, 3)] = src[31] ;
  dst[index(1, 2, 2)] = src[32] ;
  dst[index(2, 1, 2)] = src[33] ;
  dst[index(2, 2, 1)] = src[34] ;
  dst[index(3, 1, 1)] = src[35] ;
  dst[index(1, 3, 1)] = src[36] ;
  dst[index(1, 1, 3)] = src[37] ;
  dst[index(3, 2, 0)] = src[38] ;
  dst[index(3, 0, 2)] = src[39] ;
  dst[index(0, 3, 2)] = src[40] ;
  dst[index(2, 3, 0)] = src[41] ;
  dst[index(2, 0, 3)] = src[42] ;
  dst[index(0, 2, 3)] = src[43] ;
  dst[index(2, 2, 2)] = src[44] ;
  dst[index(3, 2, 1)] = src[45] ;
  dst[index(3, 1, 2)] = src[46] ;
  dst[index(1, 3, 2)] = src[47] ;
  dst[index(2, 3, 1)] = src[48] ;
  dst[index(2, 1, 3)] = src[49] ;
  dst[index(1, 2, 3)] = src[50] ;
  dst[index(0, 3, 3)] = src[51] ;
  dst[index(3, 0, 3)] = src[52] ;
  dst[index(3, 3, 0)] = src[53] ;
  dst[index(3, 2, 2)] = src[54] ;
  dst[index(2, 3, 2)] = src[55] ;
  dst[index(2, 2, 3)] = src[56] ;
  dst[index(1, 3, 3)] = src[57] ;
  dst[index(3, 1, 3)] = src[58] ;
  dst[index(3, 3, 1)] = src[59] ;
  dst[index(2, 3, 3)] = src[60] ;
  dst[index(3, 2, 3)] = src[61] ;
  dst[index(3, 3, 2)] = src[62] ;
  dst[index(3, 3, 3)] = src[63] ;
}
#undef index

// reversible inverse lifting transform, in place, 1 point
void zfpx_inv_lift_1(int32_t* x, int32_t *y, int32_t* z, int32_t* w){
  w[0] += z[0];
  z[0] += y[0]; w[0] += z[0];
  y[0] += x[0]; z[0] += y[0]; w[0] += z[0];
}

// reversible inverse lifting transform, in place, 4 points
void zfpx_inv_lift_4(int32_t* x, int32_t* y, int32_t* z, int32_t* w){
  int i;
  for(i=0 ; i<4 ; i++){
    w[i] += z[i];
    z[i] += y[i]; w[i] += z[i];
    y[i] += x[i]; z[i] += y[i]; w[i] += z[i];
  }
}

// reversible inverse lifting transform, in place, 16 points
void zfpx_inv_lift_16(int32_t* x, int32_t* y, int32_t* z, int32_t* w){
  int i;
  for(i=0 ; i<16 ; i++){
    w[i] += z[i];
    z[i] += y[i]; w[i] += z[i];
    y[i] += x[i]; z[i] += y[i]; w[i] += z[i];
  }
}

// reversible forward lifting transform, in place, 1 point
void zfpx_fwd_lift_1(int32_t* x, int32_t* y, int32_t* z, int32_t* w){
  w[0] -= z[0]; z[0] -= y[0]; y[0] -= x[0];
  w[0] -= z[0]; z[0] -= y[0];
  w[0] -= z[0];
}

// reversible forward lifting transform, in place, 4 points
void zfpx_fwd_lift_4(int32_t* x, int32_t* y, int32_t* z, int32_t* w){
  int i;
  for(i=0 ; i<4 ; i++){
    w[i] -= z[i]; z[i] -= y[i]; y[i] -= x[i];
    w[i] -= z[i]; z[i] -= y[i];
    w[i] -= z[i];
  }
}

// reversible forward lifting transform, in place, 16 points
void zfpx_fwd_lift_16(int32_t* x, int32_t* y, int32_t* z, int32_t* w){
  int i;
  for(i=0 ; i<16 ; i++){
    w[i] -= z[i]; z[i] -= y[i]; y[i] -= x[i];
    w[i] -= z[i]; z[i] -= y[i];
    w[i] -= z[i];
  }
}

// reversible inverse lifting transform, in place, 1 dimensional data (4)
void zfpx_inv_xform_1d(int32_t *t){
  zfpx_inv_lift_1(t, t+1, t+2, t+3) ;
}

// reversible inverse lifting transform, in place, 2 dimensional data (4,4)
void zfpx_inv_xform_2d(int32_t *t0){
  int32_t x[4], y[4], z[4], w[4], t[16] ;

  zfpx_unshuffle_2d(t0, t) ;  // unshuffle input
//   for(i=0 ; i<16 ; i++) t[zfpx_perm_2[i]] = t0[i] ;  // unshuffle input
  // inverse transform along y
  zfpx_inv_lift_4(t, t+4, t+8, t+12) ;
  // inverse transform along x
  x[0] = t[ 0] ; x[1] = t[ 4] ; x[2] = t[ 8] ; x[3] = t[12] ;
  y[0] = t[ 1] ; y[1] = t[ 5] ; y[2] = t[ 9] ; y[3] = t[13] ;
  z[0] = t[ 2] ; z[1] = t[ 6] ; z[2] = t[10] ; z[3] = t[14] ;
  w[0] = t[ 3] ; w[1] = t[ 7] ; w[2] = t[11] ; w[3] = t[15] ;
  zfpx_inv_lift_4(x, y, z, w) ;
  t0[ 0] = x[0] ; t0[ 4] = x[1] ; t0[ 8] = x[2] ; t0[12] = x[3] ;
  t0[ 1] = y[0] ; t0[ 5] = y[1] ; t0[ 9] = y[2] ; t0[13] = y[3] ;
  t0[ 2] = z[0] ; t0[ 6] = z[1] ; t0[10] = z[2] ; t0[14] = z[3] ;
  t0[ 3] = w[0] ; t0[ 7] = w[1] ; t0[11] = w[2] ; t0[15] = w[3] ;
}

// reversible inverse lifting transform, in place, 3 dimensional data (4,4,4)
void zfpx_inv_xform_3d(int32_t *t0){
  int32_t x[16], y[16], z[16], w[16], t[64] ;
  int i ;
  zfpx_unshuffle_3d(t0, t) ;  // unshuffle input
//   for(i=0 ; i<64 ; i++) t[zfpx_perm_3[i]] = t0[i] ;  // unshuffle input
  // inverse transform along z (stride 1, length = 16)
  zfpx_inv_lift_16(t, t+16, t+32, t+48) ;
  // inverse transform along y (stride 1, lenth = 4)
  zfpx_inv_lift_4(t   , t+ 4, t+ 8, t+12) ;  // plane 0
  zfpx_inv_lift_4(t+16, t+20, t+24, t+28) ;  // plane 1
  zfpx_inv_lift_4(t+32, t+36, t+40, t+44) ;  // plane 2
  zfpx_inv_lift_4(t+48, t+52, t+56, t+60) ;  // plane 3
  // inverse transform along x (stride 4, length = 16)
  for(i=0 ; i<16 ; i++){
    x[i] = t[4*i] ; y[i] = t[1+4*i] ; z[i] = t[2+4*i] ; w[i] = t[3+4*i] ;
  }
  zfpx_inv_lift_16(x, y, z, w) ;
  for(i=0 ; i<16 ; i++){
    t0[4*i] = x[i] ; t0[1+4*i] = y[i] ; t0[2+4*i] = z[i] ; t0[3+4*i] = w[i] ;
  }
}

// reversible forward lifting transform, in place, 1 dimensional data (4)
void zfpx_fwd_xform_1d(int32_t *t){
  zfpx_fwd_lift_1(t, t+1, t+2, t+3) ;
}

static int xxx = 0 ;
// reversible forward lifting transform, in place, 2 dimensional data (4,4)
void zfpx_fwd_xform_2d(int32_t *t){
  int32_t x[4], y[4], z[4], w[4], t1[16] ;
  // transform along y
  zfpx_fwd_lift_4(t, t+4, t+8, t+12) ;
  // transform along x
  x[0] = t[ 0] ; x[1] = t[ 4] ; x[2] = t[ 8] ; x[3] = t[12] ;
  y[0] = t[ 1] ; y[1] = t[ 5] ; y[2] = t[ 9] ; y[3] = t[13] ;
  z[0] = t[ 2] ; z[1] = t[ 6] ; z[2] = t[10] ; z[3] = t[14] ;
  w[0] = t[ 3] ; w[1] = t[ 7] ; w[2] = t[11] ; w[3] = t[15] ;
  zfpx_fwd_lift_4(x, y, z, w) ;
  t1[ 0] = x[0] ; t1[ 4] = x[1] ; t1[ 8] = x[2] ; t1[12] = x[3] ;
  t1[ 1] = y[0] ; t1[ 5] = y[1] ; t1[ 9] = y[2] ; t1[13] = y[3] ;
  t1[ 2] = z[0] ; t1[ 6] = z[1] ; t1[10] = z[2] ; t1[14] = z[3] ;
  t1[ 3] = w[0] ; t1[ 7] = w[1] ; t1[11] = w[2] ; t1[15] = w[3] ;
  zfpx_shuffle_2d(t1, t) ; // shuffle output
//   for(int i=0 ; i<16 ; i++) t[i] = t1[i] ;  // do not shuffle
// if(xxx <4){
//   int i ;
//   for(i=0 ; i<16 ; i++) printf("%d,",t1[i]) ; printf("\n") ;
//   for(i=0 ; i<16 ; i++) printf("%d,",t[i]) ; printf("\n") ;
//   xxx++ ;
// }
//   for(int i=0 ; i<16 ; i++) t[i] = t1[zfpx_perm_2[i]] ; // shuffle output
}

// reversible forward lifting transform, in place, 3 dimensional data (4,4,4)
void zfpx_fwd_xform_3d(int32_t *t){
  int32_t x[16], y[16], z[16], w[16], t1[64] ;
  int i ;
// printf("--");
// for(i=0 ; i<64 ; i++) printf("%d ",t[i]) ; printf("\n");
  // transform along z (stride 1, length = 16)
  zfpx_fwd_lift_16(t, t+16, t+32, t+48) ;
  // transform along y (stride 1, lenth = 4)
  zfpx_fwd_lift_4(t   , t+ 4, t+ 8, t+12) ;  // plane 0
  zfpx_fwd_lift_4(t+16, t+20, t+24, t+28) ;  // plane 1
  zfpx_fwd_lift_4(t+32, t+36, t+40, t+44) ;  // plane 2
  zfpx_fwd_lift_4(t+48, t+52, t+56, t+60) ;  // plane 3
  // transform along x (stride 4, length = 16)
  for(i=0 ; i<16 ; i++){
    x[i] = t[4*i] ; y[i] = t[1+4*i] ; z[i] = t[2+4*i] ; w[i] = t[3+4*i] ;
  }
  zfpx_fwd_lift_16(x, y, z, w) ;
  for(i=0 ; i<16 ; i++){
    t1[4*i] = x[i] ; t1[1+4*i] = y[i] ; t1[2+4*i] = z[i] ; t1[3+4*i] = w[i] ;
  }
  zfpx_shuffle_3d(t1,t) ; // shuffle output
//   for(i=0 ; i<64 ; i++) t[i] = t1[i] ; // do not shuffle
// printf("++");
// for(i=0 ; i<64 ; i++) printf("%d ",t[i]) ; printf("\n");
//   for(i=0 ; i<64 ; i++) t[i] = t1[zfpx_perm_3[i]] ; // shuffle output
}

#if 0
#define NBMASK 0xaaaaaaaau /* negabinary<-> 2's complement binary conversion mask */

// signed integer (2's complement) to negabinary (base -2) conversion
uint32_t int_to_negabinary(int32_t x)
{
  return ((uint32_t)x + NBMASK) ^ NBMASK;
}

void v_int_to_negabinary(int32_t x, int32_t n)
{
  int i ;
  for(i=0 ; i<n ; i++) x = int_to_negabinary(x) ;
}

// negabinary (base -2) to signed integer (2's complement) conversion
int32_t negabinary_to_int(uint32_t x)
{
  return (int32_t)((x ^ NBMASK) - NBMASK);
}

void v_negabinary_to_int(uint32_t x, int32_t n)
{
  int i ;
  for(i=0 ; i<n ; i++) x = negabinary_to_int(x) ;
}
#endif

// get the IEEE exponent of the largest float (absolute value) in float array f
uint32_t get_ieee_emax(float *f, int n){
  int i ;
  uint32_t *uf = (uint32_t *) f ;
  uint32_t t, emax = 0 ;

  for(i=0 ; i<n ; i++){
    t = uf[i] << 1 ;
    emax = (t > emax) ? t : emax ;
  }
  return emax >> 24 ;  // return exponent WITH IEEE bias (127)
}

// quantize float array f by normalizing all exponents to the largest exponent emax
// if emax == 0, compute said exponent from float array f
// the sign is preserved in the integer quantized value
// the largest float (inmagnitude) will have a 30 bit "mantissa"
// the "hidden 1" is restored
float zfpx_quantize(float *f, int32_t *fi, int n, uint32_t emax){
  int i ;
  union{
    float f    ;
    uint32_t u ;
  } factor ;
  float rfactor ;
  if(emax == 0) emax = get_ieee_emax(f, n) ;

  factor.u = (157 - emax + 127) << 23 ;
  rfactor = 1.0 / factor.f ;
  for(i=0 ; i<n ; i++) fi[i] = f[i] * factor.f ;

  return rfactor ;
}

// find restore factor associated with max IEEE exponent emax (emax includes IEEE bias of 127)
float zfpx_rfactor(uint32_t emax){
  union{
    float f    ;
    uint32_t u ;
  } factor ;
  factor.u = (157 - emax + 127) << 23 ;
  return 1.0 / factor.f ;
}

// the following 6 functions extract/insert "fake cubes" from/to a 2 dimensional
// array for compression purposes
//    ^
// (J)|
//    +-------------+
// 15 |03         33|
//    |             |
//    |    k = 3    |
// 12 |00         30|
//    +-------------+
// 11 |03         33|
//    |             |
//    |    k = 2    |
//  8 |00         30|
//    +-------------+
//  7 |03         33|
//    |             |
//    |    k = 1    |
//  4 |00         30|
//    +-------------+
//  3 |03         33|
//    |             |
//    |    k = 0    |
//  0 |00         30|
//    +-------------+------->
//     0           3     (I)
// get 4 elements from 16 rows of f and shape them as a 4x4x4 cube
// 4x16 -> 4x4x4
void zfpx_get_x4y4y4_block(int32_t *f, int32_t lni, int32_t *block){
  int i, j ;
  for(j=0 ; j<16 ; j++){                     // 16 rows
    for(i=0 ; i<4 ; i++) block[i] = f[i] ;   // 4 elements from each row
    f += lni ;
    block += 4 ;
  }
}
// inverse of zfpx_get_x4y4y4_block
void zfpx_put_x4y4y4_block(int32_t *f, int32_t lni, int32_t *block){
  int i, j ;
  for(j=0 ; j<16 ; j++){                     // 16 rows
    for(i=0 ; i<4 ; i++) f[i] = block[i] ;   // 4 elements from each row
    f += lni ;
    block += 4 ;
  }
}

//   ^
//(J)|
//   +-------------+-------------+-------------+-------------+
// 3 |03         33|03         33|03         33|           33| f3[0-15]
//   |             |             |             |             | f2[0-15]
//   |    k = 0    |    k = 1    |    k = 2    |    k = 3    | f1[0-15]
// 0 |00         30|00           |00         30|00         30| f0[0-15]
//   +-------------+-------------+-------------+-------------+------->
//    0           3 4           7 8          11 12         15    (I)
// get 16 elements from 4 rows of f and shape them as a 4x4x4 cube
// 16x4 -> 4x4x4 (4 adjacent blocks of 4x4
void zfpx_get_x4y4x4_block(int32_t *f, int32_t lni, int32_t *block){
  int i, j ;
  int32_t *f0, *f1, *f2, *f3 ;
  f0 = f ; f1 = f0 + lni ; f2 = f1 + lni ; f3 = f2 + lni ;
  for(j=0 ; j<4 ; j++) {                     // 4 blocks
    for(i=0 ; i<4 ; i++) {                   // 16 elements from each block
      block[i   ] = f0[i+4*j] ;
      block[i+ 4] = f1[i+4*j] ;
      block[i+ 8] = f2[i+4*j] ;
      block[i+12] = f3[i+4*j] ;
    }
    block += 16 ;                             // next block
  }
}
// inverse of zfpx_get_x4y4x4_block
void zfpx_put_x4y4x4_block(int32_t *f, int32_t lni, int32_t *block){
  int i, j ;
  int32_t *f0, *f1, *f2, *f3 ;
  f0 = f ; f1 = f0 + lni ; f2 = f1 + lni ; f3 = f2 + lni ;
  for(j=0 ; j<4 ; j++) {                     // 4 blocks
    for(i=0 ; i<4 ; i++) {                   // 16 elements from each block
      f0[i+4*j] = block[i   ] ;
      f1[i+4*j] = block[i+ 4] ;
      f2[i+4*j] = block[i+ 8] ;
      f3[i+4*j] = block[i+12] ;
    }
    block += 16 ;                             // next block
  }
}

//   ^
//(J)|
//   +-------------+-------------+-------------+-------------+
// 3 |003       303|013       313|023       323|033       333| plane 3
//   |002       302|012       312|022       322|032       332| plane 2
//   |001       301|011       311|021       321|031       331| plane 1
// 0 |000       300|010       310|020       320|030       330| plane 0
//   +-------------+-------------+-------------+-------------+------->
//    0           3 4           7 8          11 12         15    (I)
// get 16 elements from 4 rows of f and shape them as a 4x4x4 cube
// the 16 elements from a row are shaped as a 4x4 plane
// the 4 rows provide the 4 planes
void zfpx_get_x4x4y4_block(int32_t *f, int32_t lni, int32_t *block){
  int i, j ;
  for(j=0 ; j<4 ; j++){                       // 4 rows
    for(i=0 ; i<16 ; i++) block[i] = f[i] ;   // 16 elements from each row (4x4 plane)
    f += lni ;
    block += 16 ;
  }
}
// inverse of zfpx_get_x4x4y4_block
void zfpx_put_x4x4y4_block(int32_t *f, int32_t lni, int32_t *block){
  int i, j ;
  for(j=0 ; j<4 ; j++){                       // 4 rows
    for(i=0 ; i<16 ; i++) f[i] = block[i] ;   // 16 elements from each row (4x4 plane)
    f += lni ;
    block += 16 ;
  }
}

// transform a 64 x 64 block into 64 4x4x4 slices
void zfpx_gather_64_64(int32_t *f, int32_t lni, int32_t *blocks, int transform){
  int i0, j0;
  // 4*lni = storage length for 4 rows
  // 16 slices of 4 rows, 16 elements in each row => a 4x4x4 block
  for(j0=0 ; j0<61*lni ; j0+=4*lni){    // storage length of rows is lni
    for(i0=0 ; i0<49 ; i0+=16){
//       zfpx_get_x4x4y4_block(f+i0+j0, lni, blocks) ;
      zfpx_get_x4y4x4_block(f+i0+j0, lni, blocks) ;
      if(transform == 1) zfpx_fwd_xform_3d(blocks) ; // forward 3D transform and shuffle
      if(transform == 4) {
        zfpx_fwd_xform_2d(blocks) ;
        zfpx_fwd_xform_2d(blocks+16) ;
        zfpx_fwd_xform_2d(blocks+32) ;
        zfpx_fwd_xform_2d(blocks+48) ;
      }
      blocks += 64 ;
    }
  }
}
// inverse of zfpx_gather_64_64
void zfpx_scatter_64_64(int32_t *f, int32_t lni, int32_t *blocks, int transform){
  int i0, j0;
  // 4*lni = storage length for 4 rows
  // 16 slices of 4 rows, 16 elements in each row => a 4x4x4 block
  for(j0=0 ; j0<61*lni ; j0+=4*lni){    // storage length of rows is lni
    for(i0=0 ; i0<49 ; i0+=16){
      if(transform) zfpx_inv_xform_3d(blocks) ; // inverse 3D transform and unshuffle
//       zfpx_put_x4x4y4_block(f+i0+j0, lni, blocks) ;
      zfpx_put_x4y4x4_block(f+i0+j0, lni, blocks) ;
      blocks += 64 ;
    }
  }
}

// 16 32 bit integers -> 32 bit planes (lower 16 bits in 64)
// the first element of src ( 0) will end up in the LSB column in plabes (bit 0)
// element N will end up in column N                                     (bit N)
// the last  element of src (15) will end up in the MSB column in planes (bit 15)
// planes[0 ] is the Most  Significant Bit Plane
// planes[31] is the Least Significant Bit Plane
// bits are numbered right to left
void zfpx_bit_plane_32_16(uint32_t *src, uint16_t *planes){
  int i ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  uint64_t plane ;
  __m256i v0, v1 ;
  v0 = _mm256_loadu_si256((__m256i const *)(src   )) ; // load 16 values
  v1 = _mm256_loadu_si256((__m256i const *)(src+ 8)) ;
  for(i=0 ; i<32 ; i++){
    plane = 0 ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v1) ; v1 = _mm256_slli_epi32(v1, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v0) ; v0 = _mm256_slli_epi32(v0, 1) ;
    planes[i] = plane ;
  }
#else
  int j ;
  uint32_t bits ;
  for(i=0 ; i<32 ; i++) planes[i] = 0l ;
  for(j=15 ; j>=0 ; j--){
    bits = src[j] ;
    for(i=31 ; i>=0 ; i--){
      planes[i] = (planes[i] << 1) | (bits & 1) ; bits >>= 1 ;
    }
  }
//   printf("BEEP16 !!\n");
#endif
}

// 32 32 bit integers -> 32 bit planes
// the first element of src ( 0) will end up in the LSB column in plabes (bit 0)
// element N will end up in column N                                     (bit N)
// the last  element of src (31) will end up in the MSB column in planes (bit 31)
// planes[0 ] is the Most  Significant Bit Plane
// planes[31] is the Least Significant Bit Plane
// bits are numbered right to left
void zfpx_bit_plane_32_32(uint32_t *src, uint32_t *planes){
  int i ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  uint64_t plane ;
  __m256i v0, v1 ;
  v0 = _mm256_loadu_si256((__m256i const *)(src   )) ; // load 16 values
  v1 = _mm256_loadu_si256((__m256i const *)(src+ 8)) ;
  v2 = _mm256_loadu_si256((__m256i const *)(src+16)) ;
  v3 = _mm256_loadu_si256((__m256i const *)(src+24)) ;
  for(i=0 ; i<32 ; i++){
    plane = 0 ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v3) ; v3 = _mm256_slli_epi32(v3, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v2) ; v2 = _mm256_slli_epi32(v2, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v1) ; v1 = _mm256_slli_epi32(v1, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v0) ; v0 = _mm256_slli_epi32(v0, 1) ;
    planes[i] = plane ;
  }
#else
  int j ;
  uint32_t bits ;
  for(i=0 ; i<32 ; i++) planes[i] = 0l ;
  for(j=31 ; j>=0 ; j--){
    bits = src[j] ;
    for(i=31 ; i>=0 ; i--){
      planes[i] = (planes[i] << 1) | (bits & 1) ; bits >>= 1 ;
    }
  }
//   printf("BEEP16 !!\n");
#endif
}

// 64 32 bit integers -> 32 bit planes (64 bits)
// the first element of src (index  0) will end up in the LSB column in planes (bit 0)
// element N (index N) will end up in column N                                 (bit N)
// the last  element of src (index 63) will end up in the MSB column in planes (bit 63)
// planes[0 ] is the Most  Significant Bit Plane
// planes[63] is the Least Significant Bit Plane
// bits are numbered right to left
void zfpx_bit_plane_32_64(uint32_t *src, uint64_t *planes){
  int i ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  uint64_t plane ;
  __m256i v0, v1, v2, v3, v4, v5, v6, v7 ;
  v0 = _mm256_loadu_si256((__m256i const *)(src   )) ; // load 64 values
  v1 = _mm256_loadu_si256((__m256i const *)(src+ 8)) ;
  v2 = _mm256_loadu_si256((__m256i const *)(src+16)) ;
  v3 = _mm256_loadu_si256((__m256i const *)(src+24)) ;
  v4 = _mm256_loadu_si256((__m256i const *)(src+32)) ;
  v5 = _mm256_loadu_si256((__m256i const *)(src+40)) ;
  v6 = _mm256_loadu_si256((__m256i const *)(src+48)) ;
  v7 = _mm256_loadu_si256((__m256i const *)(src+56)) ;
  for(i=0 ; i<32 ; i++){
    plane = 0 ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v7) ; v7 = _mm256_slli_epi32(v7, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v6) ; v6 = _mm256_slli_epi32(v6, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v5) ; v5 = _mm256_slli_epi32(v5, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v4) ; v4 = _mm256_slli_epi32(v4, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v3) ; v3 = _mm256_slli_epi32(v3, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v2) ; v2 = _mm256_slli_epi32(v2, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v1) ; v1 = _mm256_slli_epi32(v1, 1) ;
    plane <<= 8 ; plane |= _mm256_movemask_ps((__m256) v0) ; v0 = _mm256_slli_epi32(v0, 1) ;
    planes[i] = plane ;
  }
#else
  int j ;
  uint32_t bits ;
  for(i=0 ; i<32 ; i++) planes[i] = 0l ;
  for(j=63 ; j>=0 ; j--){
    bits = src[j] ;
    for(i=31 ; i>=0 ; i--){
      planes[i] = (planes[i] << 1) | (bits & 1) ; bits >>= 1 ;
    }
  }
//   printf("BEEP 64!!\n");
#endif
}

#if defined(SELF_TEST)
int main(int argc, char **argv){
  int32_t fi[17] ;
  float ff[17] ;
  uint32_t *fu = (uint32_t *) ff ;
  int i, j, k, errors ;
  float factor ;
  int mask, offset, error, errabs, errmax = 0, errneg = 0, errpos = 0, nbits = 0 ;
  int64_t bias, absbias, errtot = 0, errtotabs = 0, deltaerr ;
  uint32_t ineg, ipos, npts = 0 ;
  uint64_t errplus = 0, errminus = 0 ;
  int32_t t2d[16] ;
  int32_t t3d[64] ;
  uint32_t stream[64], stream0 ;
  uint16_t planes[32], planes16[32] ;
  char string[65] ;
  int32_t coord[32][32] ;
  int32_t *coord1 = &coord[0][0] ;
  int32_t block[64] ;
#define BLKSIDE 67
  int32_t blocks0[BLKSIDE*BLKSIDE] ; // original
  int32_t blocks1[64*64] ;           // gathered
  int32_t blocks2[BLKSIDE*BLKSIDE] ; // restored

  for(i=0 ; i<sizeof(blocks0)/sizeof(int32_t) ; i++) { blocks0[i] = i*i+11 ; blocks2[i] = blocks0[i] ; }
  for(i=0 ; i<sizeof(blocks1)/sizeof(int32_t) ; i++) blocks1[i] = 0 ;
  zfpx_gather_64_64(blocks0, BLKSIDE, blocks1, 1) ;
  for(j=0 ; j<64 ; j+=4){
    for(i=0 ; i< 64 ; i++) printf("%d ", blocks1[i+64*j]) ; printf("\n") ;
  }
  for(j=0 ; j<64 ; j++){                                // set part of blocks2 to 0
    for(i=0 ; i<64 ; i++){
      blocks2[i+j*BLKSIDE] = 0;
    }
  }
  zfpx_scatter_64_64(blocks2, BLKSIDE, blocks1, 1) ;    // re-create same part of block2
  errors = 0 ;
  for(i=0 ; i<sizeof(blocks0)/sizeof(int32_t) ; i++){
    if(blocks0[i] != blocks2[i]) errors++ ;
  }
  printf("block0->block1->block2 : errors = %d\n", errors) ;
return 0 ;
#if 1
  for(j=0 ; j<32 ; j++){
    for(i=0 ; i<32 ; i++){
      coord[j][i] = (i << 8) + j ;
    }
  }
  for(j=15 ; j>=0 ; j--){
    for(i=0 ; i<16 ; i++) printf("%5.4x", coord1[i + j*32]) ; printf("\n") ;
  }
  printf("\n") ;
  printf("===== x4y4x4 (nx,4,ny/4) =====\n") ;
  zfpx_get_x4y4x4_block(coord1, 32, block) ;
  for(k = 3 ; k >= 0 ; k--){
    for(j=3 ; j>=0 ; j--){
      for(i=0 ; i<4 ; i++){
        printf("%5.4x", block[i + j*4 + k*16]) ;
      }
      printf("\n") ;
    }
    printf("\n") ;
  }
  printf("===== x4y4y4 (4,nx/4,4,ny/4) =====\n") ;
  zfpx_get_x4y4y4_block(coord1, 32, block) ;
  for(k = 3 ; k >= 0 ; k--){
    for(j=3 ; j>=0 ; j--){
      for(i=0 ; i<4 ; i++){
        printf("%5.4x", block[i + j*4 + k*16]) ;
      }
      printf("\n") ;
    }
    printf("\n") ;
  }
  printf("===== x4x4y4 (4,nx/4,ny) =====\n") ;
  zfpx_get_x4x4y4_block(coord1, 32, block) ;
  for(k = 3 ; k >= 0 ; k--){
    for(j=3 ; j>=0 ; j--){
      for(i=0 ; i<4 ; i++){
        printf("%5.4x", block[i + j*4 + k*16]) ;
      }
      printf("\n") ;
    }
    printf("\n") ;
  }
  return 0 ;
#endif
  if(argc > 1){
    sscanf(argv[1], "%d", &nbits) ;
    mask = (-1) << nbits ;
    offset = ~mask >> 1 ;
    printf("mask = %8.8x, offset = %8.8x\n", mask, offset) ;
    for(i=-0x7FFFFFFF ; i < 0x7FFFFFFF-10 ; i+=11) {
      if(argc > 2) {
        ineg = int_to_negabinary(i) ;
        ineg &= mask ;
        ipos = negabinary_to_int(ineg) ;
      }else{
        ipos = (i + offset) & mask ;
      }
      error = ipos - i ; errabs = (error <0) ? -error : error ;
      errtot += error ;
      errplus  += ((error > 0) ? 1 : 0) ;
      errminus += ((error < 0) ? 1 : 0) ;
      errtotabs += errabs ;
      errmax = (errabs > errmax) ? errabs : errmax ;
      errneg = (error < errneg) ? error : errneg ;
      errpos = (error > errpos) ? error : errpos ;
      npts ++ ;
    }
    bias = errtot / npts ; absbias = errtotabs / npts ;
    deltaerr = errplus ; deltaerr -= errminus ;
    printf("npts = %u, errors(>/</>-<)(%lu/%lu/%ld), errmax(abs/pos/neg) = (%d/%d/%d), bias = %ld, avgabserr = %ld\n",
           npts, errplus, errminus, deltaerr, errmax, errpos, errneg, bias, absbias) ;
  }else{
    for(i=0 ; i<17 ; i++) ff[i] = 1.0001 * (i - 8) ;
    uint32_t emax = get_ieee_emax(ff, 17) ;
    printf("emax = %u\n", emax) ;
    printf("ff ");
    for(i=0 ; i<17 ; i++) printf(" %8.5f", ff[i]) ; printf("\n") ;
    printf("ff ");
    for(i=0 ; i<17 ; i++) printf(" %8.8x", fu[i]) ; printf("\n") ;
    factor = zfpx_quantize(ff, fi, 17, emax) ;
    printf("fi ");
    for(i=0 ; i<17 ; i++) printf(" %8.8x", fi[i]) ; printf("\n") ;
    for(i=0 ; i<17 ; i++) fi[i] = int_to_negabinary(fi[i]) ;
    printf("fn ");
    for(i=0 ; i<17 ; i++) printf(" %8.8x", fi[i]) ; printf("\n") ;
    for(i=0 ; i<17 ; i++) fi[i] &= 0xFFFFC000u ;
    printf("fn ");
    for(i=0 ; i<17 ; i++) printf(" %8.8x", fi[i]) ; printf("\n") ;
    for(i=0 ; i<17 ; i++) fi[i] = negabinary_to_int(fi[i]) ;
    printf("fi ");
    for(i=0 ; i<17 ; i++) printf(" %8.8x", fi[i]) ; printf("\n") ;
    printf("ff ");
    for(i=0 ; i<17 ; i++) printf(" %8.5f", fi[i]*factor) ; printf("\n") ;
    for(i=0 ; i<17 ; i++) ff[i] = fi[i]*factor ;
    printf("ff ");
    for(i=0 ; i<17 ; i++) printf(" %8.8x", fu[i]) ; printf("\n") ;
  }
  printf("t1d  = ");
  for(i=0 ; i<4 ; i++) { t2d[i] = i ; printf("%3d", t2d[i]) ; }
  printf("\n");
//   rev_fwd_lift_Int(t2d, 1);
  zfpx_fwd_xform_1d(t2d) ;
  printf("t1d' = ");
  for(i=0 ; i<4 ; i++) { printf("%3d", t2d[i]) ; }
  printf("\n");
//   rev_inv_lift_Int(t2d, 1);
  zfpx_inv_xform_1d(t2d) ;
  printf("t1d  = ");
  for(i=0 ; i<4 ; i++) { printf("%3d", t2d[i]) ; }
  printf("\n");

  printf("t2d  = ");
  for(j=0 ; j<4 ; j++){
    for(i=0 ; i<4 ; i++){
      t2d[j*4+i] = i+j+10 ;
      printf("%3d", t2d[j*4+i]) ;
    }
  }
  printf("\n");
  zfpx_fwd_xform_2d(t2d) ;
  printf("t2d' = ");
  for(j=0 ; j<4 ; j++){
    for(i=0 ; i<4 ; i++){
      printf("%3d", t2d[j*4+i]) ;
    }
  }
  printf("\n");
  zfpx_inv_xform_2d(t2d) ;
  printf("t2d  = ");
  for(j=0 ; j<4 ; j++){
    for(i=0 ; i<4 ; i++){
      printf("%3d", t2d[j*4+i]) ;
    }
  }
  printf("\n");

  printf("t3d  :\n");
  for(k=0 ; k<4 ; k++){
    for(j=0 ; j<4 ; j++){
      for(i=0 ; i<4 ; i++){
//         t3d[k*16+j*4+i] = i+j+k+10 ;
        t3d[k*16+j*4+i] = (i+1)*(j+1)*(k+1)+20+i+j+k ;
        printf("%3d", t3d[k*16+j*4+i]) ;
      }
      printf("  ");
    }
    printf("\n");
  }

  zfpx_fwd_xform_3d(t3d) ;
  printf("t3d' :\n");
  for(k=0 ; k<4 ; k++){
    for(j=0 ; j<4 ; j++){
      for(i=0 ; i<4 ; i++){
        printf("%3d", t3d[k*16+j*4+i]) ;
      }
      printf("  ");
    }
    printf("\n");
  }
//   for(i=0 ; i<64 ; i++){ printf("%3d", t3d[i]) ; } printf("\n");

  zfpx_inv_xform_3d(t3d) ;
  printf("t3d  :\n");
  for(k=0 ; k<4 ; k++){
    for(j=0 ; j<4 ; j++){
      for(i=0 ; i<4 ; i++){
        printf("%3d", t3d[k*16+j*4+i]) ;
      }
      printf("  ");
    }
    printf("\n");
  }
#if 1
  stream0 = 0xCC000000u ;
  for(i=0 ; i<32 ; i++){
    stream[i] = stream0 ; stream[63-i] = stream0 ; stream0 >>= 1 ;
  }
  stream[15] = 0x55555555 ;
  stream[63] = 0xAAAAAAAA ;
  for(i=0 ; i<31 ; i++) {planes[i] = 0 ; planes16[i] = 0 ; } ;
  zfpx_bit_plane_32_64(stream, planes) ;    // 4x4x4 array
  zfpx_bit_plane_32_16(stream, planes16) ;  // 4x4 array
  for(i=0 ; i< 64 ; i++) {
    string[32] = 0;
    for(j=0 ; j<32 ; j++) string[31-j] = (stream[i] & (1 << j)) ? '1' : '0' ;
    printf(" %8.8x %s", stream[i], string) ;
    if(i<32){
      printf(" |");
      string[64] = 0;
      for(j=0 ; j<64 ; j++) string[63-j] = (planes[i] & (1l << j)) ? '1' : '0' ;
      printf(" %s", string) ;
    }
    if(i<32){
      printf(" |");
      string[16] = 0;
      for(j=0 ; j<16 ; j++) string[63-j] = (planes16[i] & (1l << j)) ? '1' : '0' ;
      printf(" %s", string) ;
    }
    printf("\n") ;
  }
#endif
  for(i=0 ; i<16 ; i++) stream[i] = 1 << i ;
}
#endif
