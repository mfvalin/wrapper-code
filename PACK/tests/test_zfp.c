//  Copyright (C) 2022  Environnement Canada
//
//  This is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation,
//  version 2.1 of the License.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <rmn/misc_timers.h>
#include <rmn/misc_types.h>
#include <misc_pack.h>
#include <rmn/misc_operators.h>
#include <rmn/lorenzo.h>
#include <misc_zfp.h>
#include <misc_analyze.h>

static inline uint32_t index2(uint32_t i, uint32_t j, uint32_t ni, uint32_t nj){
  return i + ni * j ;
}

static inline uint32_t index3(uint32_t i, uint32_t j, uint32_t k, uint32_t ni, uint32_t nj, uint32_t nk){
  return i + ni * j ;
}

int main(int argc, char **argv){
  /* allocate array of floats */
  int nx = 2048; int nx2 = nx/2 ;
  int ny = 1024; int ny2 = ny/2 ;
  int nz = 1; int nz2 = 4 ;
  int nxy2 = nx2 * ny2;
  float* array = malloc(nx * ny * nz * sizeof(float));
  int32_t* qarray = malloc(nx * ny * nz * sizeof(int32_t));  // quantized floats
  int32_t* larray = malloc(nx * ny * nz * sizeof(int32_t));  // lorenzo prediction
  float* Array = malloc(nx * ny * nz * sizeof(float));
  float* brray = malloc(nx * ny * nz * sizeof(float));
  float maxval, minval ;
  int i, j, k;
  float toler = 0.0f, maxerror = 0.0f ;
  int precision = 0 ;
  void *Stream1, *Stream2, *Stream3, *Stream4;
  int Ssize, Ssize2, errors ;
  uint64_t t ;
  int d[4], s[4], nsize ;
  double t_ns ;
  float epsilon = .00000111f ;
  float noise = 0.00f ;
  float alpha = 1.0f ;
  float delta = .002f ;
  int kind = 8 ;
  IntPair tq ;
  FloatPair tf ;
  int nbits ;
  float quantum, floatbits ;
  int32_t vbits[64] ;

  if(argc > 2) sscanf(argv[2],"%f",&alpha);       // get multiplier
  if(argc > 3) sscanf(argv[3],"%f",&delta);       // get offset
  if(argc > 4) sscanf(argv[4],"%f",&noise);       // get noise level
  if(argc > 5) sscanf(argv[5],"%d",&kind);        // get kind of data fill
  printf("\n===== alpha = %f, delta = %f, noise = %f, kind = %d =====\n", alpha, delta, noise, kind) ;

  ZfpCompress_set_debug(1) ;
  /* initialize array to be compressed */
  minval = 1.0E37 ; maxval = -minval ;
  FloatTestData_2d(array, nx, nx, ny, alpha, delta, noise, kind, &minval, &maxval) ; // sum of squares

  if(argc > 1) sscanf(argv[1],"%f",&toler);                                          // get precision control
  if(toler > 0) maxerror = toler ;
  else          precision = -toler ;
  fprintf(stderr,"zfp_codec_version = %u, ZFP_CODEC = %u, codec is %s\n", 
          get_zfp_codec_version(), get_ZFP_CODEC(), zfp_codec_consistent() ? "coherent" : "DIFFERENT") ;
  fprintf(stderr,"precision = %d, max tolerance = %f, min = %f, max = %f, noise = +-%f\n", 
          precision, maxerror, minval, maxval, noise) ;
  fprintf(stderr,"\n") ;

  fprintf(stderr,"============= zfp compress/decompress 2D array =============\n");
  t  = elapsed_cycles() ;
  Stream1 = ZfpCompress(array, nx,   ny,  nz, maxerror, precision , &Ssize) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  nsize = ZfpArrayDims(Stream1, d, s, Ssize) ;
  memset(brray, 0, sizeof(float)*nx*ny*nz) ;
  t  = elapsed_cycles() ;
  ZfpExpand(brray, nx, ny, nz, Stream1, Ssize) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  AnalyzeCompressionErrors(array, brray, nx*ny*nz, epsilon, "(nx,ny,1)");     // evaluate errors
  fprintf(stderr,"\n") ;
#if 0
#if 0
  int i0, ii0, j0, jj0, ii, jj, ix, iix ;
  fprintf(stderr,"========== zfp compress/decompress 2D as (4*nx/8,4*ny/8,4) array (4x4 blocks)============\n");
  // Array is array disguised as a 3D array, 4 levels
  for(j0=0, jj0=0 ; j0<ny-7 ; j0+=8, jj0+=4){
    for(i0=0, ii0 = 0 ; i0<nx-7 ; i0+=8, ii0+=4){
      for(j=0 ; j<4 ; j++){
        for(i=0 ; i<4 ; i++){
          Array[index3(ii0+i, jj0+j, 0, nx/2, ny/2, 4)] = array[index2(i0+0+i, j0+0+j, nx, ny)] ;
          Array[index3(ii0+i, jj0+j, 1, nx/2, ny/2, 4)] = array[index2(i0+4+i, j0+0+j, nx, ny)] ;
          Array[index3(ii0+i, jj0+j, 2, nx/2, ny/2, 4)] = array[index2(i0+0+i, j0+4+j, nx, ny)] ;
          Array[index3(ii0+i, jj0+j, 3, nx/2, ny/2, 4)] = array[index2(i0+4+i, j0+4+j, nx, ny)] ;
        }
      }
    }
  }
  if(precision >0) precision += 1 ;
  t  = elapsed_cycles() ;
  Stream2 = ZfpCompress(Array, nx2, ny2, nz2, maxerror, precision , &Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  nsize = ZfpArrayDims(Stream2, d, s, Ssize) ;
  memset(brray, 0, sizeof(float)*nx*ny*nz) ;
  t  = elapsed_cycles() ;
  ZfpExpand(brray, nx2, ny2, nz2, Stream2, Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  AnalyzeCompressionErrors(Array, brray, nx*ny*nz, epsilon, "(nx/2,ny/2,4)");     // evaluate errors
  fprintf(stderr,"\n") ;
#endif
  fprintf(stderr,"========== zfp compress/decompress 2D as [1x1](nx/2,ny/2,4) array (1x1 blocks)============\n");
  // Array is array disguised as a 3D array, 4 levels
  // disguise array as a 3D (nx/2, ny/2, 4) array
  for (j = 0; j < ny2; j++){
    for (i = 0; i < nx2; i++){
      Array[i + j*nx2 + 0*nxy2] = array[2*i + 2*j*nx + 1] ;    // odd i, even j
      Array[i + j*nx2 + 1*nxy2] = array[2*i + 2*j*nx + 0] ;    // even i, even j
      Array[i + j*nx2 + 2*nxy2] = array[2*i + 2*j*nx + nx] ;   // even i, odd j
      Array[i + j*nx2 + 3*nxy2] = array[2*i + 2*j*nx + nx+1] ; // odd i, odd j
    }
  }
  // Array is array disguised as a 3D array, 4 levels
  // level 0  odd i, even j elements
  // level 1 even i, even j elements
  // level 2 even i,  odd j elements
  // level 3  odd i,  odd j elements
  if(precision >0) precision += 1 ;
  t  = elapsed_cycles() ;
  Stream2 = ZfpCompress(Array, nx2, ny2, nz2, maxerror, precision , &Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  nsize = ZfpArrayDims(Stream2, d, s, Ssize) ;
  memset(brray, 0, sizeof(float)*nx*ny*nz) ;
  t  = elapsed_cycles() ;
  ZfpExpand(brray, nx2, ny2, nz2, Stream2, Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  AnalyzeCompressionErrors(Array, brray, nx*ny*nz, epsilon, "(nx/2,ny/2,4)");     // evaluate errors
  fprintf(stderr,"\n") ;
#endif
  // array disguised as a 3D array, 4 levels
  fprintf(stderr,"========== compress/decompress 2D as (nx,4,ny/4) array ============\n");
  if(precision >0) precision += 1 ;
  t  = elapsed_cycles() ;
  Stream3 = ZfpCompress(array, nx, 4, ny/4, maxerror, precision , &Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  nsize = ZfpArrayDims(Stream3, d, s, Ssize) ;
  memset(brray, 0, sizeof(float)*nx*ny*nz) ;
  t  = elapsed_cycles() ;
  ZfpExpand(brray, nx, 4, ny/4, Stream3, Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  AnalyzeCompressionErrors(array, brray, nx*ny*nz, epsilon, "(nx,4,ny/4)");     // evaluate errors
  fprintf(stderr,"\n") ;
if(toler <0) {
  // array disguised as a 3D array, 4 levels
  fprintf(stderr,"========== compress/decompress 2D as (4,nx/4,ny) array ============\n");
  if(precision >0) precision += 1 ;
  t  = elapsed_cycles() ;
  Stream4 = ZfpCompress(array, 4, nx/4, ny, maxerror, precision , &Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  nsize = ZfpArrayDims(Stream4, d, s, Ssize) ;
  memset(brray, 0, sizeof(float)*nx*ny*nz) ;
  t  = elapsed_cycles() ;
  ZfpExpand(brray, 4, nx/4, ny, Stream4, Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  AnalyzeCompressionErrors(array, brray, nx*ny*nz, epsilon, "(4,nx/4,ny)");     // evaluate errors
  fprintf(stderr,"\n") ;
}
  if(toler < 0) return 0 ;  // bit planes mode
  fprintf(stderr,"========== lin+lorenzo compress/decompress 2D array ============\n");
//   float_quantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, IntPair *t);
  quantum = maxerror ;
  nbits = float_quantize_simple(array, qarray, nx, nx, nx, ny, quantum, &tq);
  fprintf(stderr,"quantized bits = %d, min = %d, max= %d\n", nbits, tq.t[0], tq.t[1]) ;
//   LorenzoPredict_c(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);
  LorenzoPredict_c(qarray, larray, nx, nx, nx, ny);
//   float BitEntropy(int32_t *bitstream, int npts, int nbits, int rshift)
  floatbits = BitEntropy(qarray, nx*ny, 18, 0) ;
  fprintf(stderr,"entropy per raw point = %f,", floatbits);
  floatbits = BitEntropy(larray, nx*ny, 18, 0) ;
  fprintf(stderr," entropy per predicted point = %f\n", floatbits);
//   int32_t vBitsNeeded_32(int32_t * restrict what, int32_t * restrict bits, int nn)
  for(i=0 ; i<16 ; i++) vbits[i] = 0 ;
  nbits = vBitsNeeded_32(larray, vbits, nx*ny) ;
  fprintf(stderr,"nbits = %d, ", nbits) ;
  for(i=0 ; i<16 ; i++) fprintf(stderr,"%7d", vbits[i]) ; fprintf(stderr,"\n");
//   LorenzoUnpredict_c(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);
  LorenzoUnpredict_c(qarray, larray, nx, nx, nx, ny);
//   float_unquantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, FloatPair *t);
  float_unquantize_simple(brray, qarray, nx, nx, nx, ny, quantum, &tf);
  AnalyzeCompressionErrors(array, brray, nx*ny*nz, epsilon, "lin+lorenzo");
}
