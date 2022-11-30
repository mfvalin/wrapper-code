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
#include <string.h>
#include <math.h>
#include <misc_timers.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <misc_zfp.h>
#include <misc_analyze.h>

int main(int argc, char **argv){
  /* allocate array of floats */
  int nx = 2048; int nx2 = nx/2 ;
  int ny = 1024; int ny2 = ny/2 ;
  int nz = 1; int nz2 = 4 ;
  int nxy2 = nx2 * ny2;
  float* array = malloc(nx * ny * nz * sizeof(float));
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

  if(argc > 2) sscanf(argv[2],"%f",&alpha);       // get noise level
  if(argc > 3) sscanf(argv[3],"%f",&delta);       // get noise level
  if(argc > 4) sscanf(argv[4],"%f",&noise);       // get noise level
  if(argc > 5) sscanf(argv[5],"%d",&kind);        // get kind of data fill
  printf("\n===== alpha = %f, delta = %f, noise = %f, kind = %d =====\n", alpha, delta, noise, kind) ;

  ZfpCompress_set_debug(1) ;
  /* initialize array to be compressed */
  minval = 1.0E37 ; maxval = -minval ;
  FloatTestData_2d(array, nx, nx, ny, alpha, delta, noise, kind, &minval, &maxval) ; // sum of squares

  // disguise array as a 3D (nx/2, ny/2, 4) array
  for (j = 0; j < ny2; j++){
    for (i = 0; i < nx2; i++){
      Array[i + j*nx2 + 0*nxy2] = array[2*i + 2*j*nx + 1] ;
      Array[i + j*nx2 + 1*nxy2] = array[2*i + 2*j*nx + 0] ;
      Array[i + j*nx2 + 2*nxy2] = array[2*i + 2*j*nx + nx] ;
      Array[i + j*nx2 + 3*nxy2] = array[2*i + 2*j*nx + nx+1] ;
    }
  }

  if(argc > 1) sscanf(argv[1],"%f",&toler);                                          // get precision control
  if(toler > 0) maxerror = toler ;
  else          precision = -toler ;
  fprintf(stderr,"zfp_codec_version = %u, ZFP_CODEC = %u, codec is %s\n", 
          get_zfp_codec_version(), get_ZFP_CODEC(), zfp_codec_consistent() ? "coherent" : "DIFFERENT") ;
  fprintf(stderr,"precision = %d, max tolerance = %f, min = %f, max = %f, noise = +-%f\n", 
          precision, maxerror, minval, maxval, noise) ;
  fprintf(stderr,"\n") ;

  fprintf(stderr,"============= compress/decompress 2D array =============\n");
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
  // Array is array disguised as a 3D array, 4 levels
  // level 0  odd i, even j elements
  // level 1 even i, even j elements
  // level 2 even i,  odd j elements
  // level 3  odd i,  odd j elements
  fprintf(stderr,"========== compress/decompress 2D as (nx/2,ny/2,4) array ============\n");
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
#if 0
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
#endif
}
