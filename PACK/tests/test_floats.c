//
// Copyright (C) 2022  Environnement Canada
//
// This is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// Author:
//     M. Valin,   Recherche en Prevision Numerique, august 2022
//

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <rmn/misc_types.h>
#include <rmn/misc_timers.h>
#include <misc_floats.h>
// int verify_upper(uint32_t *a, uint32_t *b, int n);
// int verify_lower(uint32_t *a, uint32_t *b, int n);
int verify_masked(void *a, void *b, uint32_t mask, int n);
void apply_round(void *a, void *b, uint32_t round, int n);

// #define NPTS 1031
#define NPTS 4103
#define NTIMES 10000

int main(int argc, char **argv){
  int i, j, nt, errors,  n = 7 ;
  static uint32_t fp32[8] = { 0x00810203, 0x04850607, 0x08890A0B, 0x0C8D0E0F, \
                              0x10911213, 0x14951617, 0x18991A1B, 0x1C9D1E1F  } ;
  static uint32_t uf24[8]  ;
  static uint32_t fp32b [8] ;
  uint32_t lfp32[NPTS+1], luf24[NPTS+1], lfp32b[NPTS+1] ;
  float float32[NPTS], float32b[NPTS], float32m[NPTS] ;
  uint32_t *ifloat32 = (uint32_t *) float32 ;
  uint32_t *ifloat32b = (uint32_t *) float32b ;
  uint32_t *ifloat32m = (uint32_t *) float32m ;
  uint64_t t0, t1, tmin, tmax, tavg, freq ;
  double nano ;
  uint64_t t[NTIMES] ;
  double bias, avg ;
  int nexpbits ;
  FloatInt fmax ;

  freq = cycles_counter_freq() ;
  nano = 1000000000 ;
  nano /= freq ;
  for(i=0 ; i<NPTS ; i++) { lfp32[i] = i ; luf24[i] = 0 ; lfp32b[i] = 0 ; }
//   for(i=0 ; i<NPTS ; i++) { float32[i] = NPTS + (i - .111)/NPTS ; float32b[i] = -1 ; }
//   for(i=0 ; i<NPTS ; i++) { float32[i] = (i + .111)/NPTS ; float32b[i] = -1 ; }
  for(i=0 ; i<NPTS ; i++) { float32[i] = (i + .111)/NPTS*40.0f ; float32b[i] = -1 ; }

// ========================== test "brainfloat 24" ==========================
  memset(uf24, 0xFF, sizeof(uf24));
  fp32_bf24((void *) fp32, uf24, n) ;
  printf("fp32      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;

  printf("uf24      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;

  memset(fp32b, 0xFF, sizeof(fp32b));
  bf24_fp32(fp32b, uf24, n) ;
//   errors = verify_upper((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  errors = verify_masked((uint32_t *) fp32, (uint32_t *) fp32b, 0xFFFFFF00u, n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;

  printf("fp32>bf24 : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    fp32_bf24((void *) float32, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("bf24>fp32 : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    bf24_fp32(float32b, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_upper((uint32_t *) float32, (uint32_t *) float32b, NPTS) ;
  apply_round(float32, float32m, 0x80, NPTS);
  errors = verify_masked((uint32_t *) float32m, (uint32_t *) float32b, 0xFFFFFF00u,  NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;
  bias = 0.0 ; avg = 0.0 ;
  for(i=0 ; i<NPTS ; i++) {bias += (float32[i] - float32b[i]) ; avg += float32[i] ; }
  bias /= NPTS ; avg /= NPTS ;
  printf("sample float : %f %f, relative error = %8.3g, avg = %f, bias = %8.3g\n", 
         float32[NPTS/2], float32b[NPTS/2], (float32[NPTS/2]-float32b[NPTS/2])/float32[NPTS/2], avg, bias) ;
  printf("===============================================================\n") ;

// ========================== test "brainfloat 16" ==========================
  memset(uf24, 0xFF, sizeof(uf24));
  fp32_bf16((void *) fp32, uf24, n) ;
  printf("fp32      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;

  printf("uf16      : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;

  memset(fp32b, 0xFF, sizeof(fp32b));
  bf16_fp32(fp32b, uf24, n) ;
//   errors = verify_upper((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  errors = verify_masked((uint32_t *) fp32, (uint32_t *) fp32b, 0xFFFF0000u, n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;

  printf("fp32      : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    fp32_bf16((void *) float32, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  apply_round(float32, float32m, 0x8000, NPTS);
  for(i=0 ; i<16 ; i++) printf("%8.8x ",ifloat32m[i]) ; printf("\n") ;
  printf("luf16     : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x          ",luf24[i])     ; printf("\n") ;
//   fp32_nf16((void *) float32, luf24, NPTS, 8) ;
//   printf("luf16b    : ");
//   for(i=0 ; i<8 ; i++) printf("%8.8x          ",luf24[i])     ; printf("\n") ;
  printf("fp32>bf16 : ");
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("bf16>fp32 : ");
  bf16_fp32(float32b, luf24, NPTS) ;
  for(i=0 ; i<16 ; i++) printf("%8.8x ",ifloat32b[i]) ; printf("\n") ;
  printf("bf16>fp32 : ");
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    bf16_fp32(float32b, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_upper((uint32_t *) float32, (uint32_t *) float32b, NPTS) ;
//   for(i=0 ; i<16 ; i++) printf("%8.8x ",ifloat32b[i]) ; printf("\n") ;
  apply_round(float32, float32m, 0x8000, NPTS);
  errors = verify_masked((uint32_t *) float32m, (uint32_t *) float32b, 0xFFFF0000u,  NPTS) ;
//   printf(">fp32     : ");
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;
  bias = 0.0 ; avg = 0.0 ;
  for(i=0 ; i<NPTS ; i++) {bias += (float32[i] - float32b[i]) ; avg += float32[i] ; }
  bias /= NPTS ; avg /= NPTS ;
  printf("sample float : %f %f, relative error = %8.3e, avg = %f, bias = %8.3g\n", 
         float32[NPTS/2], float32b[NPTS/2], (float32[NPTS/2]-float32b[NPTS/2])/float32[NPTS/2], avg, bias) ;

// ========================== test float 16 variants ==========================
  for(nexpbits = 8 ; nexpbits > 3 ; nexpbits--){
  printf("===============================================================\n") ;
  printf(">nf16-%d   : ", nexpbits);
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    fp32_nf16((void *) float32, luf24, NPTS, nexpbits) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  apply_round(float32, float32m, 0x8000, NPTS);
  for(i=0 ; i<16 ; i++) printf("%8.8x ",ifloat32m[i]) ; printf("\n") ;
  printf("luf16     : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x          ",luf24[i])     ; printf("\n") ;
//   fp32_nf16((void *) float32, luf24, NPTS, nexpbits) ;
//   printf("luf16b    : ");
//   for(i=0 ; i<8 ; i++) printf("%8.8x          ",luf24[i])     ; printf("\n") ;
  printf("fp32>f16-%d: ", nexpbits);
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("f16-%d>fp32: ", nexpbits);
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    nf16_fp32(float32b, luf24, NPTS, nexpbits) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
    tmax = (t[j] > tmax) ? t[j] : tmax ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.5*tmin) {
      tavg = tavg + t[j] ;
      nt++ ;
    }
  }
//   errors = verify_upper((uint32_t *) float32, (uint32_t *) float32b, NPTS) ;
  apply_round(float32, float32m, 0x8000, NPTS);
  errors = verify_masked((uint32_t *) float32m, (uint32_t *) float32b, 0xFFFF0000u,  NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;
  bias = 0.0 ; avg = 0.0 ;
  for(i=0 ; i<NPTS ; i++) {bias += (float32[i] - float32b[i]) ; avg += float32[i] ; }
  bias /= NPTS ; avg /= NPTS ;
  fmax.f = nf16_floatmax(nexpbits) ;
  printf("sample float : %12.7f %12.7f, relative error = %8.3e, avg = %f, bias = %8.3g\n", 
         float32[NPTS/2], float32b[NPTS/2], (float32[NPTS/2]-float32b[NPTS/2])/float32[NPTS/2], avg, bias) ;
  printf("sample float : %12.7f %12.7f, relative error = %8.3e, maxfloat = %11.7e (%8.8x)\n", 
         float32[NPTS-1], float32b[NPTS-1], (float32[NPTS-1]-float32b[NPTS-1])/float32[NPTS-1], fmax.f, fmax.u) ;
  }
return 0 ;
// ========================== test unsigned 24 bit integers ==========================
  printf("===============================================================\n") ;
  memset(uf24, 0xFF, sizeof(uf24));
  u32_u24(fp32, uf24, n) ;
  printf("u32       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    u32_u24(lfp32, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("u24       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;
  memset(fp32b, 0xFF, sizeof(fp32b));
  u24u32(fp32b, uf24, n) ;
//   errors = verify_lower((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  errors = verify_masked((uint32_t *) fp32, (uint32_t *) fp32b, 0xFFFFFF, n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    u24u32(lfp32b, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_lower((uint32_t *) lfp32, (uint32_t *) lfp32b, NPTS) ;
  errors = verify_masked((uint32_t *) lfp32, (uint32_t *) lfp32b, 0xFFFFFF, NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;

// ========================== test signed 24 bit integers ==========================
  printf("===============================================================\n") ;

  memset(uf24, 0xFF, sizeof(uf24));
  u32_u24(fp32, uf24, n) ;
  printf("i32       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32[i]) ; printf("\n") ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    u32_u24(lfp32, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano) ;

  printf("u24       : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",uf24[i]) ; printf("\n") ;
  memset(fp32b, 0xFF, sizeof(fp32b));
  u24i32((int32_t *) fp32b, uf24, n) ;
//   errors = verify_lower((uint32_t *) fp32, (uint32_t *) fp32b, n) ;
  errors = verify_masked((uint32_t *) fp32, (uint32_t *) fp32b, 0xFFFFFF,  n) ;
  printf("restored  : ");
  for(i=0 ; i<8 ; i++) printf("%8.8x ",fp32b[i]) ; printf(" errors = %d\n", errors) ;
  tmin = 1000000000 ;  nt = 0 ; tavg = 0.0f ; tmax = 0 ;
  for(j=0 ; j < NTIMES ; j++) {
    t0 = elapsed_cycles() ;
    u24i32((int32_t *) lfp32b, luf24, NPTS) ;
    t1 = elapsed_cycles() ;
    t[j] = t1 - t0 ;
    tmin = (t[j] < tmin) ? t[j] : tmin ;
  }
  for(j=0 ; j < NTIMES ; j++) {
    if(t[j] < 1.25*tmin) {
      tavg = tavg + t[j] ;
      tmax = (t[j] > tmax) ? t[j] : tmax ;
      nt++ ;
    }
  }
//   errors = verify_lower((uint32_t *) lfp32, (uint32_t *) lfp32b, NPTS) ;
  errors = verify_masked((uint32_t *) lfp32, (uint32_t *) lfp32b, 0xFFFFFF,  NPTS) ;
  printf("NPTS = %d, ns = %6.1f ->%8.1f [avg = %6.1f, %6.3f ns/pt] (%6d) %f, errors = %d\n",
         NPTS, tmin*nano, tmax*nano, tavg/nt*nano, tavg/nt*nano/NPTS, nt, nano, errors) ;

  printf("\n\n\n\n\n");
}
