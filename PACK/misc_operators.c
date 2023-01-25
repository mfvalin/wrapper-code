// Hopefully useful code for C
// Copyright (C) 2022  Recherche en Prevision Numerique
//
// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//

#if ! defined(__INTEL_COMPILER_UPDATE)
#pragma GCC optimize "tree-vectorize"
#endif

#define STATIC extern
#include <rmn/misc_operators.h>

#include <math.h>
#include <stdlib.h>

// get the IEEE exponent of the largest float (absolute value)
// and the smallest non zero float (absolute value)  in float array f
// if all values of f are 0.0, 255 will be returned for the minimum exponent
uint32_t ieee_minmax_exponent(float *f, int n){
  int i ;
  uint32_t *uf = (uint32_t *) f ;
  uint32_t t, emax = 0 , emin = 0xFF << 24 ;

  for(i=0 ; i<n ; i++){
    t = uf[i] << 1 ;                  // get rid of sign, exponent in upper 8 bits
    emax = (t > emax) ? t : emax ;
    t = (t == 0) ? emin : t ;         // ignore values of 0
    emin = (t < emin) ? t : emin ;
  }
  emax >>= 24 ;                       // exponent WITH IEEE bias (127)
  emin >>= 24 ;                       // exponent WITH IEEE bias (127)
  return emax | (emin << 8) ;         // emax and emin
}

// get the IEEE exponent of the largest float (absolute value) in float array f
uint32_t ieee_max_exponent(float *f, int n){
  int i ;
  uint32_t *uf = (uint32_t *) f ;
  uint32_t t, emax = 0 ;

  for(i=0 ; i<n ; i++){
    t = uf[i] << 1 ;                  // get rid of sign, exponent in upper 8 bits
    emax = (t > emax) ? t : emax ;
  }
  emax >>= 24 ;                       // exponent WITH IEEE bias (127)
  return emax ;
}

// get the IEEE exponent of the smallest non zero float (absolute value) in float array f
// if all values of f are 0.0, 255 will be returned
uint32_t ieee_min_exponent(float *f, int n){
  int i ;
  uint32_t *uf = (uint32_t *) f ;
  uint32_t t, emin = 0xFF << 24 ;     // largest possible exponent

  for(i=0 ; i<n ; i++){
    t = uf[i] << 1 ;                  // get rid of sign, exponent in upper 8 bits
    t = (t == 0) ? emin : t ;         // ignore values of 0
    emin = (t < emin) ? t : emin ;
  }
  emin >>= 24 ;                       // exponent WITH IEEE bias (127)
  return emin ;
}

int32_t vBitsNeeded_u32(uint32_t * restrict what, int32_t * restrict bits, int nn){
  int i, needed=0 ;
  int n = (nn < 0) ? -nn : nn ;
  for(i=0 ; i<n ; i++) {
    int32_t nbits = BitsNeeded_u32(what[i]) ;
    if(nn >= 0) bits[nbits]++ ;
    needed = (nbits > needed) ? nbits : needed ;
  }
  return needed ;
}

int32_t vBitsNeeded_32(int32_t * restrict what, int32_t * restrict bits, int nn){
  int i, needed=0 ;
  int n = (nn < 0) ? -nn : nn ;
  for(i=0 ; i<n ; i++) {
    int32_t nbits = BitsNeeded_u32(what[i]) ;
    if(nn >= 0) bits[nbits]++ ;
    needed = (nbits > needed) ? nbits : needed ;
  }
  return needed ;
}

int32_t vBitsNeeded_u64(uint64_t * restrict what, int32_t * restrict bits, int nn){
  int i, needed=0 ;
  int n = (nn < 0) ? -nn : nn ;
  for(i=0 ; i<n ; i++) {
    int32_t nbits = BitsNeeded_u64(what[i]) ;
    if(nn >= 0) bits[nbits]++ ;
    needed = (nbits > needed) ? nbits : needed ;
  }
  return needed ;
}

int32_t vBitsNeeded_64(int64_t * restrict what, int32_t * restrict bits, int nn){
  int i, needed=0 ;
  int n = (nn < 0) ? -nn : nn ;
  for(i=0 ; i<n ; i++) {
    int32_t nbits = BitsNeeded_64(what[i]) ;
    if(nn >= 0) bits[nbits]++ ;
    needed = (nbits > needed) ? nbits : needed ;
  }
  return needed ;
}

void BitEntropy4(float entropy[4], uint32_t *bitstream, int npts, int nbits, int rshift) {
  uint32_t bins[4][16] ;
  uint32_t mask = 0xF ;
  int i, j ;
  uint32_t t ;
  float e[4] ;
  float prob, scale, ovlog2 ;

  scale = 1.0f / (float)(npts) ;
  ovlog2 = 1.0f / logf(2.0f) ;
  for(j=0 ; j<4 ; j++) for(i=0 ; i<16 ; i++) bins[j][i] = 0 ;
  for(j=0 ; j<4 ; j++) e[j] = 0.0f ;
  for (i = 0; i < npts; i++) {
    t = bitstream[i] ;
    bins[0][t & mask]++ ; t >>= 4 ;
    bins[1][t & mask]++ ; t >>= 4 ;
    bins[2][t & mask]++ ; t >>= 4 ;
    bins[3][t & mask]++ ;
  }
  for(j=0 ; j<4 ; j++){
    for(i=0 ; i<16 ; i++){
      if (bins[j][i] != 0) {
        prob = (float)bins[j][i] * scale ;
        e[j] += (prob * logf(prob) * ovlog2) ;
      }
    }
  }
  for(j=0 ; j<4 ; j++) entropy[j] = (-e[j]) ;
}

float BitEntropy(int32_t *bitstream, int npts, int nbits, int rshift) {
    int lbin, sizebins;
//     int imin, imax, range, nbits_local ;
    float prob, entropie, scale, ovlog2 ;
    uint32_t *bins;
    uint32_t mask ;
    uint32_t samples = 0 ;

    entropie = 0.0f;
    mask = RMASK32(nbits) ;
    // printf("imin : %d imax : %d range: %d, nbits : %d\n", imin, imax, imax-imin, nbits_local);

    sizebins = 1 << (nbits-rshift) ;
    bins = (uint32_t *) calloc(sizebins,sizeof(uint32_t));

    for (int i = 0; i < npts; i++) {
        lbin = (bitstream[i] & mask) >> rshift ;
        bins[lbin]++;
    }

    scale = 1.0f / (float)(npts) ;
    ovlog2 = 1.0f / logf(2.0f) ;
    for (int i = 0; i < sizebins; i++) {
        if (bins[i] != 0) {
            prob = (float)(bins[i]) * scale ;
            // printf("i: %d count: %d prob: %f contrib : %f \n", i, bins[i], prob, (prob * log(prob)/log(2.0)));
            entropie += prob * logf(prob) * ovlog2 ;
            samples++ ;
        }
    }

    free(bins) ;
// printf("sizebins = %d, mask = %8.8x, rshift = %d, samples = %d, entropie = %g\n", sizebins, mask, rshift, samples, entropie);
    return -entropie ;
}

