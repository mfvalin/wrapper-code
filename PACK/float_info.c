//  LIBRMN - Library of hopefully useful routines for C and FORTRAN
//  Copyright (C) 2022  Recherche en Prevision Numerique
//                      Environnement Canada
// 
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation,
//  version 2.1 of the License.
// 
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Library General Public License for more details.
//

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <float.h>

#if defined(__GNUC__)
// #if defined(__INTEL_COMPILER) || defined(__clang__) || defined(__PGIC__)
// #undef __GCC_IS_COMPILER__
// #else
#define __GCC_IS_COMPILER__
// #endif
#endif

// use Intel compatible intrinsics for SIMD instructions
#include <with_simd.h>

#define VL 8
#define VLMASK (VL - 1)

// float_info_simple is mostly used to check C compiler vectorization capabilities
int float_info_simple(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs){
  int i, j ;
  float zmax, zmin, abs, amin ; ;

  zmax = zmin = zz[0] ;
  amin = FLT_MAX ;
  if(lni == ni) { ni = ni * nj ; nj = 1 ; }

  for(j = 0 ; j < nj ; j++){                       // loop over rows
    for(i = 0 ; i < ni ; i++){                     // loop for a row
      zmax = (zz[i] > zmax) ? zz[i]  : zmax ;      // highest value
      zmin = (zz[i] < zmin) ? zz[i]  : zmin ;      // lowest value
      abs  = (zz[i] < 0.0f) ? -zz[i] : zz[i] ;     // absolute value
      abs  = (abs == 0.0f)  ? amin   : abs ;       // set to amin if zero
      amin = (abs < amin)   ? abs    : amin ;      // smallest NON ZERO absolute value
    }
    zz += lni ;
  }
  *maxval = zmax ;
  *minval = zmin ;
  *minabs = amin ;
  return 0 ;
}

// Fortran layout(lni,nj) is assumed, i index varying first
// zz     [IN]  : 32 bit floating point array
// ni     [IN]  : number of useful points along i
// lni    [IN]  : actual row length ( >= ni )
// nj     [IN]  : number of rows (along j)
// maxval [OUT] : highest value in zz
// minval [OUT] : lowest value in zz
// minabs [OUT] : smallest NON ZERO absolute value in zz
// return value : 0 (to be consistent with float_info_missing
int float_info_no_missing(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs){
  int i0, i, j ;
  float z[VL] ;
  float zmax[VL], zma ;
  float zmin[VL], zmi ;
  float abs[VL], zabs ;
  float amin[VL], ami ;
  int nfold, incr;

  if(lni == ni) { ni = ni * nj ; nj = 1 ; }

  if(ni < VL) {        // less than VL elements along i
    zma = zz[0] ;
    zmi = zma ;
    ami = FLT_MAX ;
    for(j = 0 ; j < nj ; j++){
      for(i = 0 ; i < ni ; i++){
        zma   = (zz[i] > zma)  ? zz[i]  : zma ;
        zmi   = (zz[i] < zmi)  ? zz[i]  : zmi ;
        zabs  = (zz[i] < 0.0f) ? -zz[i] : zz[i] ;
        zabs  = (zabs == 0.0f) ? ami    : zabs ;
        ami   = (zabs < ami)   ? zabs   : ami ;
      }
      zz += lni ;
    }
    *minabs = ami ;    // smallest NON ZERO absolute value
    *maxval = zma ;    // highest value
    *minval = zmi ;    // lowest value

  }else{
    for(i=0 ; i<VL ; i++){
      zmin[i] = FLT_MAX ;
      zmax[i] = -zmin[i] ;
      amin[i] = FLT_MAX ;
    }
    nfold = ni & VLMASK ;                                       // ni modulo VL

    for(j = 0 ; j < nj ; j++){                                  // loop over rows
      incr = (nfold > 0) ? nfold : VL ;                         // first increment is shorter if ni not a multiple of VL
      // there may be some overlap between pass 0 and pass 1
      // elements nfold -> VL -1 will be processed twice
      // but it does not matter  for min/max/abs/absmin
      // one SIMD pass is cheaper than a scalar loop
      for(i0 = 0 ; i0 < ni-VL+1 ; i0 += incr){                  // loop for a row
        for(i=0 ; i<VL ; i++){                                  // blocks of VL values
          z[i]      = zz[i0+i] ;
          zmax[i]   = (z[i] > zmax[i]) ? z[i]  : zmax[i] ;      // highest value of z
          zmin[i]   = (z[i] < zmin[i]) ? z[i]  : zmin[i] ;      // lowest value of z 
          abs[i]    = (z[i] < 0.0f)      ? -z[i]   : z[i] ;     // absolute value of z
          abs[i]    = (abs[i] == 0.0f)   ? amin[i] : abs[i] ;   // if zero, set absolute value to current minimum
          amin[i]   = (abs[i] < amin[i]) ? abs[i]  : amin[i] ;  // smallest NON ZERO absolute value
        }
        incr = VL ;
      }
      zz += lni ;
    }

    for(i=1 ; i < VL ; i++){                               // fold temporary vectors into their first element
      zmin[0] = (zmin[i] < zmin[0]) ? zmin[i] : zmin[0] ;  // fold lowest value into zmin[0]
      zmax[0] = (zmax[i] > zmax[0]) ? zmax[i] : zmax[0] ;  // fold highest value into zmax[0]
      amin[0] = (amin[i] < amin[0]) ? amin[i] : amin[0] ;  // fold smallest NON ZERO absolute value into amin[0]
    }
    *minabs = amin[0] ;    // smallest NON ZERO absolute value
    *maxval = zmax[0] ;    // highest value
    *minval = zmin[0] ;    // lowest value
  }

  return 0 ;         // number of values equal to special value
}

// as some compilers generate code with poor performance, an alternative version
// using Intel instrinsics for SIMD instructions is provided
// it is activated by adding -DWITH_SIMD to the compiler options
// Fortran layout(lni,nj) is assumed, i index varying first
// zz     [IN]  : 32 bit floating point array
// ni     [IN]  : number of useful points along i
// lni    [IN]  : actual row length ( >= ni )
// nj     [IN]  : number of rows (along j)
// maxval [OUT] : highest value in zz
// minval [OUT] : lowest value in zz
// minabs [OUT] : smallest NON ZERO absolute value in zz
// spval  [IN]  : any (value & spmask) in zz that is equal to (spval & spmask) will be IGNORED
//                should spval be a NULL pointer, there is no such value
// spmask  [IN] : mask for spval
// return value : number of elements in zz that are equal to the special value
int float_info_missing(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs, float *spval, uint32_t spmask){
  int i0, i, j ;
  uint32_t *iz ;
  float sz ;
  float zmax[VL], zma ;
  float zmin[VL], zmi ;
  float zabs ;
  float amin[VL], ami ;
  int count[VL], fixcount ;
  int nfold, incr, cnt ;
  uint32_t *ispval = (uint32_t *) spval ;
  uint32_t missing ;
#if defined(__AVX2__) && defined(__x86_64__) && defined(__GCC_IS_COMPILER__) && defined(WITH_SIMD)
  __m256i vcnt, vs1, vmissing, vmask, vti ;
  __m256  v, vt, vmax, vmin, vamin, vma, vmi, vsign, vs2, vzero ;
#endif

  // if no special value, call simpler, faster function
  if(spval == NULL) return float_info_no_missing(zz, ni, lni, nj, maxval, minval, minabs) ;
  spmask = (~spmask) ;            // complement spmask

  if(lni == ni) { ni = ni * nj ; nj = 1 ; }
  missing = (*ispval) & spmask ;  // special value & mask
  fixcount = 0 ;

  if(ni < VL) {        // less than VL elements along i
    zma = -FLT_MAX ;
    zmi = FLT_MAX ;
    ami = FLT_MAX ;
    cnt = 0 ;
    for(j = 0 ; j < nj ; j++){
      iz = (uint32_t *) zz ;
      for(i = 0 ; i < ni ; i++){
        sz    = ((iz[i]  & spmask) == missing) ? zma : zz[i] ;
        zma   = (sz > zma)  ? sz  : zma ;
        sz    = ((iz[i]  & spmask) == missing) ? zmi : zz[i] ;
        zmi   = (sz < zmi)  ? sz  : zmi ;
        sz    = ((iz[i]  & spmask) == missing) ? ami : zz[i] ;
        zabs  = (sz < 0.0f)    ? -sz  : sz ;
        zabs  = (zabs == 0.0f) ? ami  : zabs ;
        ami   = (zabs < ami)   ? zabs : ami ;
        cnt  += ((iz[i]  & spmask) == missing) ? 0 : 1 ;
      }
      zz += lni ;
    }
    *minabs = ami ;    // smallest NON ZERO absolute value (FLT_MAX if all values are missing)
    *maxval = zma ;    // highest value (-FLT_MAX if all values are missing)
    *minval = zmi ;    // lowest value (FLT_MAX if all values are missing)
    count[0] = cnt ;

  }else{
    for(i=0 ; i<VL ; i++){
      zmin[i] = FLT_MAX ;
      zmax[i] = -FLT_MAX ;
      amin[i] = FLT_MAX ;
      count[i] = 0 ;
    }
    nfold = ni & VLMASK ;                                       // ni modulo VL

#if defined(__AVX2__) && defined(__x86_64__) && defined(__GCC_IS_COMPILER__) && defined(WITH_SIMD)
    vmin  = _mm256_loadu_ps(zmin) ;
    vmax  = _mm256_loadu_ps(zmax) ;
    vamin = _mm256_loadu_ps(amin) ;
    vcnt  = _mm256_loadu_si256 ((__m256i *) count) ;
    vmissing = _mm256_set1_epi32(missing) ;                     // special value & mask
    vmask = _mm256_set1_epi32(spmask) ;                         // mask for special value
    vsign = _mm256_set1_ps(-0.0f) ;
    vzero = _mm256_xor_ps(vsign, vsign) ;                       // set to 0
#endif
    for(j = 0 ; j < nj ; j++){                                  // loop over rows
      incr = (nfold > 0) ? nfold : VL ;                         // first increment is shorter if ni not a multiple of VL
      iz = (uint32_t *) zz ;
      // there may be some overlap between pass 0 and pass 1
      // elements incr -> VL -1 will be processed twice
      // but it does not matter  for min/max/abs/absmin
      // one SIMD pass is cheaper than a scalar loop
      for(i = incr ; i < VL ; i++)                              // some missing values could be counted twice
        if((iz[i]  & spmask) == missing) {                                  // scan iz[incr:VL-1] to find and count them
          fixcount++ ;
        }
      for(i0 = 0 ; i0 < ni-VL+1 ; i0 += incr){                       // loop over a row
#if defined(__AVX2__) && defined(__x86_64__) && defined(__GCC_IS_COMPILER__) && defined(WITH_SIMD)
        v = _mm256_loadu_ps(zz + i0) ;                               // get values
        vti = _mm256_and_si256((__m256i) v, vmask) ;                  // apply special value mask
        vs1 = _mm256_cmpeq_epi32(vmissing, vti) ;                     // compare with masked missing value

        vma = _mm256_blendv_ps(v, vmax, (__m256) vs1) ;              // vmax if equal to missing value
        vmi = _mm256_blendv_ps(v, vmin, (__m256) vs1) ;              // vmin if equal to missing value
        vt = _mm256_blendv_ps(v, vamin, (__m256) vs1) ;              // vamin if equal to missing value

        vmax = _mm256_max_ps(vmax, vma) ;                            // max(vma, vmax)
        vmin = _mm256_min_ps(vmin, vmi) ;                            // min(vmi, vmin)

        vt = _mm256_andnot_ps(vsign, vt) ;                           // flush sign (absolute value)
        vs2 = _mm256_cmp_ps(vt, vzero, 0) ;                          // compare for equality
        vcnt  = _mm256_sub_epi32(vcnt, vs1) ;                        // count++ where equal to missing value
        vt = _mm256_blendv_ps(vt, vamin, (__m256) vs2) ;             // vamin if 0
        vamin = _mm256_min_ps(vamin, vt) ;                           // min of absolute value
#else
        float z[VL], z1[VL], z2[VL], abs[VL] ;
        uint32_t zi[VL] ;
        for(i=0 ; i<VL ; i++){                                  // blocks of VL values
          z[i]      = zz[i0+i] ;                                // next set of values
          zi[i]     = iz[i0+i] & spmask;

          // replace zi[i] == missing with (z[i] & missing_mask) == missing
          z1[i]     = (zi[i] == missing) ? zmax[i] : z[i] ;     // if missing, set to current highest value
          zmax[i]   = (z1[i] > zmax[i])  ? z1[i]   : zmax[i] ;  // highest value of z

          z2[i]     = (zi[i] == missing) ? zmin[i] : z[i] ;     // if missing, set to current lowest value
          zmin[i]   = (z2[i] < zmin[i])  ? z2[i]   : zmin[i] ;  // lowest value of z 

          abs[i]    = (zi[i] == missing) ? amin[i] : z[i] ;     // if missing, set to current minimum
          abs[i]    = (abs[i] < 0.0f)    ? -abs[i] : abs[i] ;   // absolute value of z
          abs[i]    = (abs[i] == 0.0f)   ? amin[i] : abs[i] ;   // if zero, set absolute valkue to current minimum
          amin[i]   = (abs[i] < amin[i]) ? abs[i]  : amin[i] ;  // smallest NON ZERO absolute value

          count[i] += ( (zi[i] == missing) ? 1 : 0 ) ;
        }
#endif
        incr = VL ;
      }
      zz += lni ;
    }

#if defined(__AVX2__) && defined(__x86_64__) && defined(__GCC_IS_COMPILER__) && defined(WITH_SIMD)
    _mm256_storeu_ps(zmin, vmin) ;                         // store SIMD registers
    _mm256_storeu_ps(zmax, vmax) ;
    _mm256_storeu_ps(amin, vamin) ;
    _mm256_storeu_si256((__m256i *) count, vcnt) ;
#endif
    for(i=1 ; i < VL ; i++){                               // fold temporary vectors into their first element
      zmin[0] = (zmin[i] < zmin[0]) ? zmin[i] : zmin[0] ;  // fold lowest value into zmin[0]
      zmax[0] = (zmax[i] > zmax[0]) ? zmax[i] : zmax[0] ;  // fold highest value into zmax[0]
      amin[0] = (amin[i] < amin[0]) ? amin[i] : amin[0] ;  // fold smallest NON ZERO absolute value into amin[0]
      count[0] += count[i] ;
    }
    *minabs = amin[0] ;    // smallest NON ZERO absolute value (FLT_MAX if all values are missing)
    *maxval = zmax[0] ;    // highest value (-FLT_MAX if all values are missing)
    *minval = zmin[0] ;    // lowest value (FLT_MAX if all values are missing)
  }
  return count[0] - fixcount ;         // number of values not equal to special value
}

// Fortran layout(lni,nj) is assumed, i index varying first
// zz     [IN]  : 32 bit floating point array
// ni     [IN]  : number of useful points along i
// lni    [IN]  : actual row length ( >= ni )
// nj     [IN]  : number of rows (along j)
// maxval [OUT] : highest value in zz
// minval [OUT] : lowest value in zz
// minabs [OUT] : smallest NON ZERO absolute value in zz
// spval  [IN]  : any value in zz that is equal to spval will be IGNORED
//                should spval be a NULL pointer, there is no such value
// return value : number of elements in zz that are equal to the special value
int float_info(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs, float *spval, uint32_t spmask){
  if(spval == NULL){
    return float_info_no_missing(zz, ni, lni, nj, maxval, minval, minabs) ;
  }else{
    return float_info_missing(zz, ni, lni, nj, maxval, minval, minabs, spval, spmask) ;
  }
}

