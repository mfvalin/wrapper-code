/*
 * Hopefully useful code for C and Fortran
 * Copyright (C) 2022  Recherche en Prevision Numerique
 *
 * This code is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 * 
 * Author : M.Valin 2022/09
 */
#include <stdint.h>

#include <with_simd.h>

// in place 1-2-1 / 4 smoothing
extern inline void Ismooth_124_inplace(int32_t *f, int n){
  int32_t a, b, c ;
  int i ;
  b = f[0] ; c = f[1] ;
//   f[0] = (b + c) / 2 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
#else
#endif
  for(i=1 ; i<n-1 ; i++){
    a = b ; b = c ; c = f[i+1] ;
    f[i] = (a + 2*b + c)/4 ;
  }
//   f[n-1] = (b + c) / 2 ;
}

// in place 1-2-1 / 4 smoothing
extern inline void Fsmooth_124_inplace(float *f, int n){
  float a, b, c ;
  int i ;
  b = f[0] ; c = f[1] ;
  for(i=1 ; i<n-1 ; i++){        // first and last element left untouched
    a = b ; b = c ; c = f[i+1] ;
    f[i] = (a + 2.0f*b + c) * .25f ;
  }
}

// 1-2-1 / 4 smoothing
extern inline void Fsmooth_124(float *s, float *d, int n){
  int i ;
  d[0] = s[0] ;            // first element left untouched
  for(i=1 ; i<n-1 ; i++){
    d[i] = (s[i-1] + 2.0f*s[i] + s[i+1]) * .25f ;
  }
  d[n-1] = s[n-1] ;        // last element left untouched
}

// in place 2 dimensional 1-2-1 / 4 smoothing
void Fsmooth_124_2D_inplace(float *f, int ni, int lni, int nj){
  float tr[ni], t ;
  int i, j ;

  Fsmooth_124_inplace(f, ni) ;                   // smooth row 0
  for(i=0 ; i<ni ; i++) tr[i] = f[i] ;
  Fsmooth_124_inplace(f+lni, ni) ;               // smooth row 1
  for(j=1 ; j<nj-1 ; j++){                       // only 1D smoothing for first and last row
    f += lni ;                                   // point to row j
    Fsmooth_124_inplace(f+lni, ni) ;             // smooth row j+1
    for(i=0 ; i<ni ; i++) {
      t = (tr[i] + 2.0f*f[i] + f[i+lni]) * .25f ;     // new value for row j
      tr[i] = f[i] ;                             // save row j for next pass
      f[i] = t ;                                 // store updated row j
    }
  }
}

// 2 dimensional 1-2-1 / 4 smoothing
void Fsmooth_124_2D(float *s, float *d, int ni, int lnis, int lnid, int nj){
  float tr[ni], t ;
  int i, j ;

  Fsmooth_124(s,      d,      ni) ;              // smoth row 0 along i
  for(i=0 ; i<ni ; i++) tr[i] = d[i] ;           // save row 0
  Fsmooth_124(s+lnis, d+lnid, ni) ;              // smoth row 1 along i
  for(j=1 ; j<nj-1 ; j++){
    s += lnis ;                                  // point to row j
    d += lnid ;                                  // point to row j
    Fsmooth_124(s+lnis, d+lnid, ni) ;            // smoth row j+1 along i
    for(i=0 ; i<ni ; i++) {
      t = (tr[i] + 2.0f*d[i] + d[i+lnid]) * .25f ;     // new value for row j
      tr[i] = d[i] ;                             // save row j for next pass
      d[i] = t ;                                 // store updated row j
    }
  }
}
