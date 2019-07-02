/* 
 * Copyright (C) 2019  Recherche en Prevision Numerique
 *                     Environnement Canada
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 */

// discrete forward linear lifted wavelet transform (scalar, possibly in place no copy form)
void DFWT53_1d(float *f1, float *f, int n){
  int i ;
  float olo, ohi, elo, ehi ;

  elo = f[0] ;
  olo = f[1] - 0.5f * (f[2] + elo) ;         // first predicted odd term
  for(i = 0 ; i < n - 2 ; i += 2){           
    ehi = f[i+2] ;                           // upper even term
    ohi = f[i+1]  - 0.5f * (elo  + ehi ) ;   // predict odd term
    f1[i+1] = ohi ;                          // store odd term
    f1[i  ] = elo + .25f * (olo + ohi) ;     // update even term and store it
    elo = ehi ;                              
    olo = ohi ;                              
  }                                          
  if(n & 1) {                                // n is odd
    f1[n-1] = f[n-1] + 0.5f * ohi ;          // update and store last even term
  }else{                                     
    ohi = f[n-1] - ehi ;                     // predict last odd term
    f1[n-1] = ohi ;                          // store last odd term
    f1[n-2] = f[n-2] + .25f * (olo + ohi) ;  // update and store last even term
  }
}

// discrete forward linear lifted wavelet transform (scalar form)
// with split even/odd components
// number of points n is neven + nodd ( nodd <= neven <=nodd+1 )
void DFWT53_1d_spliteo(float *even, int neven, float *odd, int nodd, float *f){
  int i, j ;
  float olo, ohi, elo, ehi ;
  int n = neven + nodd;

  elo = f[0] ;
  olo = f[1] - 0.5f * (f[2] + elo) ;         // first predicted odd term
  for(i = 0, j=0 ; i < n - 2 ; i += 2, j++){           
    ehi = f[i+2] ;                           // upper even term
    ohi = f[i+1]  - 0.5f * (elo  + ehi ) ;   // predict odd term
    odd[j]  = ohi ;                          // store odd term
    even[j] = elo + .25f * (olo + ohi) ;     // update even term and store it
    elo = ehi ;                              
    olo = ohi ;                              
  }                                          
  if(n & 1) {                                // n is odd
    even[neven-1] = f[n-1] + 0.5f * ohi ;          // update and store last even term
  }else{                                     
    ohi = f[n-1] - ehi ;                     // predict last odd term
    odd[nodd]     = ohi ;                          // store last odd term
    even[neven-1] = f[n-2] + .25f * (olo + ohi) ;  // update and store last even term
  }
}

// vectorized version of above
static void DFWT53_1d_spliteo_v(float *even, int neven, float *odd, int nodd, float *f){
  int i, j ;
  float olo, ohi, elo, ehi ;
  int n = neven + nodd;

  for(i=0, j=0 ; i<nodd ; i++, j+=2) { even[i] = f[j] ; odd[i]  = f[j+1] ;}
  if(n & 1) {               // n is odd, last value is even term
    even[neven-1] = f[n-1];

    for(i=0 ; i<nodd ; i++) { odd[i] = odd[i] - 0.5f * (even[i] + even[i+1]) ; }

    even[0] = even[0] + 0.5f * odd[0] ;                    // first even term is special case
    for(i=1 ; i<neven-1 ; i++) { even[i] = even[i] + 0.25f * (odd[i-1] + odd[i]) ; }
    even[neven-1] = even[neven-1] + 0.5f * odd[nodd -1] ;  // last even term is special case
  }else{                   // n is even, last value is odd term
    for(i=0 ; i<nodd-1 ; i++) { odd[i] = odd[i] - 0.5f * (even[i] + even[i+1]) ; }
    odd[nodd -1] = odd[nodd -1] - even[neven-1] ;           // last odd term is special case

    even[0] = even[0] + 0.5f * odd[0] ;                    // first even term is special case
    for(i=1 ; i<neven ; i++) { even[i] = even[i] + 0.25f * (odd[i-1] + odd[i]) ; }
  }
}

// discrete inverse linear lifted wavelet transform (scalar, possibly in place no copy form)
void DIWT53_1d(float *f1, float *f, int n){
  int i ;
  float olo, ohi, elo, ehi ;

  olo = f1[1] ;
  elo = f1[0] - 0.5f * olo ;
  f[0] = elo ;
  for(i = 2 ; i < n - 2 ; i += 2){
    ohi = f1[i+1] ;
    ehi = f1[i] - 0.25f * (olo + ohi) ;
    f[i-1] = olo  + 0.5f * (elo + ehi) ;
    f[i] = ehi ;
    olo = ohi ;
    elo = ehi ;
  }
  if(n & 1) {                                 // n is odd
    ehi    = f1[n-1] - 0.5f * f1[n-2] ;       // unupdate last even
    f[n-1] = ehi                      ;       // store last even
    f[n-2] = f1[n-2] + 0.5f * (elo + ehi) ;   // unpredict and store last odd
  }else{
    ohi = f1[n-1] ;
    ehi = f1[n-2] - 0.25f * (olo + ohi) ;     // unupdate last even
    f[n-2] = ehi ;                            // store last even
    f[n-3] = f1[n-3] + 0.5f * (elo + ehi ) ;  // unpredict and store next to last odd
    f[n-1] = f1[n-1] + ehi ;                  // unpredict and store last odd
  }
}

// discrete inverse linear lifted wavelet transform (scalar form)
// with split even/odd components
// number of points n is neven + nodd ( nodd <= neven <=nodd+1 )
void DIWT53_1d_spliteo(float *even, int neven, float *odd, int nodd, float *f){
  int i, j ;
  float olo, ohi, elo, ehi ;
  int n = neven + nodd;

  olo = odd[0] ;
  elo = even[0] - 0.5f * olo ;
  f[0] = elo ;
  for(i = 2, j = 1 ; i < n - 2 ; i += 2, j++){
    ohi = odd[j] ;
    ehi = even[j] - 0.25f * (olo + ohi) ;
    f[i-1] = olo  + 0.5f * (elo + ehi) ;
    f[i] = ehi ;
    olo = ohi ;
    elo = ehi ;
  }
  if(n & 1) {                                      // n is odd
    ehi    = even[neven-1] - 0.5f * odd[nodd-1] ;  // unupdate last even
    f[n-1] = ehi                      ;            // store last even
    f[n-2] = odd[nodd-1] + 0.5f * (elo + ehi) ;    // unpredict and store last odd
  }else{
    ohi = odd[nodd-1] ;
    ehi = even[neven-1] - 0.25f * (olo + ohi) ;    // unupdate last even
    f[n-2] = ehi ;                                 // store last even
    f[n-3] = odd[nodd-2] + 0.5f * (elo + ehi ) ;   // unpredict and store next to last odd
    f[n-1] = odd[nodd-1] + ehi ;                   // unpredict and store last odd
  }
}

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>

#define NP 12
#define NPP 13

int main(){
  float f[NP], f1[NP], even[NP], odd[NP] ;
  float g[NPP], g1[NPP] ;
  int i ;
  float delta = 1.0f ;
  float start;

  start = 1.0f; f[0] = start ; f[1] = start + delta ;
  for(i=2 ; i<NP ; i+=2){ delta = delta + delta ; f[i] = f[i-1] + delta ; f[i+1] = f[i] + delta ; }

  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n");
  DFWT53_1d(f, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n");
  DIWT53_1d(f, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n\n");
  
  delta = 1.0f ; start = 1.0; f[0] = start ; f[1] = start + delta ;
  for(i=2 ; i<NP ; i+=2){ delta = delta + delta ; f[i] = f[i-1] + delta ; f[i+1] = f[i] + delta ; }

  DFWT53_1d(f1, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f1[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NP ; i++){ f[i] = -1.0f ; }
  DIWT53_1d(f1, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n\n");

  DFWT53_1d_spliteo_v(f1, (NP+1)/2, f1+(NP+1)/2, NP/2, f) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f1[i]) ; } fprintf(stderr,"\n");
  DIWT53_1d_spliteo(f1, (NP+1)/2, f1+(NP+1)/2, NP/2, f) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n\n");

  DFWT53_1d(f1, f, NP) ;
  for(i=3 ; i<NP-2 ; i+=2) f1[i] = 0.0 ;   // zero out odd terms
  DIWT53_1d(f1, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f1[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n\n=======================\n\n");

  start = 1.0f; g[0] = start ; delta = 1.0f; g[1] = start + delta ;
//   for(i=0 ; i<NPP ; i++){ g[i] = i + 1 ; }
  for(i=2 ; i<NP ; i+=2){ delta = delta + delta ; g[i] = g[i-1] + delta ; g[i+1] = g[i] + delta ; }
  g[NPP-1] = g[NPP-2] + delta;

  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n");
  DFWT53_1d(g, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n");
  DIWT53_1d(g, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n\n");

  start = 1.0f; g[0] = start ; delta = 1.0f; g[1] = start + delta ;
  for(i=2 ; i<NP ; i+=2){ delta = delta + delta ; g[i] = g[i-1] + delta ; g[i+1] = g[i] + delta ; }
  g[NPP-1] = g[NPP-2] + delta;

  DFWT53_1d(g1, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g1[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NPP ; i++){ g[i] = -1.0f ; }
  DIWT53_1d(g1, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n\n");

  DFWT53_1d_spliteo_v(g1, (NPP+1)/2, g1+(NPP+1)/2, NPP/2, g) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g1[i]) ; } fprintf(stderr,"\n");
  DIWT53_1d_spliteo(g1, (NPP+1)/2, g1+(NPP+1)/2, NPP/2, g) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n\n");

  DFWT53_1d(g1, g, NPP) ;
  for(i=3 ; i<NPP-2 ; i+=2) g1[i] = 0.0 ;   // zero out odd terms
  DIWT53_1d(g1, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g1[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n");

  return 0;
}
#endif
