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

// discrete forward linear lifted wavelet transform (scalar in place no copy form)
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

// discrete inverse linear lifted wavelet transform (scalar in place no copy form)
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

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>

#define NP 12
#define NPP 13

int main(){
  float f[NP], f1[NP] ;
  float g[NPP], g1[NPP] ;
  int i ;
  float delta = 1.0f ;
  float start;

  start = 1.0; f[0] = start ; f[1] = start + delta ;
//   for(i=0 ; i<NP ; i++){ f[i] = i + 1 ; }
  for(i=2 ; i<NP ; i+=2){ delta = delta + delta ; f[i] = f[i-1] + delta ; f[i+1] = f[i] + delta ; }
  DFWT53_1d(f1, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f1[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NP ; i++){ f[i] = -1.0f ; }
  DIWT53_1d(f1, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n\n");
  for(i=3 ; i<NP-2 ; i+=2) f1[i] = 0.0 ;   // zero out odd terms
  DIWT53_1d(f1, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f1[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n\n\n");

  start = 1.0; g[0] = start ; delta = 1.0; g[1] = start + delta ;
//   for(i=0 ; i<NPP ; i++){ g[i] = i + 1 ; }
  for(i=2 ; i<NP ; i+=2){ delta = delta + delta ; g[i] = g[i-1] + delta ; g[i+1] = g[i] + delta ; }
  g[NPP-1] = g[NPP-2] + delta;
  DFWT53_1d(g1, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g1[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NPP ; i++){ g[i] = -1.0f ; }
  DIWT53_1d(g1, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n\n");
  DIWT53_1d(g1, g, NPP) ;
  for(i=3 ; i<NPP-2 ; i+=2) g1[i] = 0.0 ;   // zero out odd terms
  DIWT53_1d(g1, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g1[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n");

  return 0;
}
#endif
