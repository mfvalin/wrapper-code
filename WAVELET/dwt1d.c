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

// discrete forward linear lifted wavelet transform (scalar no copy form)
void DFWT53_1d(float *f1, float *f, int n){
  int i ;
  float olo, ohi, elo, ehi ;

  elo = f[0] ;
  olo = f[1] - 0.5f * (f[2] + elo) ;       // f1rst predicted odd term
  for(i = 0 ; i < n - 2 ; i += 2){
    ehi = f[i+2] ;                         // upper even term
    ohi = f[i+1]  - 0.5f * (elo  + ehi ) ; // predict odd term
    f1[i+1] = ohi ;                        // store odd term
    f1[i  ] = elo + .25f * (olo + ohi) ;   // update even term and store it
    elo = ehi ;
    olo = ohi ;
  }
  if(n & 1) {                              // n is odd
    f1[n-1] = f[n-1] + 0.5f * f1[n-2] ;    // update and store last even term
  }else{
    f1[n-1] = f[n-1] - f[n-2] ;                      // predict and store last odd term
    f1[n-2] = f[n-2] - .25f * (f1[n-3] + f1[n-1]) ;  // update and store last even term
  }
}

// discrete inverse linear lifted wavelet transform (scalar no copy form)
void DIWT_1d_full(float *f1, float *f, int n){
  int i ;
  float olo, ohi, elo, ehi ;

}

#if defined(SELF_TEST)
#include <stdio.h>
#include <stdlib.h>

#define NP 8
#define NPP 9

int main(){
  float f[NP], f1[NP] ;
  float g[NPP], g1[NPP] ;
  int i ;
  float delta = 1.0f ;

  for(i=0 ; i<NP ; i++){ f[i] = i + 1 ; }
  DFWT53_1d(f1, f, NP) ;
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NP ; i++){ fprintf(stderr,"%8f ",f1[i]) ; } fprintf(stderr,"\n\n");

  for(i=0 ; i<NPP ; i++){ g[i] = i + 1 ; }
  DFWT53_1d(g1, g, NPP) ;
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g[i]) ; } fprintf(stderr,"\n");
  for(i=0 ; i<NPP ; i++){ fprintf(stderr,"%8f ",g1[i]) ; } fprintf(stderr,"\n");

  return 0;
}
#endif
