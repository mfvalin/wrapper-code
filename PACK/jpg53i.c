/* 
 * Copyright (C) 2020  Recherche en Prevision Numerique
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
/*
   integer forward and inverse wavelet transforms (as in jpeg 2000)

   this code is using a lifting implementation
   https://en.wikipedia.org/wiki/Lifting_scheme

   1 dimensional transform, "in place", "even/odd split", or "in place with even/odd split"
     original data
     +--------------------------------------------------------+
     |                  N data                                |
     +--------------------------------------------------------+

     transformed data (in place, no split, even number of data)
     +--------------------------------------------------------+
     |   N data, even/odd, even/odd, ..... , even/odd         +
     +--------------------------------------------------------+

     transformed data (in place, no split, odd number of data)
     +--------------------------------------------------------+
     |   N data, even/odd, even/odd, ..... , even/odd, even   +
     +--------------------------------------------------------+

     transformed data, in place with even/odd split
     +--------------------------------------------------------+
     | (N+1)/2 even data            |    (N/2) odd data       |
     +--------------------------------------------------------+

     original data                     transformed data (2 output arrays)
     +------------------------------+  +-------------------+  +------------------+
     |             N data           |  | (N+1)/2 even data |  |   N/2 odd data   |
     +------------------------------+  +-------------------+  +------------------+

     even data are the "approximation" terms ("low frequency" terms)
     odd data are the "detail" terms         ("high frequency" terms)

   2 dimensional in place with 2 D split
     original data                               transformed data (in same array)
     +------------------------------------+      +-------------------+----------------+
     |                  ^                 |      +                   |                |
     |                  |                 |      +   even i/odd j    |  odd i/odd j   |
     |                  |                 |      +                   |                |
     |                  |                 |      +                   |                |
     |                  |                 |      +-------------------+----------------+
     |               NJ data              |      +                   |                |
     |                  |                 |      +                   |                |
     |                  |                 |      +   even i/even j   |  odd i/even j  |
     |<----- NI data ---|---------------->|      +                   |                |
     |                  v                 |      +                   |                |
     +------------------------------------+      +-------------------+----------------+
     the process can be applied again to the even/even transformed part to achieve a multi level transform
*/
#include <stdio.h>
#include <stdint.h>
#include <math.h>

void F_dwt53i_1D_inplace(int32_t *x, int32_t n){
  int i ;

  if(n<3) return ;           // no transform for 2 points

  x[0] = x[0] - (2 + x[1] + x[1])/4 ;
  if(n & 1){          // n is odd
    for(i=2 ; i<n-1 ; i+=2) x[i] = x[i] - (2 + x[i-1] + x[i+1])/4 ;
    x[n-1] = x[n-1] - (2 + x[n-2] + x[n-2])/4 ;
    for(i=1 ; i<n ; i+=2) x[i] = x[i] + (x[i-1] + x[i+1])/2 ;
  }else{              // n is even
    for(i=2 ; i<n ; i+=2) x[i] = x[i] - (2 + x[i-1] + x[i+1])/4 ;
    for(i=1 ; i<n-1 ; i+=2) x[i] = x[i] + (x[i-1] + x[i+1])/2 ;
    x[n-1] = x[n-1] + (x[n-2] + x[n-2])/2 ;
  }
}

void F_dwt53i_1D_split(int32_t *x, int32_t *e, int32_t *o, int32_t n){
  int i, j ;
  int neven = (n+1) >> 1 ;
  int nodd = n >> 1 ;

  if(n<3) {
    e[0] = x[0] ;
    if(n == 2) o[0] = x[1] ;
    return ;           // no transform for 2 points
  }

  for(i=0 , j=0 ; j<nodd ; j++, i+=2){
    e[j] = x[i] ;
    o[j] = x[i+1] ;
  }
  e[0] = e[0] - (2 + o[0] + o[0])/4 ;
  if(neven > nodd) {         // n is odd
    e[neven-1] = x[n-1] ;
    for(i = 1 ; i < neven-1 ; i++) e[i] = e[i] - (2 + o[i-1] + o[i])/4 ;
    e[neven-1] = e[neven-1] - (2 + o[nodd-1] + o[nodd-1])/4 ;
    for(i = 0 ; i < nodd ; i++) o[i] = o[i] + (e[i] + e[i+1])/2 ;
  }else{                     // n is even
    for(i = 1 ; i < neven ; i++) e[i] = e[i] - (2 + o[i-1] + o[i])/4 ;
    for(i = 0 ; i < nodd-1 ; i++) o[i] = o[i] + (e[i] + e[i+1])/2 ;
    o[nodd-1] = o[nodd-1] + (e[neven-1] + e[neven-1])/2 ;
  }
}

void F_dwt53i_1D_split_inplace(int32_t *x, int32_t n){
  int i, j ;
  int neven = (n+1) >> 1 ;
  int nodd = n >> 1 ;
  int32_t o[n] ;             // local temporary for odd terms
  int32_t *e = x ;

  if(n<3) return ;           // no transform for 2 points

  for(i=0 , j=0 ; j<nodd ; j++, i+=2){
    e[j] = x[i] ;                       // collect even terms
    o[j] = x[i+1] ;                     // collect odd terms
  }
  e[0] = e[0] - (2 + o[0] + o[0])/4 ;                                      // predict first even term
  if(neven > nodd) {         // n is odd
    e[neven-1] = x[n-1] ;               // collect last even term
    for(i = 1 ; i < neven-1 ; i++) e[i] = e[i] - (2 + o[i-1] + o[i])/4 ;   // predict even terms
    e[neven-1] = e[neven-1] - (2 + o[nodd-1] + o[nodd-1])/4 ;
    for(i = 0 ; i < nodd ; i++) o[i] = o[i] + (e[i] + e[i+1])/2 ;          // update odd terms
  }else{                     // n is even
    for(i = 1 ; i < neven ; i++) e[i] = e[i] - (2 + o[i-1] + o[i])/4 ;     // predict even terms
    for(i = 0 ; i < nodd-1 ; i++) o[i] = o[i] + (e[i] + e[i+1])/2 ;        // update odd terms
    o[nodd-1] = o[nodd-1] + (e[neven-1] + e[neven-1])/2 ;                  // last term is odd
  }
  for(j=0 ; j<nodd ; j++) x[neven+j] = o[j] ;
}

void F_dwt53i_2D_split_inplace(int32_t *x, int ni, int lni, int nj){
  int i, j ;
  int32_t *s = x ;
  int nieven = (ni+1) >> 1 ;  // number of even terms in row
  int niodd  = ni >> 1 ;      // number of odd terms in row
  int njeven = (nj+1) >> 1 ;  // number of even rows
  int njodd  = nj >> 1 ;      // number of odd rows
  int32_t o[ni*njodd] ;       // local temporary for odd rows
  int32_t *erow, *orow ;      // pointers to even row and odd row
  int32_t *orow2 ;            // where odd rows will end up after transform

  // 1D transform on rows along i
  s    = x ;                           // pointer to first source row
  erow = x ;                           // pointer to first even row
  orow = o ;                           // pointer to first odd row
  F_dwt53i_1D_split_inplace(s, ni) ;   // first even row
  for(j = 1 ; j < nj ; j++){
    s+= lni ;                          // next source row
    if(j & 1) {                        // odd row, move to temporary storage
      F_dwt53i_1D_split(s, orow, orow + nieven, ni) ;
      orow += ni ;                     // next odd row storage
    }else{                             // even row
      erow += lni ;                    // next even row storage
      F_dwt53i_1D_split(s, erow, erow + nieven, ni) ;
    }
  }
  erow = x ;
  orow = o ;
  // transforms along i done, odd rows are in o[njodd][ni], even rows are in x[njeven][lni]
  // predict even rows
  for(i = 0 ; i < ni ; i++) erow[i] = erow[i] - (2 + orow[i] + orow[i]) / 4 ;  // first even row
  orow += ni ;
  for(j = 1 ; j < njeven-1 ; j++){  // predict even rows (all but first and last)
    erow += lni ;
    for(i = 0 ; i < ni ; i++) erow[i] = erow[i] - (2 + orow[i-ni] + orow[i]) / 4 ;
    orow += ni ;
  }
  erow += lni ;
  if(njeven == njodd) {         // even mumber of rows, last row is odd
    for(i = 0 ; i < ni ; i++) erow[i] = erow[i] - (2 + orow[i-ni] + orow[i]) / 4 ;
  }else{                        // odd number of rows, last row is even
    for(i = 0 ; i < ni ; i++) erow[i] = erow[i] - (2 + orow[i-ni] + orow[i-ni]) / 4 ;
  }
  // update odd rows
  erow  = x ;
  orow  = o ;
  orow2 = x + lni * njeven ;     // store odd rows on top of even rows
  for(j = 0 ; j < njodd-1 ; j++){   // update all odd rows but last
    for(i = 0 ; i < ni ; i++) orow2[i] = orow[i] + ( erow[i] +erow[i+lni]  ) / 2 ;
    erow  += lni ;
    orow  += ni ;
    orow2 += lni ;
  }
  if(njeven == njodd) {         // even mumber of rows, last row is odd
    for(i = 0 ; i < ni ; i++) orow2[i] = orow[i] + ( erow[i] +erow[i]  ) / 2 ;
  }else{                        // odd number of rows, last row is even
    for(i = 0 ; i < ni ; i++) orow2[i] = orow[i] + ( erow[i] +erow[i+lni]  ) / 2 ;
  }
}

void F_dwt53i_2D_split_inplace_n(int32_t *x, int ni, int lni, int nj, int levels){
  F_dwt53i_2D_split_inplace(x, ni, lni, nj) ;
printf("forward level %d done (%d x %d)\n", levels, ni, nj);
  if(levels > 1){
    F_dwt53i_2D_split_inplace_n(x + (ni+1)/2 + lni *((nj+1)/2), (ni)/2, lni, (nj)/2, levels-1) ;
//     F_dwt53i_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels-1) ;
  }
}

void I_dwt53i_1D_inplace(int32_t *x, int32_t n){
  int i ;

  if(n<3) return ;           // no transform for 2 points

  if(n & 1){          // n is odd
    for(i=1 ; i<n ; i+=2) x[i] = x[i] - (x[i-1] + x[i+1])/2 ;
    x[0] = x[0] + (2 + x[1] + x[1])/4 ;
    for(i=2 ; i<n-1 ; i+=2) x[i] = x[i] + (2 + x[i-1] + x[i+1])/4 ;
    x[n-1] = x[n-1] + (2 + x[n-2] + x[n-2])/4 ;
  }else{              // n is even
    for(i=1 ; i<n-1 ; i+=2) x[i] = x[i] - (x[i-1] + x[i+1])/2 ;
    x[n-1] = x[n-1] - (x[n-2] + x[n-2])/2 ;
    x[0] = x[0] + (2 + x[1] + x[1])/4 ;
    for(i=2 ; i<n ; i+=2) x[i] = x[i] + (2 + x[i-1] + x[i+1])/4 ;
  }
}

void I_dwt53i_1D_split(int32_t *x, int32_t *e, int32_t *o, int32_t n){
  int i ;
  int neven = (n+1) >> 1 ;
  int nodd = n >> 1 ;

  if(n<3) {
    x[0] = e[0] ;
    if(n == 2) x[1] = o[1] ;
    return ;           // no transform for 2 points
  }

  if(neven > nodd){          // n is odd
    for(i = 0 ; i < nodd ; i++) x[i+i+1] = o[i] - (e[i] + e[i+1])/2 ;
    x[0] = e[0] + (2 + x[1] + x[1])/4 ;
    for(i = 1 ; i < neven-1 ; i++) x[i+i] = e[i] + (2 + x[i+i-1] + x[i+i])/4 ;
    x[n-1] = e[neven-1] + (2 + x[n-2] + x[n-2])/4 ;
  }else{
    for(i = 0 ; i < nodd-1 ; i++) x[i+i+1] = o[i] - (e[i] + e[i+1])/2 ;
    x[n-1] = o[nodd-1] - (e[neven-1] + e[neven-1])/2 ;
    x[0] = e[0] + (2 + x[1] + x[1])/4 ;
    for(i = 1 ; i < neven ; i++) x[i+i] = e[i] + (2 + x[i+i-1] + x[i+i+1])/4 ;
  }
}

void I_dwt53i_1D_split_inplace(int32_t *x, int32_t n){
  int i ;
  int neven = (n+1) >> 1 ;
  int nodd = n >> 1 ;
  int32_t e[n] ;             // local temporary for even terms
  int32_t *o ;

  if(n<3) return ;           // no transform for 2 points

  o = x + neven ;
  for(i = 0 ; i < neven ; i++) e[i] = x[i] ;   // copy even terms

  if(neven > nodd){          // n is odd
    for(i = 0 ; i < nodd ; i++) x[i+i+1] = o[i] - (e[i] + e[i+1])/2 ;
    x[0] = e[0] + (2 + x[1] + x[1])/4 ;
    for(i = 1 ; i < neven-1 ; i++) x[i+i] = e[i] + (2 + x[i+i-1] + x[i+i+1])/4 ;
    x[n-1] = e[neven-1] + (2 + x[n-2] + x[n-2])/4 ;
  }else{
    for(i = 0 ; i < nodd-1 ; i++) x[i+i+1] = o[i] - (e[i] + e[i+1])/2 ;
    x[n-1] = o[nodd-1] - (e[neven-1] + e[neven-1])/2 ;
    x[0] = e[0] + (2 + x[1] + x[1])/4 ;
    for(i = 1 ; i < neven ; i++) x[i+i] = e[i] + (2 + x[i+i-1] + x[i+i+1])/4 ;
  }
}

void I_dwt53i_2D_split_inplace(int32_t *x, int ni, int lni, int nj){
  int i, j ;
  int nieven = (ni+1) >> 1 ;  // number of even terms in row
  int niodd  = ni >> 1 ;      // number of odd terms in row
  int njeven = (nj+1) >> 1 ;  // number of even rows
  int njodd  = nj >> 1 ;      // number of odd rows
  int32_t e[ni*njeven] ;      // local temporary for even rows
  int32_t *erow, *orow ;      // pointers to even row and odd row
  int32_t *erow2 ;            // temporary pointer to even row storage
  int32_t *orow2 ;            // temporary pointer to odd row storage
  int32_t *row ;              // temporary pointer to row storage

  // "save" even rows
  erow = x ;
  erow2 = e ;
  for(j=0 ; j<njeven ; j++){
    for(i=0 ; i<ni ; i++) erow2[i] = erow[i] ;
    erow2 += ni ;   // next row in e
    erow  += lni ;  // next row in x
  }
  // unupdate odd rows
  orow  = x + njeven*lni ;       // odd rows in x (top part)
  orow2 = x + lni ;              // first odd row in x (from bottom)
  erow  = e ;                    // first even row
  for(j=0 ; j<njodd-1 ; j++){    // unupdate all but last odd row
    for(i=0 ; i<ni ; i++) orow2[i] = orow[i] - (erow[i] + erow[i+ni]) / 2 ;
    orow  += lni ;
    orow2 += (lni*2) ;
    erow  += ni ;
  }
  if(nj & 1){                     // nj is odd, top row is an even row
    for(i=0 ; i<ni ; i++) orow2[i] = orow[i] - (erow[i] + erow[i+ni]) / 2 ;
  }else{                          // nj is even, top row is an odd row
    for(i=0 ; i<ni ; i++) orow2[i] = orow[i] - (erow[i] + erow[i]) / 2 ;
  }
  // unpredict enven rows
  erow2 = x ;                     // first even row in x
  erow  = e ;                     // first even row in e
  orow  = x + lni ;               // first odd row in x
  for(i=0 ; i<ni ; i++) erow2[i] = erow[i] + (2 + orow[i] + orow[i])/4 ;
  for(j=1 ; j<njeven-1 ; j++){    // all but last even row
    erow  += ni ;                 // next even row in e
    erow2 += (lni*2) ;            // next even row in x
    orow  += (lni*2) ;            // next odd row in x
    for(i=0 ; i<ni ; i++) erow2[i] = erow[i] + (2 + orow[i] + orow[i-lni*2])/4 ;
  }
  erow2 += (lni*2) ;              // next even row in x
  erow  += ni ;                   // next even row in e
  if(nj & 1){                     // nj is odd, top row is an even row
    for(i=0 ; i<ni ; i++) erow2[i] = erow[i] + (2 + orow[i] + orow[i])/4 ;
  }else{                          // nj is even, top row is an odd row
    for(i=0 ; i<ni ; i++) erow2[i] = erow[i] + (2 + orow[i] + orow[i+lni*2])/4  ;
  }
  // inverse 1D transform on rows along i
  row = x ;
  for(j=0 ; j<nj ; j++){
    I_dwt53i_1D_split_inplace(row, ni) ;
    row += lni ;
  }
}

void I_dwt53i_2D_split_inplace_n(int32_t *x, int ni, int lni, int nj, int levels){
  if(levels > 1){
    I_dwt53i_2D_split_inplace_n(x + (ni+1)/2 + lni *((nj+1)/2), (ni)/2, lni, (nj)/2, levels-1) ;
//     I_dwt53i_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels-1) ;
  }
  I_dwt53i_2D_split_inplace(x, ni, lni, nj) ;
printf("inverse level %d done (%d x %d)\n", levels, ni, nj);
}


#if defined(SELF_TEST)

#if ! defined(NPTS)
#define NPTS 4
#endif

#if ! defined(LEVELS)
#define LEVELS 1
#endif

#define NI NPTS
#define LNI (NPTS+1)
#define NJ NPTS

int main(int argc, char **argv){
  int32_t x[NPTS+1], e[NPTS+1], o[NPTS+1], r[NPTS+1] ;
  int32_t x2d[NJ+1][LNI] ;
  int32_t r2d[NJ+1][LNI] ;
  int i, j, errors ;
  int neven = (NPTS+1)/2 ;
  int nodd = NPTS/2 ;
#if defined(TEST1D)
  for(i=0 ; i<NPTS ; i++) {
//     x[i] = i ;
    x[i] =  2 * sqrtf( (i-NPTS/2.0f)*(i-NPTS/2.0f) ) ;
    r[i] = x[i] ;
  }
  x[NPTS] = -1 ;
  r[NPTS] = -1 ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]) ;
  printf(" original\n\n") ;
  F_dwt53i_1D_split(x, e, o, NPTS) ;
  for(i=0 ; i<NPTS/2 ; i++) printf("%5d%5d", e[i], o[i]) ;
  if(NPTS&1) printf("%5d", e[NPTS/2]);
  printf("      F_dwt53i_1D_split (even/odd)\n") ;

  I_dwt53i_1D_split(x, e, o, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]-r[i]) ;
  printf(" I_dwt53i_1D_split (restored)\n\n") ;

  for(i=0 ; i<NPTS ; i++) x[i] = r[i] ; x[NPTS] = -1 ;
  F_dwt53i_1D_inplace(x, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]) ;
  printf(" F_dwt53i_1D_inplace\n") ;
  for(i=0 ; i<neven ; i++) printf("%5d", x[i+i]) ;
  for(i=0 ; i<nodd  ; i++) printf("%5d", x[i+i+1]) ;
  printf("%5d F_dwt53i_1D_inplace (%d even, %d odd)\n",x[NPTS], neven, nodd) ;

  I_dwt53i_1D_inplace(x, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", r[i]-x[i]) ;
  printf(" I_dwt53i_1D_inplace (restored)\n\n") ;

  for(i=0 ; i<NPTS ; i++) x[i] = r[i] ; x[NPTS] = -1 ;
  F_dwt53i_1D_split_inplace(x, NPTS) ;
  for(i=0 ; i<nodd ; i++) printf("%5d%5d", x[i], x[neven+i]) ;
  if(NPTS&1) printf("%5d", x[NPTS/2]);
  printf("%5d F_dwt53i_1D_split_inplace (interlaced)\n",x[NPTS]) ;
  for(i=0 ; i<NPTS ; i++) printf("%5d", x[i]);
  printf("%5d F_dwt53i_1D_split_inplace (straight)\n",x[NPTS]) ;

  I_dwt53i_1D_split_inplace(x, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]-r[i]) ;
  printf(" I_dwt53i_1D_split_inplace (restored)\n\n") ;
#endif
  for(j = NJ-1 ; j >= 0 ; j--){
    for(i = 0 ; i < NI ; i++) {
//       x2d[j][i] = j + i ;
      x2d[j][i] =  2 * sqrtf( (j-NJ/2.0f)*(j-NJ/2.0f) + (i-NI/2.0f)*(i-NI/2.0f)) ;
//       x2d[j][i] =  2 * sqrtf( (j-NJ/2.0f)*(j-NJ/2.0f) ) ;
      r2d[j][i] = x2d[j][i] ;
//       printf("%5d", x2d[j][i]) ;
    }
//     printf("\n");
  }
//   printf(" original data\n\n");
  for(i = 0 ; i < NI ; i++) {
    x2d[NJ][i] = -1 ;
  }

  F_dwt53i_2D_split_inplace_n((int32_t *) x2d, NI, LNI, NJ, LEVELS) ;
  for(j = NJ-1 ; j >= 0 ; j--){
    for(i=0 ; i<NI ; i++) printf("%5d", x2d[j][i]) ;
    printf("\n");
  }
  printf(" F_dwt53i_2D_split_inplace\n\n");

  I_dwt53i_2D_split_inplace_n((int32_t *) x2d, NI, LNI, NJ, LEVELS) ;
  errors = 0 ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      if(x2d[j][i] != r2d[j][i]) errors++ ;
    }
  }
  if(errors == 0){
//     for(j = NJ-1 ; j >= 0 ; j--){
//       for(i=0 ; i<NI ; i++) printf("%5d", x2d[j][i]) ;
//       printf("\n");
//     }
  }else{
    for(j = NJ-1 ; j >= 0 ; j--){
      for(i=0 ; i<NI ; i++) printf("%5d", x2d[j][i]-r2d[j][i]) ;
      printf("\n");
    }
  }
  printf(" I_dwt53i_2D_split_inplace\n");
  printf("restore errors = %d\n", errors);
}
#endif
