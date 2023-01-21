/* Hopefully useful routines for C and FORTRAN
 * Copyright (C) 2021  Recherche en Prevision Numerique
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
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

// forward transform for integer data
// step 1 : predict odd points [-= (x[i-1] + x[i+1])/2]      [-= (even[i] + even[i+1])/2
// step 2 : update even points [+= (2 + x[i-1] + x[i+1])/4]  [+= (2 + odd[i-1] + odd[i])/4]
// mirror symmetry condition : x[-1] = x[1], x[n] = x[n-2]
//                             even[-1] = even[0], even[neven] = even[neven-1]
//                             odd[-1]  = odd[0],  odd[nodd]   = odd[nodd-1]

#define STATIC extern

// full transform, in place
STATIC inline void FDWT53i_1D_inplace_full(int32_t *x, uint32_t n){
  int32_t i, em, ep, om, op ;

  if(n == 2){                // special case for 2 points
    x[1] -= x[0] ;
    x[0] += (1 + x[1]) / 2 ;
  }
  if(n<3) return ;           // no transform for 1 point

  ep = x[2] ;
  em = x[0] ;
  x[1] = om = x[1] - (em + ep) / 2 ;             // predict first odd term
  x[0] = x[0] + (1 + om) / 2 ;                   // update first even term
  for(i=3 ; i<n-1 ; i+=2){
    em = ep ;
    ep = x[i+1] ;
    x[i] = op = x[i] - (em + ep)/2 ;             // predict next odd term
    x[i-1] = em + (2 + om + op) / 4 ;            // update next even term
    om = op ;
  }
  if(n&1){   // n is odd, last term is even
    x[n-1] = ep + (1 + op) / 2 ;                 // update last even term
  }else{     // n is even, last term is odd
    x[n-1] = op = x[n-1] - (ep + ep) / 2 ;       // predict last odd term
    x[n-2] = ep + (2 + om + op) / 4 ;            // update last even term
  }
}

// full transform, split output
STATIC inline void FDWT53i_1D_split_full(int32_t *x, int32_t *e, int32_t *o, uint32_t n){
  int32_t i, em, ep, om, op ;

  if(n == 1) {
    e[0] = x[0] ;
    return ;           // no transform for 1 point
  }
  if(n == 2){                // special case for 2 points
    *o = x[1] - x[0] ;
    *e = x[0] + (1 + *o) / 2 ;
    return ;
  }

  ep = x[2] ;
  em = x[0] ;
  *o = om = x[1] - (em + ep) / 2 ; o++ ;       // predict first odd term
  *e = x[0] + (1 + om) / 2 ;       e++ ;       // update first even term
  for(i=3 ; i<n-1 ; i+=2){
    em = ep ;
    ep = x[i+1] ;
    *o = op = x[i] - (em + ep)/2 ;  o++ ;      // predict next odd term
    *e = em + (2 + om + op) / 4 ;   e++ ;      // update next even term
    om = op ;
  }
  if(n&1){   // n is odd, last term is even
    *e = ep + (1 + op) / 2 ;                   // update last even term
  }else{     // n is even, last term is odd
    *o = op = x[n-1] - (ep + ep) / 2 ;         // predict last odd term
    *e = ep + (2 + om + op) / 4 ;              // update last even term
  }
}

// full transform, in place, split
STATIC inline void FDWT53i_1D_inplace_split_full(int32_t *x, uint32_t n){
  int i ;
  int nodd  = n/2 ;
  int neven = (n+1)/2 ;
  int32_t t[nodd] ;

  if(n == 2){                // special case for 2 points
    x[1] -= x[0] ;
    x[0] += (1 + x[1]) / 2 ;
  }
  if(n<3) return ;           // no transform for 1 point

  FDWT53i_1D_split_full(x, x, t, n) ;
  for(i=0 ; i<nodd ; i++) x[neven+i] = t[i] ;
}

// inverse transform for integer data
// step 1 : unupdate even points [-= (2 + x[i-1] + x[i+1])/4]  [-= (2 + odd[i-1] + odd[i])/4]
// step 2 : unpredict odd points [+= (x[i-1] + x[i+1])/2]      [+= (even[i] + even[i+1])/2

// full inverse transform, in place
STATIC inline void IDWT53i_1D_inplace_full(int32_t *x, uint32_t n){
  int32_t i, em, ep, om, op ;

  if(n == 2){
    x[0] -= (1 + x[1]) / 2 ;
    x[1] += x[0] ;
  }
  if(n<3) return ;           // no transform for 2 points

  om = x[1] ;
  x[0] = em = x[0] - (1 + om) / 2 ;              // un update first even term
  for(i=2 ; i<n-1 ; i+=2){
    op = x[i+1] ;
    x[i] = ep = x[i] - (2 + om + op) / 4 ;       // un update next even term
    x[i-1] = om + (em + ep) / 2 ;                // un predict next odd term
    om = op ;
    em = ep ;
  }
  if(n&1){   // n is odd, last term is even
    x[n-1] = ep = x[n-1] - (1 + op) / 2 ;        // un update last even term
    x[n-2] = x[n-2] + (em + ep) / 2 ;            // un predict last odd term
  }else{     // n is even, last term is odd
    x[n-1] = op + ep ;                           // un predict last odd term
  }
}

// full inverse transform, split input
STATIC inline void IDWT53i_1D_split_full(int32_t *x, int32_t *e, int32_t *o, uint32_t n){
  int32_t i, em, ep, om, op ;

  if(n == 1) {
    x[0] = e[0] ;
    return ;           // no transform for 1 point
  }
  if(n == 2){                // special case for 2 points
    x[0] = *e - (1 + *o) / 2 ;
    x[1] = *o + x[0] ;
    return ;
  }

  om = *o ;                                o++ ;
  x[0] = em = *e - (1 + om) / 2 ;          e++ ; // un update first even term
  for(i=2 ; i<n-1 ; i+=2){
    op = *o ;                              o++ ;
    x[i] = ep = *e - (2 + om + op) / 4 ;   e++ ; // un update next even term
    x[i-1] = om + (em + ep) / 2 ;                // un predict next odd term
    om = op ;
    em = ep ;
  }
  if(n&1){   // n is odd, last term is even
    ep = *e - (1 + op) / 2 ;                     // un update last even term
    o-- ;
    x[n-2] = *o + (em + ep) / 2 ;                // un predict last odd term
    x[n-1] = ep ;                                // delayed store last even term
  }else{     // n is even, last term is odd
    x[n-1] = op + ep ;                           // un predict last odd term
  }
}

// full inverse transform, in place, split
STATIC inline void IDWT53i_1D_inplace_split_full(int32_t *x, uint32_t n){
  int i ;
  int nodd  = n/2 ;
  int neven = (n+1)/2 ;
  int32_t te[neven] ;

  if(n == 2){
    x[0] -= (1 + x[1]) / 2 ;
    x[1] += x[0] ;
  }
  if(n<3) return ;           // no transform for 2 points

  for(i=0 ; i<neven ; i++) te[i] = x[i] ;
  IDWT53i_1D_split_full(x, te, x+neven, n) ;
}

// analyze the four quadrants for min and max values
void IDWT53i_quadrants(int32_t *z, int ni, int lni, int nj, int32_t *min, int32_t *max){
  int i, j, pts[4] ; ;
for (i=0;i<4;i++) pts[i] = 0 ;
  min[0] = min[1] = min[2] = min[3] =  9999999 ;
  max[0] = max[1] = max[2] = max[3] = -9999999 ;
  for (j=0;j<(nj+1)/2;j++) {
    for (i=0;i<(ni+1)/2;i++) {          // LL quadrant
      min[0] = (z[i+lni*j] < min[0]) ? z[i+lni*j] : min[0] ;
      max[0] = (z[i+lni*j] > max[0]) ? z[i+lni*j] : max[0] ;
      pts[0]++ ;
    }
  }
  for (j=0;j<(nj+1)/2;j++) {
    for (i=(ni+1)/2;i<ni;i++) {       // HL quadrant
      min[1] = (z[i+lni*j] < min[1]) ? z[i+lni*j] : min[1] ;
      max[1] = (z[i+lni*j] > max[1]) ? z[i+lni*j] : max[1] ;
      pts[1]++ ;
    }
  }
  for (j=(nj+1)/2;j<nj;j++) {
    for (i=0;i<(ni+1)/2;i++) {          // LH quadrant
      min[2] = (z[i+lni*j] < min[2]) ? z[i+lni*j] : min[2] ;
      max[2] = (z[i+lni*j] > max[2]) ? z[i+lni*j] : max[2] ;
      pts[2]++ ;
    }
  }
  for (j=(nj+1)/2;j<nj;j++) {
    for (i=(ni+1)/2;i<ni;i++) {       // HH quadrant
      min[3] = (z[i+lni*j] < min[3]) ? z[i+lni*j] : min[3] ;
      max[3] = (z[i+lni*j] > max[3]) ? z[i+lni*j] : max[3] ;
      pts[3]++ ;
    }
  }
}

void FDWT53i_2D_split_inplace(int32_t *x, int ni, int lni, int nj){
  int i, j ;
  int32_t *s = x ;
  int nieven = (ni+1) >> 1 ;  // number of even terms in row
  int niodd  = ni >> 1 ;      // number of odd terms in row
  int njeven = (nj+1) >> 1 ;  // number of even rows
  int njodd  = nj >> 1 ;      // number of odd rows
  int32_t o[ni*njodd] ;       // local temporary for odd rows
  int32_t *erow, *orow ;      // pointers to even row and odd row
  int32_t *orow2 ;            // where odd rows will end up after transform

  if(nj == 2){                // special case for 2 rows
    FDWT53i_1D_inplace_split_full(s    , ni) ;   // first even row
    FDWT53i_1D_inplace_split_full(s+lni, ni) ;   // first odd row
    for(i=0 ; i<ni ; i++) {
      s[lni+i] -= s[i] ;                         // predict odd row
      s[i]     += (1 + s[lni+i]) / 2 ;           // update even row
    }
    return ;
  }
  // 1D transform on rows along i
  s    = x ;                           // pointer to first source row
  erow = x ;                           // pointer to first even row
  orow = o ;                           // pointer to first odd row
  FDWT53i_1D_inplace_split_full(s, ni) ;   // first even row done in place
  for(j = 1 ; j < nj ; j++){
    s+= lni ;                          // next source row
    if(j & 1) {                        // odd row, move to temporary storage
      FDWT53i_1D_split_full(s, orow, orow + nieven, ni) ;
      orow += ni ;                     // next odd row storage
    }else{                             // even row
      erow += lni ;                    // next even row storage
      FDWT53i_1D_split_full(s, erow, erow + nieven, ni) ;
    }
  }
  erow = x ;                     // even rows at the bottom of x
  orow = o ;                     // odd rows in temporary storage
  orow2 = x + lni * njeven ;     // odd rows above even rows
  // transforms along i done, odd rows are in o[njodd][ni], even rows are in x[njeven][lni]
  // predict odd rows
  for(j = 0 ; j < njodd-1 ; j++){   // all odd rows but last
    for(i = 0 ; i < ni ; i++) orow2[i] = orow[i] - (erow[i] + erow[i+lni]) / 2 ;
    erow  += lni ;    // next even row
    orow  += ni ;     // next odd row (temporary)
    orow2 += lni ;    // next odd row
  }
  if(njeven == njodd) {         // even mumber of rows, last row is odd
    for(i = 0 ; i < ni ; i++) orow2[i] = orow[i] - (erow[i] + erow[i]) / 2 ;
  }else{                        // odd mumber of rows, last row is even
    for(i = 0 ; i < ni ; i++) orow2[i] = orow[i] - (erow[i] + erow[i+lni]) / 2 ;
  }
  // update even rows
  erow = x ;                     // even rows at the bottom of x
  orow = x + lni * njeven ;      // odd rows above even rows
  for(i = 0 ; i < ni ; i++) erow[i] = erow[i] + (2 + orow[i] + orow[i])/4 ;
  for(j = 1 ; j < njeven-1 ; j++){  // all but first and last row
    erow += lni ;
    orow += lni ;
    for(i = 0 ; i < ni ; i++) erow[i] = erow[i] + (2 + orow[i-lni] + orow[i])/4 ;
  }
  if(njeven == njodd) {         // even mumber of rows, last row is odd
    erow += lni ;
    orow += lni ;
    for(i = 0 ; i < ni ; i++) erow[i] = erow[i] + (2 + orow[i-lni] + orow[i])/4 ;
  }else{                        // odd mumber of rows, last row is even
    erow += lni ;
    for(i = 0 ; i < ni ; i++) erow[i] = erow[i] + (2 + orow[i] + orow[i])/4 ;
  }
}

// multi level forward transform
// next level applied to LL (even-even) quadrant
void FDWT53i_2D_split_inplace_n(int32_t *x, int ni, int lni, int nj, int levels){
  int min[4], max[4] ;
  FDWT53i_2D_split_inplace(x, ni, lni, nj) ;
//   IDWT53i_quadrants(x, ni, lni, nj, min, max) ;
// if(levels == 1)
//   printf("(%5d %5d) (%5d %5d)\n(%5d %5d) (%5d %5d)",min[2],max[2],min[3],max[3],min[0],max[0],min[1],max[1]) ;
// if(levels > 1)
//   printf("(%5d %5d) (%5d %5d)\n(           ) (%5d %5d)",min[2],max[2],min[3],max[3],min[1],max[1]) ;
// printf("  forward level %d done (%d x %d)\n", levels, ni, nj);
  if(levels > 1){
    FDWT53i_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels-1) ;
  }
}

void FDWT53i_8x8_3level(int32_t *x, int lni){
  FDWT53i_2D_split_inplace_n(x, 8, lni, 8, 2) ;
}

// integer inverse wavelet transform
// x    : integer data to be restored from wavelet transform
// ni   : useful row length
// lni  : storage length of rows
// nj   : number of rows
void IDWT53i_2D_split_inplace(int32_t *x, int ni, int lni, int nj){
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

  if(nj == 2){        // special case for 2 rows
    for(i=0 ; i<ni ; i++){
      x[i]     -= (x[lni+i] + 1) / 2 ;         // unupdate even row
      x[lni+i] += x[i] ;                       // unpredict odd row
    }
    IDWT53i_1D_inplace_split_full(x    , ni) ; // inverse 1D transform on row
    IDWT53i_1D_inplace_split_full(x+lni, ni) ; // inverse 1D transform on row
    return ;
  }
  // "save" even rows
  erow = x ;
  erow2 = e ;
  for(j=0 ; j<njeven ; j++){
    for(i=0 ; i<ni ; i++) erow2[i] = erow[i] ;
    erow2 += ni ;   // next row in e
    erow  += lni ;  // next row in x
  }
  // unupdate even rows
  erow2 = e ;                    // saved even rows
  orow  = x + njeven*lni ;       // odd rows in x (top part)
  for(i=0 ; i<ni ; i++) erow2[i] = erow2[i] - (2 + orow[i] + orow[i]) / 4 ;
  for(j=1 ; j<njeven-1 ; j++){    // all but last even row
    erow2 += ni ;   // next row in e
    orow  += lni ;  // next row in x
    for(i=0 ; i<ni ; i++) erow2[i] = erow2[i] - (2 + orow[i] + orow[i-lni]) / 4 ;
  }
  if(njeven > njodd){          // n is odd, last row is even
    erow2 += ni ;   // next row in e
    for(i=0 ; i<ni ; i++) erow2[i] = erow2[i] - (2 + orow[i] + orow[i]) / 4 ;
  }else{                       // n is even, last row is odd
    erow2 += ni ;   // next row in e
    orow  += lni ;  // next row in x
    for(i=0 ; i<ni ; i++) erow2[i] = erow2[i] - (2 + orow[i] + orow[i-lni]) / 4 ;
  }
  // unpredict odd rows
  orow2 = x + njeven*lni ;       // odd rows in x (top part)
  orow = x + lni ;               // odd rows in x
  erow2 = e ;                    // saved even rows
  erow = x ;                     // even rows in x (bottom part)
  for(j=0 ; j<njodd-1 ; j++){    // all but last odd row
    for(i=0 ; i<ni ; i++) orow[i] = orow2[i] + (erow2[i] + erow2[i+ni])/2 ;
    for(i=0 ; i<ni ; i++) erow[i] = erow2[i] ;
    orow  += lni*2 ;
    orow2 += lni ;
    erow  += lni*2 ;
    erow2 += ni ;
  }
  if(njeven > njodd){          // n is odd, last row is even
    for(i=0 ; i<ni ; i++) orow[i] = orow2[i] + (erow2[i] + erow2[i+ni])/2 ;
    for(i=0 ; i<ni ; i++) erow[i] = erow2[i] ;
    for(i=0 ; i<ni ; i++) erow[i+lni*2] = erow2[i+ni] ;
  }else{                       // n is even, last row is odd
    for(i=0 ; i<ni ; i++) orow[i] = orow2[i] + (erow2[i] + erow2[i])/2 ;
    for(i=0 ; i<ni ; i++) erow[i] = erow2[i] ;
  }

  // inverse 1D transform on rows along i
  row = x ;
  for(j=0 ; j<nj ; j++){
    IDWT53i_1D_inplace_split_full(row, ni) ;
    row += lni ;
  }
}

// multi level inverse transform
// next level applied to LL (even-even) quadrant
void IDWT53i_2D_split_inplace_n(int32_t *x, int ni, int lni, int nj, int levels){
  if(levels > 1){
    IDWT53i_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels-1) ;
  }
  IDWT53i_2D_split_inplace(x, ni, lni, nj) ;
printf("inverse level %d done (%d x %d)\n", levels, ni, nj);
}

void IDWT53i_8x8_3level(int32_t *x, int lni){
  IDWT53i_2D_split_inplace_n(x, 8, lni, 8, 2) ;
}
