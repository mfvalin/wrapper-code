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

// full transform, in place
extern inline void FDWT53i_1D_inplace_full(int32_t *x, uint32_t n){
  int32_t i, em, ep, om, op ;

  if(n<3) return ;           // no transform for 2 points

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
extern inline void FDWT53i_1D_split_full(int32_t *x, int32_t *e, int32_t *o, uint32_t n){
  int32_t i, em, ep, om, op ;

  if(n<3) {
    e[0] = x[0] ;
    if(n == 2) o[0] = x[1] ;
    return ;           // no transform for 2 points
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
extern inline void FDWT53i_1D_inplace_split_full(int32_t *x, uint32_t n){
  int i ;
  int nodd  = n/2 ;
  int neven = (n+1)/2 ;
  int32_t t[nodd] ;

  if(n<3) return ;           // no transform for 2 points

  FDWT53i_1D_split_full(x, x, t, n) ;
  for(i=0 ; i<nodd ; i++) x[neven+i] = t[i] ;
}

// inverse transform for integer data
// step 1 : unupdate even points [-= (2 + x[i-1] + x[i+1])/4]  [-= (2 + odd[i-1] + odd[i])/4]
// step 2 : unpredict odd points [+= (x[i-1] + x[i+1])/2]      [+= (even[i] + even[i+1])/2

// full inverse transform, in place
extern inline void IDWT53i_1D_inplace_full(int32_t *x, uint32_t n){
  int32_t i, em, ep, om, op ;

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
extern inline void IDWT53i_1D_split_full(int32_t *x, int32_t *e, int32_t *o, uint32_t n){
  int32_t i, em, ep, om, op ;

  if(n<3) {
    x[0] = e[0] ;
    if(n == 2) x[1] = o[1] ;
    return ;           // no transform for 2 points
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
extern inline void IDWT53i_1D_inplace_split_full(int32_t *x, uint32_t n){
  int i ;
  int nodd  = n/2 ;
  int neven = (n+1)/2 ;
  int32_t te[neven] ;

  if(n<3) return ;           // no transform for 2 points

  for(i=0 ; i<neven ; i++) te[i] = x[i] ;
  IDWT53i_1D_split_full(x, te, x+neven, n) ;
}

// analyze the four quadrants for min and max values
static void quadrants(int32_t *z, int ni, int lni, int nj, int32_t *min, int32_t *max){
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
  quadrants(x, ni, lni, nj, min, max) ;
if(levels == 1)
  printf("(%5d %5d) (%5d %5d)\n(%5d %5d) (%5d %5d)",min[2],max[2],min[3],max[3],min[0],max[0],min[1],max[1]) ;
if(levels > 1)
  printf("(%5d %5d) (%5d %5d)\n(           ) (%5d %5d)",min[2],max[2],min[3],max[3],min[1],max[1]) ;
printf("  forward level %d done (%d x %d)\n", levels, ni, nj);
  if(levels > 1){
    FDWT53i_2D_split_inplace_n(x, (ni+1)/2, lni, (nj+1)/2, levels-1) ;
  }
}

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

#if defined(SELF_TEST)


// reference transform
static void F_DWT53i_1D_inplace(int32_t *x, int32_t n){
  int i ;

  if(n<3) return ;           // no transform for 2 points

  // predict odd points [-= (x[i-1] + x[i+1])/2]
  // update even points [+= (2 + x[i-1] + x[i+1])/4]
  if(n & 1){          // n is odd, last point is even
    // predict odd points
    for(i=1 ; i<n ; i+=2) x[i] = x[i] - (x[i-1] + x[i+1])/2 ;
    // update even points
    x[0] = x[0] + (2 + x[1] + x[1])/4 ;               // first even point
    for(i=2 ; i<n-2 ; i+=2) x[i] = x[i] + (2 + x[i-1] + x[i+1])/4 ;
    x[n-1] = x[n-1] + (2 + x[n-2] + x[n-2])/4 ;       // last even point
  }else{              // n is even, last point is odd
    // predict odd points (but last)
    for(i=1 ; i<n-1 ; i+=2) x[i] = x[i] - (x[i-1] + x[i+1])/2 ;
    x[n-1] = x[n-1] - (x[n-2] + x[n-2])/2 ;   // last odd point
    // update even points
    x[0] = x[0] + (2 + x[1] + x[1])/4 ;               // first even point
    for(i=2 ; i<n-1 ; i+=2) x[i] = x[i] + (2 + x[i-1] + x[i+1])/4 ;
  }
}

#include <math.h>

#if ! defined(NPTS)
#define NPTS 4
#endif

#if ! defined(LEVELS)
#define LEVELS 1
#endif

#if defined(TEST2D)
#define NI 2
#else
#define NI NPTS
#endif

#define LNI (NPTS+1)
#define NJ NPTS

int main(int argc, char **argv){
  int32_t x[NPTS+1], e[NPTS+1], o[NPTS+1], r[NPTS+1] ;
  int32_t x2d[NJ+1][LNI] ;
  int32_t r2d[NJ+1][LNI] ;
  int i, j, errors ;
  int neven = (NPTS+1)/2 ;
  int nodd = NPTS/2 ;

  for(i=0 ; i<NPTS ; i++) {
//     x[i] = i ;
    x[i] =  2 * sqrtf( (i-NPTS/2.0f)*(i-NPTS/2.0f) ) ;
    r[i] = x[i] ;
  }
  x[NPTS] = -1 ;
  r[NPTS] = -1 ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]) ;
  printf(" 1D original\n\n") ;

//   if(NI < 3) goto test2d ;
  if(argc > 1) goto test2d ;

  for(i=0 ; i<NPTS ; i++) x[i] = r[i] ; x[NPTS] = -1 ;
  FDWT53i_1D_split_full(x, e, o, NPTS) ;
  for(i=0 ; i<neven ; i++) printf("%5d", e[i]) ;
  for(i=0 ; i<nodd  ; i++) printf("%5d", o[i]) ;
  printf("%5d FDWT53i_1D_split_full (%d even, %d odd)\n",x[NPTS], neven, nodd) ;

  for(i=0 ; i<=NPTS ; i++) x[i] = -1 ;
  IDWT53i_1D_split_full(x, e, o, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", r[i]-x[i]) ;
  printf(" IDWT53i_1D_inplace_full (errors)\n") ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]) ;
  printf(" IDWT53i_1D_split_full (restored)\n\n") ;

  for(i=0 ; i<NPTS ; i++) x[i] = r[i] ; x[NPTS] = -1 ;
  FDWT53i_1D_inplace_split_full(x, NPTS) ;
  for(i=0 ; i<NPTS ; i++) printf("%5d", x[i]) ;
  printf("%5d FDWT53i_1D_inplace_split_full (%d even, %d odd)\n",x[NPTS], neven, nodd) ;

  IDWT53i_1D_inplace_split_full(x, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", r[i]-x[i]) ;
  printf(" IDWT53i_1D_inplace_split_full (errors)\n") ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]) ;
  printf(" IDWT53i_1D_inplace_split_full (restored)\n\n") ;

  for(i=0 ; i<NPTS ; i++) x[i] = r[i] ; x[NPTS] = -1 ;
  F_DWT53i_1D_inplace(x, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]) ;
  printf(" F_DWT53i_1D_inplace\n") ;
  for(i=0 ; i<neven ; i++) printf("%5d", x[i+i]) ;
  for(i=0 ; i<nodd  ; i++) printf("%5d", x[i+i+1]) ;
  printf("%5d F_DWT53i_1D_inplace (%d even, %d odd)\n\n",x[NPTS], neven, nodd) ;

  for(i=0 ; i<=NPTS ; i++) x[i] = r[i] ;
  FDWT53i_1D_inplace_full(x, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]) ;
  printf(" FDWT53i_1D_inplace_full\n") ;
  for(i=0 ; i<neven ; i++) printf("%5d", x[i+i]) ;
  for(i=0 ; i<nodd  ; i++) printf("%5d", x[i+i+1]) ;
  printf("%5d FDWT53i_1D_inplace_full (%d even, %d odd)\n\n",x[NPTS], neven, nodd) ;

  IDWT53i_1D_inplace_full(x, NPTS) ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", r[i]-x[i]) ;
  printf(" IDWT53i_1D_inplace_full (errors)\n") ;
  for(i=0 ; i<=NPTS ; i++) printf("%5d", x[i]) ;
  printf(" IDWT53i_1D_inplace_full (restored)\n\n") ;

test2d:
  for(j = NJ-1 ; j >= 0 ; j--){
    for(i = 0 ; i < NI ; i++) {
//       x2d[j][i] = j + i ;
      x2d[j][i] =  2 * sqrtf( (j-NJ/2.0f)*(j-NJ/2.0f) + (i-NI/2.0f)*(i-NI/2.0f)) ;
      if(NI <= 2) x2d[j][i] =  2 * sqrtf( (j-NJ/2.0f)*(j-NJ/2.0f) ) ;
      r2d[j][i] = x2d[j][i] ;
      printf("%5d", x2d[j][i]) ;
    }
    printf("\n");
  }
  printf(" 2D original data\n\n");
  for(i = 0 ; i < NI ; i++) {
    x2d[NJ][i] = -1 ;
  }

//   F_DWT53i_2D_split_inplace((int32_t *) x2d, NI, LNI, NJ) ;
  FDWT53i_2D_split_inplace_n((int32_t *) x2d, NI, LNI, NJ, LEVELS) ;
  for(j = NJ-1 ; j >= 0 ; j--){
    for(i=0 ; i<NI ; i++) printf("%5d", x2d[j][i]) ;
    printf("\n");
  }
  printf(" FDWT53i_2D_split_inplace\n\n");

//   I_DWT53i_2D_split_inplace((int32_t *) x2d, NI, LNI, NJ) ;
  IDWT53i_2D_split_inplace_n((int32_t *) x2d, NI, LNI, NJ, LEVELS) ;
  errors = 0 ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      if(x2d[j][i] != r2d[j][i]) errors++ ;
    }
  }
//   for(j = NJ-1 ; j >= 0 ; j--){
//     for(i=0 ; i<NI ; i++) printf("%5d", x2d[j][i]) ;
//     printf("\n");
//   }
//   printf(" IDWT53i_2D_split_inplace\n");
  printf("restore errors = %d\n", errors);

}
#endif
