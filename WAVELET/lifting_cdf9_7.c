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
#if defined(ORIGINAL)
#define A    (-1.58615986717275f)
#define B    (-0.05297864003258f)
#define C      0.88293362717904f
#define D      0.44350482244527f
#define S      1.14960430535816f
#define Z      0.86986452237446f
#else
#define A    (-1.586134342f)
#define B    (-0.05298011854f)
#define C      0.8829110762f
#define D      0.4435068522f
#define S      1.149604398f
#define Z      0.869864452f
#endif

// Cohen-Daubechies-Favreau 9/7 wavelets
// lifting implementation
// in place or  even/odd split 

// Forward DWT transform (analysis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[(n+1)/2]  : odd coefficients of the transform (detail)
void F_CDF97_1D_split_N_even(float *x, float *e, float *o, int n){
  int i;
  int neven = (n+1) >> 1;
  int nodd  = neven;

  for(i = 0 ; i < nodd-1 ; i++) o[i] = x[i+i+1] + A * (x[i+i] + x[i+i+2]);
  o[nodd-1] = x[n-1] + 2 * A * x[n-2];  

  e[0 ] = x[0] + 2 * B * o[0];
  for(i = 1; i < neven ; i++) e[i] = x[i+i] + B * (o[i] + o[i-1]);

  for(i = 0 ; i < nodd-1 ; i++) o[i] +=  C * (e[i] + e[i+1]);
  o[nodd-1] +=  2 * C * e[neven-1];

  e[0] = S * (e[0] + 2 * D * o[0]);
  for(i = 1; i < neven ; i++) { e[i] = S * (e[i] +  D * (o[i] + o[i-1])) ; o[i-1] *= (-Z); }
  o[nodd-1] *= (-Z);
}

// Forward DWT transform (analysis) (in place, x is overwritten)
// n           : number of data points (even)
// x[n]        : input data
void F_CDF97_1D_inplace_N_even(float *x, int n){
  int i;
  
  for (i = 1; i < n - 2; i += 2) x[i] += A * (x[i-1] + x[i+1]);  // predict odd terms #1
  x[n-1] += 2 * A * x[n-2];                                      // last term is odd
  
  x[0] += 2 * B * x[1];                                          // update even terms #1
  for (i = 2; i < n; i += 2) x[i] += B * (x[i+1] + x[i-1]);
  
  for (i = 1; i < n - 2; i += 2) x[i] += C * (x[i-1] + x[i+1]);  // predict odd terms #2
  x[n-1] += 2 * C * x[n-2];                                      // last term is odd
  
  x[0] = S * (x[0] + 2 * D * x[1]);                              // update even terms #2 and scale
  for (i = 2; i < n; i += 2) { x[i] = S * (x[i] +  D * (x[i+1] + x[i-1])) ; x[i-1] *= (-Z); }
  x[n-1] *= (-Z);                                                // scale last (odd) term
}

// Forward DWT transform (analysis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
void F_CDF97_1D_split_N_odd(float *x, float *e, float *o, int n){
  int i;
  int neven = (n+1) >> 1;
  int nodd  = n >> 1;

  for(i = 0 ; i < nodd ; i++) o[i] = x[i+i+1] + A * (x[i+i] + x[i+i+2]);

  e[0 ] = x[0] + 2 * B * o[0];
  for(i = 1; i < neven-1 ; i++) e[i] = x[i+i] + B * (o[i] + o[i-1]);
  e[neven-1] = x[n-1] + 2 * B * o[nodd-1];

  for(i = 0 ; i < nodd ; i++) o[i] +=  C * (e[i] + e[i+1]);

  e[0] = S * (e[0] + 2 * D * o[0]);
  for(i = 1; i < neven-1 ; i++) { e[i] = S * (e[i] +  D * (o[i] + o[i-1])) ; o[i-1] *= (-Z); }
  e[neven-1] = S * (e[neven-1] + 2 * D * o[nodd-1]);
  o[nodd-1] *= (-Z);
}

// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (odd)
// x[n]        : input data
void F_CDF97_1D_inplace_N_odd(float *x, int n){
  int i;

  for (i = 1; i < n - 1; i += 2) x[i] += A * (x[i-1] + x[i+1]);  // predict odd terms #1
  
  x[0] += 2 * B * x[1];                                          // update even terms #1
  for (i = 2; i < n - 2; i += 2) x[i] += B * (x[i+1] + x[i-1]);
  x[n - 1] += 2 * B * x[n - 2];                                  // last term is even

  for (i = 1; i < n - 1; i += 2) x[i] += C * (x[i-1] + x[i+1]);  // predict odd terms #2
  
  x[0] = S * (x[0] + 2 * D * x[1]);                              // update even terms #2 and scale
  for (i = 2; i < n - 2; i += 2) { x[i] = S * (x[i] + D * (x[i+1] + x[i-1]));  x[i-1] *= (-Z); }
  x[n - 1] = S * (x[n - 1] + 2 * D * x[n - 2]);                  // last term is even
  x[n - 2] *= (-Z);                                             // scale last odd term
  }
  
// Forward DWT transform (analysis)
// n           : number of data points (even or odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
void F_CDF97_1D_split(float *x, float *e, float *o, int n){
  if(n & 1){
    F_CDF97_1D_split_N_odd(x, e, o, n);
  }else{
    F_CDF97_1D_split_N_even(x, e, o, n);
  }
}
  
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
void F_CDF97_1D_inplace(float *x, int n){
  if(n & 1){
    F_CDF97_1D_inplace_N_odd(x, n);
  }else{
    F_CDF97_1D_inplace_N_even(x, n);
  }
}

// Inverse DWT transform (synthesis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[(n+1)/2]  : odd coefficients of the transform (detail)
void I_CDF97_1D_split_N_even(float *x, float *e, float *o, int n){
  int i;
  int neven = (n+1) >> 1;
  int nodd  = neven;

  for(i = 0 ; i < neven ; i++){ x[i+i] = e[i]*Z ; x[i+i+1] = o[i] * (-S) ; }  // unscale and move to x

  x[0] = x[0] - 2 * D * x[1];
  for (i = 2; i < n; i += 2) x[i] -= D * (x[i+1] + x[i-1]);         // unupdate even terms #2

  for (i = 1; i < n - 2; i += 2) x[i] -= C * (x[i-1] + x[i+1]);     // unpredict odd terms #2
  x[n - 1] -= 2 * C * x[n - 2];

  x[0] -= 2 * B * x[1];
  for (i = 2; i < n; i += 2) x[i] -= B * (x[i+1] + x[i-1]);         // unupdate even terms #1

  for (i = 1; i < n - 2; i += 2) x[i] -= A * (x[i-1] + x[i+1]);     // unpredict odd terms #1
  x[n - 1] -= 2 * A * x[n - 2];
}

// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even)
// x[n]        : input data
void I_CDF97_1D_inplace_N_even(float *x, int n){
  int i;
  
  for(i = 0 ; i < n-1 ; i+=2 ) { x[i] *= Z ; x[i+1] *= (-S) ; }
  x[0] = x[0] - 2 * D * x[1];
  for (i = 2; i < n; i += 2) x[i] -= D * (x[i+1] + x[i-1]);         // unupdate even terms #2

  for (i = 1; i < n - 2; i += 2) x[i] -= C * (x[i-1] + x[i+1]);     // unpredict odd terms #2
  x[n - 1] -= 2 * C * x[n - 2];
  
  x[0] -= 2 * B * x[1];                                             // unupdate even terms #1
  for (i = 2; i < n; i += 2) x[i] -= B * (x[i+1] + x[i-1]);
  
  for (i = 1; i < n - 2; i += 2)  x[i] -= A * (x[i-1] + x[i+1]);    // unpredict odd terms #1
  x[n - 1] -= 2 * A * x[n - 2];
}

// Inverse DWT transform (synthesis)
// n           : number of data points (odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
void I_CDF97_1D_split_N_odd(float *x, float *e, float *o, int n){
  int i;
  int neven = (n+1) >> 1;
  int nodd  = n >> 1;

  for(i = 0 ; i < nodd ; i++){ x[i+i] = e[i] * Z ; x[i+i+1] = o[i] * (-S) ; }  // unscale and move to x
  x[n-1] = e[neven-1] * Z;

  x[0] -= 2 * D * x[1];                                             // unupdate even terms #2
  for (i = 2; i < n - 2; i += 2)  x[i] -= D * (x[i+1] + x[i-1]);
  x[n - 1] -= 2 * D * x[n - 2];

  for (i = 1; i < n - 1; i += 2)  x[i] -= C * (x[i-1] + x[i+1]);    // unpredict odd terms #2

  x[0] -= 2 * B * x[1];                                             // unupdate even terms #1
  for (i = 2; i < n - 2; i += 2) x[i] -= B * (x[i+1] + x[i-1]);
  x[n - 1] -= 2 * B * x[n - 2];

  for (i = 1; i < n - 1; i += 2) x[i] -= A * (x[i-1] + x[i+1]);     // unpredict odd terms #1  
}

// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (odd)
// x[n]        : input data
void I_CDF97_1D_inplace_N_odd(float *x, int n){
  int i;
  
  for(i = 0 ; i < n-2 ; i+=2 ) { x[i] *= Z ; x[i+1] *= (-S) ; }        // unscale odd and even terms
  x[n-1] *= Z;

  x[0] -= 2 * D * x[1];
  for (i = 2; i < n - 2; i += 2) x[i] -= D * (x[i+1] + x[i-1]);     // unupdate even terms #2
  x[n - 1] -= 2 * D * x[n - 2];

  for (i = 1; i < n - 1; i += 2)  x[i] -= C * (x[i-1] + x[i+1]);    // unpredict odd terms #2

  x[0] -= 2 * B * x[1];                                             // unupdate even terms #1
  for (i = 2; i < n - 2; i += 2) x[i] -= B * (x[i+1] + x[i-1]);
  x[n - 1] -= 2 * B * x[n - 2];

  for (i = 1; i < n - 1; i += 2) x[i] -= A * (x[i-1] + x[i+1]);     // unpredict odd terms #1
}

// Inverse DWT transform (synthesis)
// n           : number of data points (even or odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
void I_CDF97_1D_split(float *x, float *e, float *o, int n){
  if(n < 3) return;   // 3 points minimum
  if(n & 1){
    I_CDF97_1D_split_N_odd(x, e, o, n);
  }else{
    I_CDF97_1D_split_N_even(x, e, o, n);
  }
}

// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
void I_CDF97_1D_inplace(float *x, int n){
  if(n < 3) return;   // 3 points minimum
  if(n & 1){
    I_CDF97_1D_inplace_N_odd(x, n);
  }else{
    I_CDF97_1D_inplace_N_even(x, n);
  }
}

#define VCONTRIB(DEST,SCALE,SRC1,SRC2,N) { int i; for(i=0 ; i<N ; i++) {DEST[i] += SCALE *(SRC1[i] + SRC2[i]) ; } }
#define VSCALE(WHAT,FAC,N)  { int i; for(i=0 ; i<N ; i++) {WHAT[i] *= FAC ; } } 

void I_CDF97_2D_inplace(float *x, int ni, int lni, int nj){
  int neven = (nj+1) >> 1;
  int nodd = nj >> 1;
  int j;
  float *rowd, *rows1, *rows2;
  int lni2 = lni+lni;
#if ! defined(COMBINED)
  rowd = x;
  for(j = 0 ; j < nj ; j++){                 // un scaling pass
    if(j & 1){    // odd rows
      VSCALE(rowd,(-S),ni)
    }else{         // even rows
      VSCALE(rowd,(Z),ni)
    }
    rowd += lni;
  }

  rowd = x; rows1 = x + lni;                 // unupdate even rows #2
  VCONTRIB(rowd, (-D), rows1, rows1, ni);    // row 0, first even row
  rowd += lni2;
  for(j = 1 ; j < nodd ; j++){
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, (-D), rows1, rows2, ni);
    rowd += lni2;
  }
  if(neven > nodd) { rows1 = rowd - lni ; VCONTRIB(rowd, (-D), rows1, rows1, ni); }  // nj odd, last row is even

  rowd = x + lni;
  for(j = 0 ; j < neven-1 ; j++){            // un predict odd rows #2
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, (-C), rows1, rows2, ni);
    rowd += lni2;
  }
  if(nodd == neven){rows1 = rowd - lni ; VCONTRIB(rowd, (-C), rows1, rows1, ni); }  // nj even, last row is odd
#else
// combined unscale, unupdate #2, unpredict #2
  rowd = x; rows1 = x + lni;                 // unupdate even rows #2, unpredict odd rows #2
  VSCALE(rowd,(Z),ni) ;                      // unscale first even row
  VSCALE(rows1,(-S),ni) ;                    // unscale first odd row
  VCONTRIB(rowd, (-D), rows1, rows1, ni);    // row 0, first even row
  rowd += lni2;
  for(j = 1 ; j < nodd ; j++){
    VSCALE(rowd,(Z),ni) ;                    // unscale even row
    rows1 = rowd - lni; rows2 = rowd + lni;
    VSCALE(rows2,(-S),ni) ;                  // unscale odd row above
    VCONTRIB(rowd, (-D), rows1, rows2, ni);  // unupdate even rows
    rows2 = rowd - lni2;
    VCONTRIB(rows1, (-C), rowd, rows2, ni);  // unpredict odd row below 
    rowd += lni2;
  }
  if(nodd == neven){                         // last row is odd, unpredict it
    rows1 = rowd - lni; rowd =  rowd - lni2;
    VCONTRIB(rows1, (-C), rowd, rowd, ni);
  }else{                                     // last row is even, unupdate it, then unpredict odd row below
    VSCALE(rowd,(Z),ni) ;                    // unscale last even row
    rows1 = rowd - lni;                      // odd row below last even row
    VCONTRIB(rowd, (-D), rows1, rows1, ni);  // unupdate even row
    rows2 = rowd - lni2;                     // even row below last odd row
    VCONTRIB(rows1, (-C), rowd, rows2, ni);
  }
  
#endif
#if ! defined(COMBINED)
  rowd = x; rows1 = x + lni;                 // un update even rows #1
  VCONTRIB(rowd, (-B), rows1, rows1, ni);    // row 0, first even row
  rowd += lni2;
  for(j = 1 ; j < nodd ; j++){
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, (-B), rows1, rows2, ni);
    rowd += lni2;
  }
  if(neven > nodd) { rows1 = rowd - lni ; VCONTRIB(rowd, (-B), rows1, rows1, ni); }  // nj odd, last row is even

  // perform the 1D transform on the last pass after row is used (odd row unprediction)
  rowd = x + lni;
  for(j = 0 ; j < neven-1 ; j++){            // un predict odd rows #1
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, (-A), rows1, rows2, ni);
    // insert 1D in place ransform for rowd
    // I_CDF97_1D_inplace(rowd, ni);
    // insert 1D in place ransform for rows1 (previous odd row)
    // I_CDF97_1D_inplace(rows1, ni);
    rowd += lni2;
  }
  if(nodd == neven){rows1 = rowd - lni ; VCONTRIB(rowd, (-A), rows1, rows1, ni); }  // nj even, last row is odd
  // if(nodd == neven) F_CDF97_1D_inplace(rowd, ni);
  // if(nodd == neven) F_CDF97_1D_inplace(rows1, ni);
#else
  rowd = x; rows1 = x + lni;                 // un update even rows #1
  VCONTRIB(rowd, (-B), rows1, rows1, ni);    // row 0, first even row
printf("unupdate  %p\n",rowd);
  rowd += lni2;
  for(j = 1 ; j < nodd ; j++){
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, (-B), rows1, rows2, ni);  // un update even row
printf("unupdate  %p\n",rowd);
    rows2 = rowd - lni2;
    VCONTRIB(rows1, (-A), rowd, rows2, ni);  // unpredict odd row below 
printf("unpredict %p\n",rows1);
//     I_CDF97_1D_inplace(rows2, ni);
printf("unadjust  %p\n",rows2);
//     I_CDF97_1D_inplace(rows1, ni);
printf("unadjust  %p\n",rows1);
    rowd += lni2;
  }
  if(nodd == neven){                         // nj even, last row is odd
    rows1 = rowd - lni; rowd =  rowd - lni2;
    VCONTRIB(rows1, (-A), rowd, rowd, ni);   // unpredict last odd row
printf("unpredict %p\n",rows1);
  }else{                                     // last row is even, unupdate it, then unpredict odd row below
    rows2 = rowd - lni2;
    rows1 = rowd - lni;                      // odd row below last even row
    VCONTRIB(rowd, (-B), rows1, rows1, ni);  // unupdate even row
printf("unupdate  %p\n",rowd);
    VCONTRIB(rows1, (-A), rowd, rows2, ni);
printf("unpredict %p\n",rows1);
//     I_CDF97_1D_inplace(rows2, ni);
printf("unadjust  %p\n",rows2);
  }
//     I_CDF97_1D_inplace(rows1, ni);
printf("unadjust  %p\n",rows1);
//     I_CDF97_1D_inplace(rowd, ni);
printf("unadjust  %p\n",rowd);
#endif
  rowd = x;
  for(j = 0 ; j < nj ; j++){          // temporary last pass for 1D transform
      I_CDF97_1D_inplace(rowd, ni);
      row += lni;
  }
}

void F_CDF97_2D_inplace(float *x, int ni, int lni, int nj){
  int neven = (nj+1) >> 1;
  int nodd = nj >> 1;
  int j;
  float *rowd, *rows1, *rows2;
  int lni2 = lni+lni;

  rowd = x;
  for(j = 0 ; j < nj ; j++){          // temporary first pass for 1D transform
      F_CDF97_1D_inplace(rowd, ni);
      row += lni;
  }
// combined pass : 1D transform, predict #1 , update #1
#if defined(COMBINED)
  // perform the 1D transform on the first pass before row is used
  rowd = x + lni ; rows1 = rowd-lni ; rows2 = rowd+lni ; 
  // F_CDF97_1D_inplace(x, ni);
  // F_CDF97_1D_inplace(rowd, ni);
  // F_CDF97_1D_inplace(rows2, ni);
// printf("adjust^ %p\nadjust^ %p\nadjust^ %p\n",x,rowd,rows2);
  VCONTRIB(rowd, A, rows1, rows2, ni) ;                                        // predict first odd row
  VCONTRIB(x, B, rowd, rowd, ni) ;                                             // update first even row
// printf("predict %p\nupdate  %p\n",rowd,x);
  rowd += lni2 ;                                                               // next odd row
  for(j = 1 ; j < neven-1 ; j++){
    rows1 = rowd-lni ; rows2 = rowd+lni ;
    // F_CDF97_1D_inplace(rowd, ni);
    // F_CDF97_1D_inplace(rows2, ni);
// printf("adjust- %p\nadjust- %p\n",rowd,rows2);
    VCONTRIB(rowd, A, rows1, rows2, ni);                                       // predict odd row
    rows2 = rowd-lni2 ; VCONTRIB(rows1, B, rowd, rows2, ni) ;                  // update even row below odd row
// printf("predict %p\nupdate  %p\n",rowd,rows1);
    rowd += lni2 ;                                                             // next odd row
  }
  rows1 = rowd-lni ;
  if(nodd == neven){                          // nj even
//     F_CDF97_1D_inplace(rowd, ni);
// printf("adjust_ %p\n",rowd);
    VCONTRIB(rowd, A, rows1, rows1, ni);      // last row is odd, predict it
// printf("predict  %p\n",rowd);
    rows2 = rowd-lni2 ; VCONTRIB(rows1, B, rowd, rows2, ni);    // update even row below last odd row
  }else{                                      // nj odd, 
    rows1 = rowd-lni ; rowd = rowd - lni2;
    VCONTRIB(rows1, B, rowd, rowd, ni);    // last row is even, update it
// printf("update  %p\n",rows1);
  }
#else
  rowd = x + lni;
  // insert 1D in place ransform for row 0
  // F_CDF97_1D_inplace(x, ni);
  for(j = 0 ; j < neven-1 ; j++){            // predict odd rows #1
    // insert 1D in place ransform for rowd
    // F_CDF97_1D_inplace(rowd, ni);
    // insert 1D in place ransform for rows2
    // F_CDF97_1D_inplace(rows2, ni);
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, A, rows1, rows2, ni);
    rowd += lni2;
  }
  // if(nodd == neven) F_CDF97_1D_inplace(rowd, ni);
  if(nodd == neven){rows1 = rowd - lni ; VCONTRIB(rowd, A, rows1, rows1, ni); }  // nj even, last row is odd
  rowd = x; rows1 = x + lni;                 // update even rows #1
  VCONTRIB(rowd, B, rows1, rows1, ni);       // row 0, first even row
  rowd += lni2;
  for(j = 1 ; j < nodd ; j++){
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, B, rows1, rows2, ni);
    rowd += lni2;
  }
  if(neven > nodd) { rows1 = rowd - lni ; VCONTRIB(rowd, B, rows1, rows1, ni); }  // nj odd, last row is even
#endif
// combined pass : predict #2 , update #2, scaling
#if defined(COMBINED)
  rowd = x + lni ; rows1 = rowd-lni ; rows2 = rowd+lni ; 
  VCONTRIB(rowd, C, rows1, rows2, ni) ;                                        // predict first odd row
  VCONTRIB(x, D, rowd, rowd, ni) ;                                             // update first even row
  VSCALE(x,(S),ni);
  rowd += lni2 ;                                                               // next odd row
  for(j = 1 ; j < neven-1 ; j++){
    rows1 = rowd-lni ; rows2 = rowd+lni ;
    VCONTRIB(rowd, C, rows1, rows2, ni);                                       // predict odd row
    rows2 = rowd-lni2 ; VCONTRIB(rows1, D, rowd, rows2, ni) ;                  // update even row below odd row
    VSCALE(rows2,(-Z),ni);
    VSCALE(rows1,(S),ni);
    rowd += lni2 ;                                                             // next odd row
  }
  rows1 = rowd-lni ;
  if(nodd == neven){                          // nj even
    VCONTRIB(rowd, C, rows1, rows1, ni);      // last row is odd, predict it
    rows2 = rowd-lni2 ; VCONTRIB(rows1, D, rowd, rows2, ni);    // update even row below last odd row
    VSCALE(rows2,(-Z),ni);
    VSCALE(rows1,(S),ni);
    VSCALE(rowd,(-Z),ni);
  }else{                                      // nj odd, 
    rows1 = rowd-lni ; rowd = rowd - lni2;
    VCONTRIB(rows1, D, rowd, rowd, ni);    // last row is even, update it
    VSCALE(rowd,(-Z),ni);
    VSCALE(rows1,(S),ni);
  }
#else
  rowd = x + lni;
  for(j = 0 ; j < neven-1 ; j++){            // predict odd rows #2
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, C, rows1, rows2, ni);
    rowd += lni2;
  }
  if(nodd == neven){rows1 = rowd - lni ; VCONTRIB(rowd, C, rows1, rows1, ni); }  // nj even, last row is odd

  rowd = x; rows1 = x + lni;                 // update even rows #2
  VCONTRIB(rowd, D, rows1, rows1, ni);       // row 0, first even row
  rowd += lni2;
  for(j = 1 ; j < nodd ; j++){
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, D, rows1, rows2, ni);
    rowd += lni2;
  }
  if(neven > nodd) { rows1 = rowd - lni ; VCONTRIB(rowd, D, rows1, rows1, ni); }  // nj odd, last row is even

  rowd = x;
  for(j = 0 ; j < nj ; j++){                 // scaling pass
    if(j & 1){    // odd rows
      VSCALE(rowd,(-Z),ni)
    }else{         // even rows
      VSCALE(rowd,(S),ni)
    }
    rowd += lni;
  }
#endif
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NPTS 9

int main() {
  float x[NPTS+1], y[NPTS+1], e[NPTS+1], o[NPTS+1], z[NPTS+1], d[NPTS+1];
  int i;

  // Makes a fancy cubic signal
  for (i=0;i<NPTS;i++) x[i]=5+i+0.4*i*i-0.02*i*i*i;
  for (i=0;i<NPTS;i++) y[i]=x[i];
  for (i=0;i<NPTS;i++) d[i]=x[i];
  for (i=0;i<NPTS;i++) { e[i] = 0 ; o[i] = 0 ; }
  
  // Prints original sigal x
  printf("Original signal:\n");
  for (i=0;i<NPTS;i++) printf("x[%2d]=%10f\n",i,x[i]);
  printf("\n");

  // Do the forward 9/7 transform
  F_CDF97_1D_split(x,e,o,NPTS);
  I_CDF97_1D_split(z,e,o,NPTS);
  F_CDF97_1D_inplace(x,NPTS);
  F_CDF97_1D_inplace(y,NPTS);
  F_CDF97_2D_inplace(d, 1, 1, NPTS);
//   for (i=1;i<NPTS;i+=2) y[i] = (fabs(y[i]) > .01) ? y[i] : 0;
  
  // Prints the wavelet coefficients
  printf("Wavelets coefficients:\n");
  for (i=0;i<NPTS;i+=2) printf("wc[%2d,%2d]=%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n",
                               i,i+1,x[i],x[i+1],y[i],y[i+1],e[i>>1],o[i>>1],d[i],d[i+1]);
  printf("\n");

  // Do the inverse 9/7 transform
  I_CDF97_1D_inplace(x,NPTS); 
  I_CDF97_1D_inplace(y,NPTS); 
  I_CDF97_2D_inplace(d, 1, 1, NPTS);

  // Prints the reconstructed signal 
  printf("Reconstructed signal:\n");
  for (i=0;i<NPTS;i++) printf("xx[%2d]=%10f,%10f,%10f,%10f,%10f,%10f,%10f\n",i,x[i]-(5+i+0.4*i*i-0.02*i*i*i), y[i]-x[i],(5+i+0.4*i*i-0.02*i*i*i),x[i],y[i],z[i],d[i]);

//   x[31] = 0.0;
//   F_CDF97_1D_inplace(x,31);
//   I_CDF97_1D_inplace(x,31);
//   printf("Reconstructed signal:\n");
//   for (i=0;i<NPTS;i++) printf("xx[%d]=%f\n",i,x[i]-(5+i+0.4*i*i-0.02*i*i*i));

}