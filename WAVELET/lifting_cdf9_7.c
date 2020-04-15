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

#define A    (-1.58615986717275f)
#define B    (-0.05297864003258f)
#define C      0.88293362717904f
#define D      0.44350482244527f
#define S      1.14960430535816f
#define Z      0.86986452237446f

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
  
  x[1] *= (-S);
  x[0] = x[0]*Z - 2 * D * x[1];                                     // unscale all, then unupdate even terms #2
  for (i = 2; i < n; i += 2){ x[i+1] *= (-S); x[i] = x[i]/S - D * (x[i+1] + x[i-1]); }

  for (i = 1; i < n - 2; i += 2) x[i] -= C * (x[i-1] + x[i+1]);     // unpredict odd terms #2
  x[n - 1] -= 2 * C * x[n - 2];
  
  x[0] -= 2 * B * x[1];                                             // unupdate even terms #1
  for (i = 2; i < n; i += 2) x[i] -= B * (x[i+1] + x[i-1]);
  x[n - 1] -= 2 * B * x[n - 2];
  
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

  x[0] -= 2 * D * x[1];                                             // unupdate even terms #2 and unscale
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
  
  for (i = 1; i < n - 1; i += 2) x[i] *= (-S);                      // unscale odd terms

  x[0] = x[0] * Z - 2 * D * x[1];                                   // unupdate even terms #2 and unscale
  for (i = 2; i < n - 2; i += 2) x[i] = x[i] * Z - D * (x[i+1] + x[i-1]);
  x[n - 1] = x[n - 1] * Z - 2 * D * x[n - 2];

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
  if(n & 1){
    I_CDF97_1D_inplace_N_odd(x, n);
  }else{
    I_CDF97_1D_inplace_N_even(x, n);
  }
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NPTS 15

int main() {
  float x[NPTS+1], y[NPTS+1], e[NPTS+1], o[NPTS+1], z[NPTS+1];
  int i;

  // Makes a fancy cubic signal
  for (i=0;i<NPTS;i++) x[i]=5+i+0.4*i*i-0.02*i*i*i;
  for (i=0;i<NPTS;i++) y[i]=x[i];
  for (i=0;i<NPTS;i++) { e[i] = 0 ; o[i] = 0 ; }
  
  // Prints original sigal x
//   printf("Original signal:\n");
//   for (i=0;i<NPTS;i++) printf("x[%2d]=%10f\n",i,x[i]);
  printf("\n");

  // Do the forward 9/7 transform
  F_CDF97_1D_split(x,e,o,NPTS);
  I_CDF97_1D_split(z,e,o,NPTS);
  F_CDF97_1D_inplace(x,NPTS);
  F_CDF97_1D_inplace(y,NPTS);
  for (i=1;i<NPTS;i+=2) y[i] = (fabs(y[i]) > .01) ? y[i] : 0;
  
  // Prints the wavelet coefficients
  printf("Wavelets coefficients:\n");
  for (i=0;i<NPTS;i+=2) printf("wc[%2d,%2d]=%10f,%10f,%10f,%10f,%10f,%10f\n",i,i+1,x[i],x[i+1],y[i],y[i+1],e[i>>1],o[i>>1]);
  printf("\n");

  // Do the inverse 9/7 transform
  I_CDF97_1D_inplace(x,NPTS); 
  I_CDF97_1D_inplace(y,NPTS); 

  // Prints the reconstructed signal 
  printf("Reconstructed signal:\n");
  for (i=0;i<NPTS;i++) printf("xx[%2d]=%10f,%10f,%10f,%10f,%10f\n",i,x[i]-(5+i+0.4*i*i-0.02*i*i*i), y[i]-x[i],(5+i+0.4*i*i-0.02*i*i*i),x[i],z[i]);

//   x[31] = 0.0;
//   F_CDF97_1D_inplace(x,31);
//   I_CDF97_1D_inplace(x,31);
//   printf("Reconstructed signal:\n");
//   for (i=0;i<NPTS;i++) printf("xx[%d]=%f\n",i,x[i]-(5+i+0.4*i*i-0.02*i*i*i));

}