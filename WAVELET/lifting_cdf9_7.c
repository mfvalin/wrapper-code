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
#include <stdio.h>

#if defined(ORIGINAL)
#define A    (-1.58615986717275f)
#define B    (-0.05297864003258f)
#define C      0.88293362717904f
#define D      0.44350482244527f
#define S      1.14960430535816f
#define Z      0.86986452237446f
#else
#define A    (-1.586134342f)
#define B    (-0.0529801185f)
#define C      0.8829110762f
#define D      0.4435068522f
#define S      1.149604398f
#define Z      0.869864452f
#endif

// interface        !InTf
//   subroutine F_CDF97_1D_split_N_even(x, e, o, n) bind(C,name='F_CDF97_1D_split_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_CDF97_1D_split_N_even                 !InTf
// end interface    !InTf

// Cohen-Daubechies-Favreau 9/7 wavelets
// lifting implementation
// in place or  even/odd split 

// Forward DWT transform (analysis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[(n+1)/2]  : odd coefficients of the transform (detail)
void F_CDF97_1D_split_N_even(float *x, float *e, float *o, int n){    // InTc
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

// interface        !InTf
//   subroutine F_CDF97_1D_inplace_N_even(x, e, o, n) bind(C,name='F_CDF97_1D_inplace_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_CDF97_1D_inplace_N_even               !InTf
// end interface    !InTf

// Forward DWT transform (analysis) (in place, x is overwritten)
// n           : number of data points (even)
// x[n]        : input data
void F_CDF97_1D_inplace_N_even(float *x, int n){    // InTc
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

// interface        !InTf
//   subroutine F_CDF97_1D_split_N_odd(x, e, o, n) bind(C,name='F_CDF97_1D_split_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_CDF97_1D_split_N_odd                  !InTf
// end interface    !InTf

// Forward DWT transform (analysis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
void F_CDF97_1D_split_N_odd(float *x, float *e, float *o, int n){    // InTc
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

// interface        !InTf
//   subroutine F_CDF97_1D_inplace_N_odd(x, e, o, n) bind(C,name='F_CDF97_1D_inplace_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_CDF97_1D_inplace_N_odd                !InTf
// end interface    !InTf

// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (odd)
// x[n]        : input data
void F_CDF97_1D_inplace_N_odd(float *x, int n){    // InTc
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
  
// interface        !InTf
//   subroutine F_CDF97_1D_split(x, e, o, n) bind(C,name='F_CDF97_1D_split')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(IN) :: x         !InTf
//     real(C_FLOAT), dimension(*), intent(OUT) :: e, o     !InTf
//   end subroutine F_CDF97_1D_split                        !InTf
// end interface    !InTf

// Forward DWT transform (analysis)
// n           : number of data points (even or odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
void F_CDF97_1D_split(float *x, float *e, float *o, int n){    // InTc
  if(n < 3) {
    if(n > 0) e[0] = x[0];
    if(n > 1) o[0] = x[1];
    return;
  }
  if(n & 1){
    F_CDF97_1D_split_N_odd(x, e, o, n);
  }else{
    F_CDF97_1D_split_N_even(x, e, o, n);
  }
}
  
// interface        !InTf
//   subroutine F_CDF97_1D_inplace(x, e, o, n) bind(C,name='F_CDF97_1D_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_CDF97_1D_inplace                      !InTf
// end interface    !InTf
  
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
void F_CDF97_1D_inplace(float *x, int n){    // InTc
  if(n < 3) return;
  if(n & 1){
    F_CDF97_1D_inplace_N_odd(x, n);
  }else{
    F_CDF97_1D_inplace_N_even(x, n);
  }
}
  
// interface        !InTf
//   subroutine F_CDF97_1D_split_inplace(x, e, o, n) bind(C,name='F_CDF97_1D_split_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine F_CDF97_1D_split_inplace                      !InTf
// end interface    !InTf
  
// Forward DWT transform  (analysis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
void F_CDF97_1D_split_inplace(float *x, int n){    // InTc
  int i;
  int neven = (n+1) >> 1;
  int lasteven = neven - 1;
  int nodd  = n >> 1;
  int lastodd = nodd - 1;
  float odd[n];

  if(n < 3) return;  
  for (i = 0; i < nodd; i ++) { odd[i] = x[i+i+1] ; x[i] = x[i+i] ; }  // copy odd terms to temporary odd array
  x[lasteven] = x[lasteven+lasteven]; x[neven] = x[lasteven] ;         // even terms to beginning of array

  for (i = 0; i < nodd; i ++) odd[i] += A * (x[i] + x[i+1]) ;          // predict odd terms #1
  odd[nodd] = odd[lastodd] ;                                           // pad end of odd to make it at least neven terms

  x[0] += 2 * B * odd[0] ;                                             // update even terms #1
  for (i = 1; i < neven; i ++) x[i] += B * (odd[i] + odd[i-1]);
  x[neven] = x[lasteven] ;                                             // pad end of even to make it at least nodd + 1 terms

  for (i = 0; i < nodd; i ++) odd[i] += C * (x[i] + x[i+1]);           // predict odd terms #2
  odd[nodd] = odd[lastodd] ;                                           // pad end of odd to make it at least neven terms

  x[0] = S * (x[0] + 2 * D * odd[0]) ;                                 // update even terms #2 and scale them
  for (i = 1; i < neven; i ++) x[i] = S * ( x[i] + D * (odd[i] + odd[i-1]));
  for (i = 0; i < nodd; i ++) x[neven+i] = (-Z) * odd[i] ;             // scale odd terms
}

// interface        !InTf
//   subroutine I_CDF97_1D_split_N_even(x, e, o, n) bind(C,name='I_CDF97_1D_split_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_CDF97_1D_split_N_even                 !InTf
// end interface    !InTf

// Inverse DWT transform (synthesis)
// n           : number of data points (even)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[(n+1)/2]  : odd coefficients of the transform (detail)
void I_CDF97_1D_split_N_even(float *x, float *e, float *o, int n){    // InTc
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

// interface        !InTf
//   subroutine I_CDF97_1D_inplace_N_even(x, e, o, n) bind(C,name='I_CDF97_1D_inplace_N_even')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_CDF97_1D_inplace_N_even               !InTf
// end interface    !InTf

// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even)
// x[n]        : input data
void I_CDF97_1D_inplace_N_even(float *x, int n){    // InTc
  int i;
  
  for(i = 0 ; i < n-1 ; i+=2 ) { x[i] *= Z ; x[i+1] *= (-S) ; }     // unscale
  x[0] = x[0] - 2 * D * x[1];
  for (i = 2; i < n; i += 2) x[i] -= D * (x[i+1] + x[i-1]);         // unupdate even terms #2

  for (i = 1; i < n - 2; i += 2) x[i] -= C * (x[i-1] + x[i+1]);     // unpredict odd terms #2
  x[n - 1] -= 2 * C * x[n - 2];
// printf("Ieven1: "); for(i=0 ; i<n ; i+=2) printf("%8.3f ",x[i]); printf("\n");
// printf("Iodd1 : "); for(i=1 ; i<n ; i+=2) printf("%8.3f ",x[i]); printf("\n");
  
  x[0] -= 2 * B * x[1];                                             // unupdate even terms #1
  for (i = 2; i < n; i += 2) x[i] -= B * (x[i+1] + x[i-1]);
  
  for (i = 1; i < n - 2; i += 2)  x[i] -= A * (x[i-1] + x[i+1]);    // unpredict odd terms #1
  x[n - 1] -= 2 * A * x[n - 2];
// printf("Ieven2: "); for(i=0 ; i<n ; i+=2) printf("%8.3f ",x[i]); printf("\n");
// printf("Iodd2 : "); for(i=1 ; i<n ; i+=2) printf("%8.3f ",x[i]); printf("\n");
}

// interface        !InTf
//   subroutine I_CDF97_1D_split_N_odd(x, e, o, n) bind(C,name='I_CDF97_1D_split_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_CDF97_1D_split_N_odd                  !InTf
// end interface    !InTf

// Inverse DWT transform (synthesis)
// n           : number of data points (odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
void I_CDF97_1D_split_N_odd(float *x, float *e, float *o, int n){    // InTc
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

// interface        !InTf
//   subroutine I_CDF97_1D_inplace_N_odd(x, e, o, n) bind(C,name='I_CDF97_1D_inplace_N_odd')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_CDF97_1D_inplace_N_odd                !InTf
// end interface    !InTf

// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (odd)
// x[n]        : input data
void I_CDF97_1D_inplace_N_odd(float *x, int n){    // InTc
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

// interface        !InTf
//   subroutine I_CDF97_1D_split(x, e, o, n) bind(C,name='I_CDF97_1D_split')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(OUT) :: x        !InTf
//     real(C_FLOAT), dimension(*), intent(IN) :: e, o      !InTf
//   end subroutine I_CDF97_1D_split                        !InTf
// end interface    !InTf

// Inverse DWT transform (synthesis)
// n           : number of data points (even or odd)
// x[n]        : input data
// e[(n+1)/2]  : even coefficients of the transform (approximation)
// o[n/2]      : odd coefficients of the transform (detail)
void I_CDF97_1D_split(float *x, float *e, float *o, int n){    // InTc
  if(n < 3) {   // 3 points minimum
    if(n > 0) x[0] = e[0] ;
    if(n > 1) x[1] = o[0] ;
    return;
  }
  if(n & 1){
    I_CDF97_1D_split_N_odd(x, e, o, n);
  }else{
    I_CDF97_1D_split_N_even(x, e, o, n);
  }
}

// interface        !InTf
//   subroutine I_CDF97_1D_split_inplace(x, e, o, n) bind(C,name='I_CDF97_1D_split_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_CDF97_1D_split_inplace                      !InTf
// end interface    !InTf

// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
void I_CDF97_1D_split_inplace(float *x, int n){    // InTc
  int i;
  int neven = (n+1) >> 1;
  int nodd  = n >> 1;
  float even[n];
  float *odd = x + neven;

  if(n < 3) return;   // 3 points minimum
  for(i = 0 ; i < neven ; i++) even[i] = x[i] ;             // copy even terms to temporary array
  for(i = 0 ; i < nodd ; i++) x[i+i+1] = odd[i] ;           // unshuffle odd terms
  for(i = 0 ; i < neven ; i++) x[i+i]  = even[i] ;          // unshuffle even terms
  if(n & 1){
    I_CDF97_1D_inplace_N_odd(x, n);
  }else{
    I_CDF97_1D_inplace_N_even(x, n);
  }
}

// interface        !InTf
//   subroutine I_CDF97_1D_inplace(x, e, o, n) bind(C,name='I_CDF97_1D_inplace')    !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: n                      !InTf
//     real(C_FLOAT), dimension(n), intent(INOUT) :: x      !InTf
//   end subroutine I_CDF97_1D_inplace                      !InTf
// end interface    !InTf

// Inverse DWT transform (synthesis) (in place, x is overwritten)
// n           : number of data points (even or odd)
// x[n]        : input data
void I_CDF97_1D_inplace(float *x, int n){    // InTc
  if(n < 3) return;   // 3 points minimum
  if(n & 1){
    I_CDF97_1D_inplace_N_odd(x, n);
  }else{
    I_CDF97_1D_inplace_N_even(x, n);
  }
}

#define VCONTRIB(DEST,SCALE,SRC1,SRC2,N) { int i; for(i=0 ; i<N ; i++) {DEST[i] += SCALE *(SRC1[i] + SRC2[i]) ; } }
#define VCONTRIB2(WHERE,DEST,SCALE,SRC1,SRC2,N) { int i; for(i=0 ; i<N ; i++) {WHERE[i] = DEST[i] + SCALE *(SRC1[i] + SRC2[i]) ; } }
#define VSCALE(WHAT,FAC,N)  { int i; for(i=0 ; i<N ; i++) {WHAT[i] *= FAC ; } } 
#define VSCALE2(WHERE,WHAT,FAC,N)  { int i; for(i=0 ; i<N ; i++) {WHERE[i] = WHAT[i] *FAC ; } } 

// interface        !InTf
// subroutine I_CDF97_2D_inplace(x, ni, lni, nj) BIND(C,name='I_CDF97_2D_inplace') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni            !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine I_CDF97_2D_inplace                        !InTf
// end interface    !InTf

void I_CDF97_2D_inplace(float *x, int ni, int lni, int nj){    // InTc
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
    // I_CDF97_1D_split_inplace(rowd, ni);
    // insert 1D in place ransform for rows1 (previous odd row)
    // I_CDF97_1D_split_inplace(rows1, ni);
    rowd += lni2;
  }
  if(nodd == neven){rows1 = rowd - lni ; VCONTRIB(rowd, (-A), rows1, rows1, ni); }  // nj even, last row is odd
  // if(nodd == neven) F_CDF97_1D_inplace(rowd, ni);
  // if(nodd == neven) F_CDF97_1D_inplace(rows1, ni);
#else
  rowd = x; rows1 = x + lni;                 // un update even rows #1
  VCONTRIB(rowd, (-B), rows1, rows1, ni);    // row 0, first even row
// printf("unupdate  %p\n",rowd);
  rowd += lni2;
  for(j = 1 ; j < nodd ; j++){
    rows1 = rowd - lni; rows2 = rowd + lni;
    VCONTRIB(rowd, (-B), rows1, rows2, ni);  // un update even row
// printf("unupdate  %p\n",rowd);
    rows2 = rowd - lni2;
    VCONTRIB(rows1, (-A), rowd, rows2, ni);  // unpredict odd row below 
// printf("unpredict %p\n",rows1);
//     I_CDF97_1D_split_inplace(rows2, ni);
// printf("unadjust  %p\n",rows2);
//     I_CDF97_1D_split_inplace(rows1, ni);
// printf("unadjust  %p\n",rows1);
    rowd += lni2;
  }
  if(nodd == neven){                         // nj even, last row is odd
    rows1 = rowd - lni; rowd =  rowd - lni2;
    VCONTRIB(rows1, (-A), rowd, rowd, ni);   // unpredict last odd row
// printf("unpredict %p\n",rows1);
  }else{                                     // last row is even, unupdate it, then unpredict odd row below
    rows2 = rowd - lni2;
    rows1 = rowd - lni;                      // odd row below last even row
    VCONTRIB(rowd, (-B), rows1, rows1, ni);  // unupdate even row
// printf("unupdate  %p\n",rowd);
    VCONTRIB(rows1, (-A), rowd, rows2, ni);
// printf("unpredict %p\n",rows1);
//     I_CDF97_1D_split_inplace(rows2, ni);
// printf("unadjust  %p\n",rows2);
  }
//     I_CDF97_1D_split_inplace(rows1, ni);
// printf("unadjust  %p\n",rows1);
//     I_CDF97_1D_split_inplace(rowd, ni);
// printf("unadjust  %p\n",rowd);
#endif
  rowd = x;
  for(j = 0 ; j < nj ; j++){          // temporary last pass for 1D transform
      I_CDF97_1D_split_inplace(rowd, ni);
      rowd += lni;
  }
}

// interface        !InTf
// subroutine I_CDF97_2D_split_inplace(x, ni, lni, nj) BIND(C,name='I_CDF97_2D_split_inplace') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni            !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine I_CDF97_2D_split_inplace                        !InTf
// end interface    !InTf

void I_CDF97_2D_split_inplace(float *x, int ni, int lni, int nj){    // InTc
  int neven = (nj+1) >> 1;
  int ieven = (ni+1) >> 1;
  int nodd = nj >> 1;
  int j;
  float *rowd, *row0, *rows, *row1, *row2;
  int lni2 = lni+lni;
  float t[(nj*ni) >> 1];
  float *rowt, *oddrow;;
  // combined unscale, unupdate #2, unpredict #2
  rows = x + neven * lni ;                  // first odd row in x
  rowt = t ;                                // first temporary odd row
  row1 = x ; row2 = row1 ;                  // first even row
  VSCALE(row1,(Z),ni) ;                     // unscale first even row
  VSCALE2(rowt,rows,(-S),ni) ;              // unscale first odd row
  VCONTRIB(row1, (-D), rowt, rowt, ni);     // unupdate #2 first even row
  for(j = 1 ; j < nodd ; j++){
    row1 = row2 ; row2 = row1 + lni ;       // next even row pair in x
    VSCALE(row2,(Z),ni) ;                   // unscale next even row
    rows += lni ; rowt += ni ;              // next odd row (x and temporary)
    row0 = rowt - ni;                       // previous (unscaled) odd row
    VSCALE2(rowt,rows,(-S),ni) ;            // unscale next odd row
    VCONTRIB(row2, (-D), row0, rowt, ni);   // unupdate next even row
    VCONTRIB(row0, (-C), row1, row2, ni);   // unpredict odd row
  }
  if(nodd == neven){                         // last row is odd, unpredict it
    VCONTRIB(rowt, (-C), row2, row2, ni);
  }else{                                     // last row is even, unupdate it, then unpredict odd row below
    row1 = row2 ; row2 = row1 + lni ;        // next even row pair in x
    VSCALE(row2,(Z),ni) ;                    // unscale last even row
    VCONTRIB(row2, (-D), rowt, rowt, ni);    // unupdate last even row
    VCONTRIB(rowt, (-C), row1, row2, ni);    // unpredict last odd row
  }
  // combined unupdate #1, unpredict #1
  rowt = t ;                                 // first temporary odd row
  row1 = x ; row2 = row1 ;                   // first even row
  VCONTRIB(row1, (-B), rowt, rowt, ni);      // unupdate #2 first even row
  for(j = 1 ; j < nodd ; j++){
    row1 = row2 ; row2 = row1 + lni ;        // next even row pair in x
    rowt += ni ; row0 = rowt - ni;           // current and previous odd rows
    VCONTRIB(row2, (-B), row0, rowt, ni);    // un update even row
    VCONTRIB(row0, (-A), row1, row2, ni);    // unpredict odd row below 
  }
  if(nodd == neven){                         // nj even, last row is odd
    VCONTRIB(rowt, (-A), row2, row2, ni);    // unpredict last odd row
  }else{                                     // last row is even, unupdate it, then unpredict odd row below
    row1 = row2 ; row2 = row1 + lni ;        // next even row pair in x
    VCONTRIB(row2, (-B), rowt, rowt, ni);    // unupdate even row
    VCONTRIB(rowt, (-A), row1, row2, ni);
  }
  // last pass, 1D split inverse transform (to be eventually combined wit unupdate/unpredict #1)
  rows = x + (nj -1 ) * lni ;                // last row in x
  rowt = t + (nodd - 1) * ni ;               // last odd row in t
  rowd = x + (neven - 1) * lni ;             // last even row in split x
  for(j = nj - 1 ; j > 0 ; j--){             // last pass for 1D transform
    if(j & 1){                               // odd row, from rowt to rows
      I_CDF97_1D_split(rows, rowt, rowt + ieven, ni) ;
      rowt -= ni ;
    }else{                                   // even row, from rowd to rows
      I_CDF97_1D_split(rows, rowd, rowd + ieven, ni);
      rowd -= lni;
    }
    rows -= lni;
  }
  I_CDF97_1D_split_inplace(x, ni);           // first even row is in place
}

// interface        !InTf
// subroutine F_CDF97_2D_inplace(x, ni, lni, nj) BIND(C,name='F_CDF97_2D_inplace') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni            !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine F_CDF97_2D_inplace                        !InTf
// end interface    !InTf

void F_CDF97_2D_inplace(float *x, int ni, int lni, int nj){    // InTc
  int neven = (nj+1) >> 1;
  int nodd = nj >> 1;
  int j;
  float *rowd, *rows1, *rows2;
  int lni2 = lni+lni;

  rowd = x;
  for(j = 0 ; j < nj ; j++){          // temporary first pass for 1D transform
      F_CDF97_1D_split_inplace(rowd, ni);
      rowd += lni;
  }
// combined pass : 1D transform, predict #1 , update #1
#if defined(COMBINED)
  // perform the 1D transform on the first pass before row is used
  rowd = x + lni ; rows1 = rowd-lni ; rows2 = rowd+lni ; 
  // F_CDF97_1D_split_inplace(x, ni);
  // F_CDF97_1D_split_inplace(rowd, ni);
  // F_CDF97_1D_split_inplace(rows2, ni);
// printf("adjust^ %p\nadjust^ %p\nadjust^ %p\n",x,rowd,rows2);
  VCONTRIB(rowd, A, rows1, rows2, ni) ;                                        // predict first odd row
  VCONTRIB(x, B, rowd, rowd, ni) ;                                             // update first even row
// printf("predict %p\nupdate  %p\n",rowd,x);
  rowd += lni2 ;                                                               // next odd row
  for(j = 1 ; j < neven-1 ; j++){
    rows1 = rowd-lni ; rows2 = rowd+lni ;
    // F_CDF97_1D_split_inplace(rowd, ni);
    // F_CDF97_1D_split_inplace(rows2, ni);
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
  // F_CDF97_1D_split_inplace(x, ni);
  for(j = 0 ; j < neven-1 ; j++){            // predict odd rows #1
    // insert 1D in place ransform for rowd
    // F_CDF97_1D_split_inplace(rowd, ni);
    // insert 1D in place ransform for rows2
    // F_CDF97_1D_split_inplace(rows2, ni);
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

// interface        !InTf
// subroutine F_CDF97_2D_split_inplace(x, ni, lni, nj) BIND(C,name='F_CDF97_2D_split_inplace') !InTf
//     import :: C_FLOAT, C_INT                             !InTf
//     integer, intent(IN), value :: ni, nj, lni            !InTf
//     real(C_FLOAT), dimension(*), intent(INOUT) :: x      !InTf
// end subroutine F_CDF97_2D_split_inplace                  !InTf
// end interface    !InTf

void F_CDF97_2D_split_inplace(float *x, int ni, int lni, int nj){    // InTc
  int neven = (nj+1) >> 1;
  int ieven = (ni+1) >> 1;
  int nodd = nj >> 1;
  int i, j;
  float *rowd, *row0, *row1, *row2;
  int lni2 = lni+lni;
  float t[(nj*ni) >> 1];
  float *rowt, *rowx, *rows, *top ;

  top = x + nj * lni ;                // top of array (one row above last, rows must be < top)
  rowx = x;                           // first even row storage
  rowt = t;                           // first odd row temporary storage
  // 1 D transform for the first 4 rows
  rows = x; F_CDF97_1D_split_inplace(rowx, ni); rowx += lni ;  // split, inplace 1D transform (first even row)
  rows += lni ; if(rows < top) F_CDF97_1D_split(rows, rowt, rowt + ieven, ni) ; rowt += ni ;  // first odd row
  rows += lni ; if(rows < top) F_CDF97_1D_split(rows, rowx, rowx + ieven, ni) ; rowx += lni ; // next even row
  rows += lni ; if(rows < top) F_CDF97_1D_split(rows, rowt, rowt + ieven, ni) ; rowt += ni ;  // next odd row
  // combined pass : predict #1 , update #1
  // perform the 1D transform on on the fly in the first pass before row is needed
  rowd = t ; row1 = x ; row2 = row1 + lni;
  VCONTRIB(rowd, A, row1, row2, ni) ;         // predict first odd row
  VCONTRIB(row1, B, rowd, rowd, ni) ;         // update first even row
  for(j = 1 ; j < neven-1 ; j++){
    rows += lni ; if(rows < top) F_CDF97_1D_split(rows, rowx, rowx + ieven, ni) ; rowx += lni ; // next even row
    rows += lni ; if(rows < top) F_CDF97_1D_split(rows, rowt, rowt + ieven, ni) ; rowt += ni ;  // next odd row
    rowd += ni ;                              // next odd row
    row1 = row2 ; row2 = row1 + lni ;         // next pair of even rows
    VCONTRIB(rowd, A, row1, row2, ni) ;       // predict odd row
    row0 = rowd - ni ; 
    VCONTRIB(row1, B, rowd, row0, ni) ;       // update even row below odd row
  }
  if(nodd == neven){                          // nj even
    rowd += ni ;                              // next (last) odd row
    VCONTRIB(rowd, A, row2, row2, ni);        // last row is odd, predict it using last even row
    row0 = rowd-ni ; 
    VCONTRIB(row2, B, rowd, row0, ni);        // update even row below last odd row
  }else{                                      // nj odd, 
    VCONTRIB(row2, B, rowd, rowd, ni);        // last row is even, update it using last odd row
  }
// combined pass : predict #2 , update #2, scaling
  rowd = t ; row1 = x ; row2 = row1 + lni;
  VCONTRIB(rowd, C, row1, row2, ni) ;         // predict first odd row
  VCONTRIB(row1, D, rowd, rowd, ni) ;         // update first even row
  VSCALE(row1,(S),ni);                        // scale first even row
  rows = x + neven * lni ;                    // first odd row in split x array 
  for(j = 1 ; j < neven-1 ; j++){
    rowd += ni ;                              // next odd row
    row1 = row2 ; row2 = row1 + lni ;         // next pair of even rows
    VCONTRIB(rowd, C, row1, row2, ni) ;       // predict odd row
    row0 = rowd - ni ; 
    VCONTRIB(row1, D, rowd, row0, ni) ;       // update even row below odd row
    VSCALE2(rows,row0,(-Z),ni);               // scale low odd row and transfer it to upper part of x
    VSCALE(row1,(S),ni) ;                     // scale lower even row
    rows += lni ;
  }
  if(nodd == neven){                          // nj even
    rowd += ni ;                              // next odd row
    VCONTRIB(rowd, C, row2, row2, ni);        // last row is odd, predict it using last even row
    row0 = rowd - ni ; 
    VCONTRIB(row2, D, rowd, row0, ni);        // update even row below last odd row
    VSCALE(row2,(S),ni);                      // scale last even row
    VSCALE2(rows,row0,(-Z),ni);               // scale last 2 odd rows
    rows += lni ;
    VSCALE2(rows,rowd,(-Z),ni);
  }else{                                      // nj odd, 
    VCONTRIB(row2, D, rowd, rowd, ni);        // last row is even, update it
    VSCALE2(rows,rowd,(-Z),ni);               // scale last odd row
    VSCALE(row2,(S),ni);                      // scale last even row
  }
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NPTS 16

int main() {
  float x[NPTS+1], y[NPTS+1], e[NPTS+1], o[NPTS+1], z[NPTS+1], d[NPTS+1];
  float xy[NPTS][NPTS] ;
  int i, j, k;
  double sum2;
  float quantum;
  int npts2 = (NPTS+1)/2;
  int npts4 = (npts2+1)/2;

  // Makes a fancy cubic signal
  for (i=0;i<NPTS;i++) x[i]=5+i+0.4*i*i-0.02*i*i*i;
  for (i=0;i<NPTS;i++) y[i]=x[i];
  for (i=0;i<NPTS;i++) d[i]=x[i];
  for (i=0;i<NPTS;i++) { e[i] = 0 ; o[i] = 0 ; }
  for (j=0;j<NPTS;j++) {
    for (i=0;i<NPTS;i++) xy[j][i] = (3+i+0.4*i*i-0.02*i*i*i) * (3+i+0.4*j*j-0.02*j*j*j);
//     for (i=0;i<NPTS;i++) xy[j][i] = sqrt((i-7.45)*(i-7.45) + (j-7.55)*(j-7.55));
  }
  
  printf("Original 2D signal:\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }
// void F_CDF97_2D_split_inplace(float *x, int ni, int lni, int nj)
  F_CDF97_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
  F_CDF97_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
  F_CDF97_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
//   F_CDF97_2D_inplace((float *)xy, NPTS, NPTS, NPTS);
//   for (j=0;j<NPTS;j++) {
//     for (i=0;i<NPTS;i++) xy[j][i] = 0;
//   }
//   xy[3*NPTS/4][3*NPTS/4] = 1.0;
  
  printf("Transformed 2D signal (before quantification):\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }
  quantum = .025;
  printf("quantum used = %8.3f\n",quantum);
  for (j=0;j<NPTS;j++) {               // quantification pass
    for (i=0;i<NPTS;i++) { k = xy[j][i] / quantum + .5f ; xy[j][i] = k * quantum ; }
  }
//   printf("Quantified 2D signal:\n");
//   for (j=NPTS-1;j>=0;j--) {
//     for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
//     printf("\n");
//   }
//   I_CDF97_2D_inplace((float *)xy, NPTS, NPTS, NPTS);
//   sum2 = 0.0;
//   for (j=0;j<NPTS;j++) {
//     for (i=0;i<NPTS;i++) sum2 += xy[j][i]*xy[j][i];
//   }
  I_CDF97_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
  I_CDF97_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
  I_CDF97_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
//   printf("Restored 2D signal:\n");
//   for (j=NPTS-1;j>=0;j--) {
//     for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
//     printf("\n");
//   }
  printf("Restored 2D signal error:\n");
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i] - ((3+i+0.4*i*i-0.02*i*i*i) * (3+i+0.4*j*j-0.02*j*j*j)));
    printf("\n");
  }
exit(0);
//   printf(" ddot, sqrt(ddot) = %f %f\n",sum2,sqrt(sum2));
  // Prints original sigal x
  printf("\nOriginal 1D signal:\n");
  for (i=0;i<NPTS;i++) printf("%8.3f ",x[i]); printf("\n");
  printf("Transformed 1D signal:\n");
  F_CDF97_1D_inplace(y,NPTS);
  for (i=0;i<NPTS;i++) printf("%8.3f ",y[i]); printf("\n");
  printf("Transformed shuffled 1D signal:\n");
  F_CDF97_1D_split_inplace(x,NPTS);
  for (i=0;i<NPTS;i++) printf("%8.3f ",x[i]); printf("\n");
  printf("Transformed shuffled 2D signal:\n");
  F_CDF97_2D_split_inplace(d, 1, 1, NPTS);
  for (i=0;i<NPTS;i++) printf("%8.3f ",d[i]); printf("\n");
  printf("Restored shuffled 1D signal:\n");
  I_CDF97_1D_split_inplace(x,NPTS);
  for (i=0;i<NPTS;i++) printf("%8.3f ",x[i]); printf("\n");
  printf("Restored shuffled 2D signal:\n");
  I_CDF97_2D_split_inplace(d, 1, 1, NPTS);
  for (i=0;i<NPTS;i++) printf("%8.3f ",d[i]); printf("\n");

exit(0);
  // Do the forward 9/7 transform
  F_CDF97_1D_split(x,e,o,NPTS);
  I_CDF97_1D_split(z,e,o,NPTS);
  F_CDF97_1D_inplace(x,NPTS);
  F_CDF97_1D_inplace(y,NPTS);
  F_CDF97_2D_inplace(d, 1, 1, NPTS);
//   for (i=1;i<NPTS;i+=2) y[i] = (fabs(y[i]) > .01) ? y[i] : 0;
  
  // Prints the wavelet coefficients
//   printf("Wavelets coefficients:\n");
//   for (i=0;i<NPTS;i+=2) printf("wc[%2d,%2d]=%10f,%10f,%10f,%10f,%10f,%10f,%10f,%10f\n",
//                                i,i+1,x[i],x[i+1],y[i],y[i+1],e[i>>1],o[i>>1],d[i],d[i+1]);
//   printf("\n");
for (i=0 ; i<NPTS ; i++) d[i] = 0 ; d[NPTS/4] = 1.0 ;
  // Do the inverse 9/7 transform
  I_CDF97_1D_inplace(x,NPTS); 
  I_CDF97_1D_inplace(y,NPTS); 
  I_CDF97_2D_inplace(d, 1, 1, NPTS);

  // Prints the reconstructed signal 
  printf("Reconstructed signal:\n");
  for (i=0;i<NPTS;i++) printf("xx[%2d]=%10f,%10f,%10f,%10f,%10f,%10f,%10f\n",i,x[i]-(5+i+0.4*i*i-0.02*i*i*i), y[i]-x[i],(5+i+0.4*i*i-0.02*i*i*i),x[i],y[i],z[i],d[i]);
sum2 = 0.0; for (i=0 ; i<NPTS ; i++) sum2 = sum2 + d[i]*d[i] ; printf(" d.d, sqrt(d.d) = %f %f\n",sum2,sqrt(sum2));

//   x[31] = 0.0;
//   F_CDF97_1D_inplace(x,31);
//   I_CDF97_1D_inplace(x,31);
//   printf("Reconstructed signal:\n");
//   for (i=0;i<NPTS;i++) printf("xx[%d]=%f\n",i,x[i]-(5+i+0.4*i*i-0.02*i*i*i));

}