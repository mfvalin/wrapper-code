//  useful routines for C and FORTRAN programming
//  Copyright (C) 2022  Environnement Canada
//
//  This is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation,
//  version 2.1 of the License.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

// quick and dirty evaluations of compression losses
void AnalyzeCompressionErrors(float *fa, float *fb, int np, float small, char *str){  // will have to add a few options
  int i;
  float maxval, minval, relerr, rdiff, snr;
  double err, errmax, errsum, errsuma, sum2, acc2, acc0;
  double suma, sumb, suma2, sumb2, sumab;
  double vara, varb, varab, avga, avgb, rab;
  double ssim;  // structural similarity
  uint32_t *ia = (uint32_t *)fa;
  uint32_t *ib = (uint32_t *)fb;
  uint32_t ierr = 0;
  int32_t idiff;
  int indx, iacc, iabs;
  int accuracy, n;
  uint64_t idif64;

  small  = fabsf(small); // in case small is negative
  err    = 0.0f;
  sum2   = 0.0;      // sum of errors**2 (for RMS)
  suma2  = 0.0;      // sum of squares, array fa
  sumb2  = 0.0;      // sum of squares, array fb
  suma   = 0.0;      // sum of terms, array fa
  sumb   = 0.0;      // sum of terms, array fb
  sumab  = 0.0;      // sum of fa*fb products (for covariance)
  indx   = 0;        // position of largest relative error
  iacc   = 0;        // position of largest bit inaccuracy
  iabs   = 0;        // position of largest absolute error
  ierr   = 0;        // largest bit inaccuracy
  idif64 = 0;        // sum of bit inaccuracies
  relerr = 0.0f;     // largest relative error
  errmax = 0.0f;     // largest absolute error
  errsum = 0.0f;     // sum of errors (for BIAS)
  errsuma = 0.0f;    // sum of absolute errors
  maxval = fa[0];    // highest signed value in array fa
  minval = fa[0];    // lowest signed value in array fa
  n      = 0;
  for(i=0 ; i < np ; i++){
    suma  += fa[i] ;
    suma2 += ( fa[i] * fa[i] ) ;
    sumb  += fb[i] ;
    sumb2 += ( fb[i] * fb[i] ) ;
    sumab += ( fb[i] * fa[i] ) ;
    maxval = (fa[i] > maxval ) ? fa[i] : maxval ;
    minval = (fa[i] < minval ) ? fa[i] : minval ;
    err = fb[i] - fa[i] ;               // signed error
    sum2 = sum2 + err * err ;           // sum of squared error (to compute RMS)
    errsum += err ;                     // sum of signed errors (to compute BIAS)
    err = fabs(err) ;                   // absolute error
    errsuma += err ;                    // sum of absolute errors
    if(err > errmax) {errmax = err ; iabs = i; } ;
    if(fabsf(fa[i]) <= small) continue ;  // ignore absolute values smaller than threshold
    if( (fa[i] > 0.0f && fb[i] < 0.0f ) || (fa[i] < 0.0f && fb[i] > 0.0f) ){
      continue ;  // opposite signs, ignore
    }
    n++;
    rdiff = fabsf(fa[i] - fb[i]) ;
    rdiff = rdiff / fa[i] ;              // fa[i] should never be zero at this point
    if(rdiff > relerr) { relerr = rdiff; indx = i ; }   // largest relative error
    idiff = (ia[i] > ib[i]) ? (ia[i] - ib[i]) : (ib[i] - ia[i]) ;
    idif64 += idiff;
    if(idiff > ierr) { ierr = idiff ; iacc = i ; }
  }
//   printf("np, n = %d %d\n",np,n);
  avga  = suma/np;
  vara  = suma2/np - avga*avga ; // variance of a
  avgb  = sumb / np;
  varb  = sumb2/np - avgb*avgb ; // variance of b
  varab = sumab/np - avga*avgb ; // covariance of a and b
  ssim  = (2.0 * avga * avgb + .01) * (2.0 * varab + .03) / ((avga*avga + avgb*avgb + .01)*(vara + varb + .03)) ; // structural similarity
  rab   = (np*sumab - suma*sumb) / (sqrt(np*suma2 - suma*suma) * sqrt(np*sumb2 - sumb*sumb)) ;  // Pearson coefficient
  idif64 = idif64/n;                                  // average ULP difference
  acc2 = log2(1.0+idif64);                     // accuracy (in agreed upon bits)
  if(acc2 < 0) acc2 = 0;
  acc0 = log2(1.0+ierr);                       // worst accuracy
  if(acc0 < 0) acc0 = 0;
  sum2 = sum2 / np;                            // average quadratic error
  sum2 = sqrt(sum2);                           // RMS
  snr  = .25 * (maxval - minval) * (maxval - minval);
  snr  = 10.0 * log10(snr / sum2);                    // Peak Signal / Noise Ratio
//   if(relerr < .000001f) relerr = .000001f;
  if(relerr == 0.0) relerr = 1.0E-10;
  printf("%s epsilon = %6.2g ",str,small);
  printf("max/avg/bias/rms/rel err = (%8.6f, %8.6f, %8.6f, %8.6f, 1/%8.2g), range = %10.6f",
         errmax, errsuma/n, errsum/n, sum2, 1.0/relerr, maxval-minval);
  printf(", worst/avg accuracy = (%6.2f,%6.2f) bits, PSNR = %5.0f",24-acc0, 24-acc2, snr);  // probably not relevant for FCST verif
//   printf(", DISSIM = %12.6g, Pearson = %12.6g", 1.0 - ssim,1.0-rab);   // not very useful for packing error analysis, maybe for FCST verif ?
  printf(" [%.0f]\n",(maxval-minval)/errmax);
  accuracy = 0;
  while(ierr >>= 1) accuracy ++;
//   printf("max error  = %8.6f, bias = %10.6f, avg error = %8.6f, rms = %8.6f, npts = %d, min,max = (%10.6f %10.6f), ",
// 	errmax, errsum/(n), errsuma/(n), sqrt(sum2), n, minval, maxval);
//   printf("range/errormax = %g\n",(maxval-minval)/errmax);
//   printf("accuracy(bits) = %9d [%8d](%8.3g,%8.3g)(%9.6f)\n",accuracy,iacc, fa[iacc], fb[iacc],fb[iacc] - fa[iacc]);
//   printf("max rel error  = %9.5g [%8d](%8.3g,%8.3g)(%9.6f)\n",relerr,indx, fa[indx], fb[indx],fb[indx] - fa[indx]);
//   printf("max abs error  = %9.6f [%8d](%8.3g,%8.3g)(%9.6f)\n",errmax,iabs, fa[iabs], fb[iabs],fb[iabs] - fa[iabs]);
//   accuracy = 24;
//   while(idif64 >>= 1) accuracy --;
//   printf("average accuracy bits = %6.2f, PSNR = %f\n",acc2, snr);
}

// create a 2D float array
void FloatTestData_2d(float *z, int ni, int lni, int nj, float alpha, float delta, float noise, int kind, float *min, float *max){
  int i, j ;
  float *t = z ;
  float x, y, v, minval, maxval ;

  minval = 1.0E37 ; maxval = -minval ;
  switch(kind)
  {
    case 1 :  // distance to center (0 -> 1.4)
      for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ; y = (y - 1.0) * (y - 1.0) ;
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ; x = (x - 1.0) * (x - 1.0) ;
          v =  sqrtf(x + y) ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          v = v * sqrtf(2.0f) * 4.0 ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
      break ;
    case 2 :  // square of distance to center (0 -> 2)
      for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ; y = (y - 1.0) * (y - 1.0) ;
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ; x = (x - 1.0) * (x - 1.0) ;
          v = x + y ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          v *= 4.0f ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
      break ;
    case 3 :  // sum of sines (-2 -> 2)
      for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ; y = y * 2.0 * 3.14159 ;   // 0 -> 4 * pi
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ; x = x * 2.0 * 3.14159 ;   // 0 -> 4 * pi
          v = sinf(x) + sinf(y) ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          v *= 2.0f ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
      break ;
    case 4 :  // sum of cosines (-2 -> 2)
      for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ; y = y * 2.0 * 3.14159 ;   // 0 -> 4 * pi
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ; x = x * 2.0 * 3.14159 ;   // 0 -> 4 * pi
          v = cosf(x) + cosf(y) ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          v *= 2.0f ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
      break ;
    case 5 :  // sine * cosine (-1 -> 1)
      for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ; y = y * 2.0 * 3.14159 ;   // 0 -> 4 * pi
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ; x = x * 2.0 * 3.14159 ;   // 0 -> 4 * pi
          v = sinf(x) * cosf(y) ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          v *= 4.0f ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
      break ;
    case 6 :  // product (0 -> 4)
      for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ;
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ;
          v = x * y ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          v *= 2.0f ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
     break ;
    case 7 :  // sum (0 -> 4)
      for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ;
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ;
          v = x + y ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          v *= 2.0f ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
     break ;
    case 8 :  // sum of squares (0 -> 8)
      for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ; y = y * y ;
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ; x = x * x ;
          v = x + y ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
     break ;
    default :   // sloped plane normalized to 1.0
       for(j=0 ; j<nj ; j++){
        y = 2.0 * j / nj ;
        for(i=0 ; i<ni ; i++){
          x = 2.0 * i / ni ;
          v = x + y ;
          v *= .25 ;
          v = v * alpha + delta + noise * (drand48() - .5) ;
          v *= 8.0f ;
          t[i] = v ;
          minval = (v < minval) ? v : minval ;
          maxval = (v > maxval) ? v : maxval ;
        }
        t += lni ;
      }
     break ;
  }
  *min = minval ;
  *max = maxval ;
}
