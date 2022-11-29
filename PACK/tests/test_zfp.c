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
#include <string.h>
#include <math.h>
#include <misc_timers.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include <misc_zfp.h>

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
    if(fa[i] < 0.0f || fb[i] < 0.0f ) continue ;  // opposite signs, ignore
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
  sum2 = sum2 / np;                                    // average quadratic error
  snr  = .25 * (maxval - minval) * (maxval - minval);
  snr  = 10.0 * log10(snr / sum2);                    // Peak Signal / Noise Ratio
//   if(relerr < .000001f) relerr = .000001f;
  if(relerr == 0.0) relerr = 1.0E-10;
  printf("%s[%6.4f] ",str,small);
  printf("max/avg/bias/rms/rel err = (%8.6f, %8.6f, %8.6f, %8.6f, 1/%8.2g), range = %10.6f",errmax, errsuma/n, errsum/n, sqrt(sum2),1.0/relerr, maxval-minval);
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

int main(int argc, char **argv){
  /* allocate array of floats */
  int nx = 128; int nx2 = nx/2 ;
  int ny = 128; int ny2 = ny/2 ;
  int nz = 1; int nz2 = 4 ;
  int nxy2 = nx2 * ny2;
  float* array = malloc(nx * ny * nz * sizeof(float));
  float* Array = malloc(nx * ny * nz * sizeof(float));
  float* brray = malloc(nx * ny * nz * sizeof(float));
  float maxval, minval ;
  int i, j, k;
  float toler = 0.0f, maxerror = 0.0f ;
  int precision = 0 ;
  void *Stream1, *Stream2;
  int Ssize, Ssize2, errors ;
  uint64_t t ;
  int d[4], s[4], nsize ;
  double t_ns ;

  ZfpCompress_set_debug(1) ;
  /* initialize array to be compressed */
  minval = 1.0E37 ; maxval = -minval ;
  // array is a 2D (nx, ny) array
  for (k = 0; k < nz; k++){
    for (j = 0; j < ny; j++){
      for (i = 0; i < nx; i++) {
        float x = 2.0 * i / nx;
        float y = 2.0 * j / ny;
        float z = 2.0 * k / nz;
        float v = (x * x + y * y + z * z) + .002f ;
        array[i + nx * (j + ny * k)] = v ;
        minval = (v < minval) ? v : minval ;
        maxval = (v > maxval) ? v : maxval ;
//         array[i + nx * (j + ny * k)] = exp(-(x * x + y * y + z * z));
      }
    }
  }
  // disguise array as a 3D (nx/2, ny/2, 4) array
  for (j = 0; j < ny2; j++){
    for (i = 0; i < nx2; i++){
      Array[i + j*nx2 + 0*nxy2] = array[2*i + 2*j*nx + 1] ;
      Array[i + j*nx2 + 1*nxy2] = array[2*i + 2*j*nx + 0] ;
      Array[i + j*nx2 + 2*nxy2] = array[2*i + 2*j*nx + nx] ;
      Array[i + j*nx2 + 3*nxy2] = array[2*i + 2*j*nx + nx+1] ;
    }
  }

  if(argc > 1) sscanf(argv[1],"%f",&toler);                                          // get precision control
  if(toler > 0) maxerror = toler ;
  else          precision = -toler ;
  fprintf(stderr,"zfp_codec_version = %u, ZFP_CODEC = %u\n", zfp_codec_version, ZFP_CODEC) ;
  fprintf(stderr,"precision = %d, max tolerance = %f, min = %f, max = %f\n", precision, maxerror, minval, maxval) ;
  fprintf(stderr,"\n") ;

  fprintf(stderr,"============= compress/decompress 2D array =============\n");
  t  = elapsed_cycles() ;
  Stream1 = ZfpCompress(array, nx,   ny,  nz, maxerror, precision , &Ssize) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  nsize = ZfpArrayDims(Stream1, d, s, Ssize) ;
  memset(brray, 0, sizeof(float)*nx*ny*nz) ;
  t  = elapsed_cycles() ;
  ZfpExpand(brray, nx, ny, nz, Stream1, Ssize) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  AnalyzeCompressionErrors(array, brray, nx*ny*nz, 0.0f, "Ctest1");     // evaluate errors
  fprintf(stderr,"\n") ;

  // Array is array disguised as a 3D array, 4 levels
  // level 0  odd i, even j elements
  // level 1 even i, even j elements
  // level 2 even i,  odd j elements
  // level 3  odd i,  odd j elements
  fprintf(stderr,"========== compress/decompress 2D as 3D array ============\n");
  if(precision >0) precision += 1 ;
  t  = elapsed_cycles() ;
  Stream2 = ZfpCompress(Array, nx2, ny2, nz2, maxerror, precision , &Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  nsize = ZfpArrayDims(Stream2, d, s, Ssize) ;
  memset(brray, 0, sizeof(float)*nx*ny*nz) ;
  t  = elapsed_cycles() ;
  ZfpExpand(brray, nx2, ny2, nz2, Stream2, Ssize2) ;
  t = elapsed_cycles() -t ;
  t_ns = cycles_to_ns(t) ;
  fprintf(stderr,"cycles = %ld, time = %8.0f (%f ns/pt), %f MB/s\n", t, t_ns, t_ns/(nx*ny*nz), 4000.0/(t_ns/(nx*ny*nz))) ;
  AnalyzeCompressionErrors(Array, brray, nx*ny*nz, 0.0f, "Ctest2");     // evaluate errors
}