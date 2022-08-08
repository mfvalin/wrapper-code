#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>

#define NI 8193
#define LNI 8200
#define NJ 10001
#define ITER 10

int float_info_simple(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs);
int float_info_no_missing(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs);
int float_info_missing(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs, float *spval, uint32_t spmask);
int float_info(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs, float *spval, uint32_t spmask);

int main(int argc, char **argv){
  float z[LNI*NJ] ;
  float maxval[ITER], minval[ITER], minabs[ITER], spval ;
  int i, j, nspecial ;
  struct timeval t1, t2 ;
  uint64_t t0 ;
  double T0 ;
  int *ispval = (int *) &spval ;

  for(i=0 ; i < LNI*NJ ; i++) z[i] = (i & 0xFFFF) ;
  z[LNI/2] = -0.12345f ;
  z[3] = -123.456f ;
//   spval = -1.234567f ;
  spval = 0xFFFF ;
  fprintf(stderr,"spval = %f (%8.8x)\n", spval, *ispval) ;
  z[6] = spval ;
  z[5] = spval ;
  z[4] = spval ;
//   z[2] = spval ;
  nspecial = 0 ;

  gettimeofday(&t1, NULL) ;
  for(j=0 ; j < ITER ; j++) {
    nspecial = float_info(z, NI, LNI, NJ, &maxval[j], &minval[j], &minabs[j], &spval, 0x3FF) ;
  }
  gettimeofday(&t2, NULL) ;
  t0 = t2.tv_usec - t1.tv_usec ;
  t0 += (t2.tv_sec - t1.tv_sec) * 1000000 ;
  T0 = t0 ;
  T0 *= 1000. ;
  T0 /= ITER ;
  fprintf(stderr,"nspecial = %d, time = %ld, %f ns/pt %f GB/s\n", nspecial, t0, T0/(NI*NJ), NI*NJ/T0*4);
  for(j=1 ; j < ITER ; j++) {
    maxval[0] = (maxval[j] > maxval[0]) ? maxval[j] : maxval[0] ;
  }
  fprintf(stderr,"maxval = %f, minval = %f, amin = %f\n\n", maxval[0], minval[0], minabs[0]);

  gettimeofday(&t1, NULL) ;
  for(j=0 ; j < ITER ; j++) {
    nspecial = float_info(z, NI, LNI, NJ, &maxval[j], &minval[j], &minabs[j], NULL, 0) ;
  }
  gettimeofday(&t2, NULL) ;
  t0 = t2.tv_usec - t1.tv_usec ;
  t0 += (t2.tv_sec - t1.tv_sec) * 1000000 ;
  T0 = t0 ;
  T0 *= 1000. ;
  T0 /= ITER ;
  fprintf(stderr,"nspecial = %d, time = %ld, %f ns/pt %f GB/s\n", nspecial, t0, T0/(NI*NJ), NI*NJ/T0*4);
  for(j=1 ; j < ITER ; j++) {
    maxval[0] = (maxval[j] > maxval[0]) ? maxval[j] : maxval[0] ;
  }
  fprintf(stderr,"maxval = %f, minval = %f, amin = %f\n\n", maxval[0], minval[0], minabs[0]);

  gettimeofday(&t1, NULL) ;
  for(j=0 ; j < ITER ; j++) {
    nspecial = float_info_simple(z, NI, LNI, NJ, &maxval[j], &minval[j], &minabs[j]) ;
  }
  gettimeofday(&t2, NULL) ;
  t0 = t2.tv_usec - t1.tv_usec ;
  t0 += (t2.tv_sec - t1.tv_sec) * 1000000 ;
  T0 = t0 ;
  T0 *= 1000. ;
  T0 /= ITER ;
  fprintf(stderr,"nspecial = %d, time = %ld, %f ns/pt %f GB/s\n", nspecial, t0, T0/(NI*NJ), NI*NJ/T0*4);
  for(j=1 ; j < ITER ; j++) {
    maxval[0] = (maxval[j] > maxval[0]) ? maxval[j] : maxval[0] ;
  }
  fprintf(stderr,"maxval = %f, minval = %f, amin = %f\n\n", maxval[0], minval[0], minabs[0]);
}

