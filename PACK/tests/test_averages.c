#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <rmn/misc_timers.h>
// #include <misc_swizzle.h>

#define NTIMES 10

// averages
void average_rows_8x2(float *a, float *c, int lrow) ;
void average_rows_16x2(float *a, float *c, int lrow) ;
void average_rows_32x2(float *a, float *c, int lrow) ;
void average_rows_64x2(float *a, float *c, int lrow) ;
void averages_64x8(void *a, void *a2, void *a4, void *a8, int lrow);
void averages_64x64(void *a, void *a2, void *a4, void *a8, int lrow);
void expand_rows_8x2(void *a, void *c, int lrow) ;
void expand_rows_16x2(void *a, void *c, int lrow) ;
void expand_rows_32x2(void *a, void *c, int lrow) ;
void expand_rows_64x2(void *a, void *c, int lrow) ;
void expand_rows_nx2(void *a, void *c, int n, int lrow) ;
void expand_rows_64x8(void *a, void *c, int lrow) ;
void expand_rows_64x64(void *a, void *c, int lrow) ;
// misc_swizzle
void SplitEvenOdd_8(void *a, void *e, void *o) ;
void ShuffleEvenOdd_8(void *a, void *e, void *o) ;
void SplitEvenOdd_16(void *a, void *e, void *o) ;
void ShuffleEvenOdd_16(void *a, void *e, void *o) ;
void SplitEvenOdd_32(void *a, void *e, void *o) ;
void ShuffleEvenOdd_32(void *a, void *e, void *o) ;
void SplitEvenOdd_64(void *a, void *e, void *o) ;
void ShuffleEvenOdd_64(void *a, void *e, void *o) ;
void SplitEvenOdd_n(void *pa, void *pe, void *po, int n) ;
void ShuffleEvenOdd_n(void *pa, void *pe, void *po, int n) ;

int main(int argc, char **argv){
  float a[64*64], a2[32*32], a4[16*16], a8[8*8], e[64*32], o[64*32] ;
  float b[64*64], b2[32*32], b4[16*16], b8[8*8] ;
  int i, j, errors ;
  uint64_t t0, t1, t[NTIMES*1000], freq, tmin ;
  double nano ;
  double tavg ;
  char *test ;

  freq = cycles_counter_freq() ;
  nano = 1000000000 ;
  nano /= freq ;
  test = "timing" ;

  if(argc > 1) test = argv[1] ;

  if(strcmp(test, "expand") == 0) {
    printf("===== EXPAND test =====\n");
    for(i=0 ; i<32*32 ; i++) a2[i] = i ;
    for(i=0 ; i<64*64 ; i++) a[i] = -1 ;

    expand_rows_8x2(a, a2, 8) ;
    for(i=0 ; i< 4 ; i++) printf("%5.1f", a2[i])   ; printf("\n") ;
    for(i=0 ; i< 8 ; i++) printf("%5.1f", a[i+ 8]) ; printf("\n") ;
    for(i=0 ; i< 8 ; i++) printf("%5.1f", a[i+ 0]) ; printf("\n\n") ;

    expand_rows_16x2(a, a2, 16) ;
    for(i=0 ; i< 8 ; i++) printf("%5.1f", a2[i])   ; printf("\n") ;
    for(i=0 ; i<16 ; i++) printf("%5.1f", a[i+16]) ; printf("\n") ;
    for(i=0 ; i<16 ; i++) printf("%5.1f", a[i+ 0]) ; printf("\n\n") ;

    expand_rows_32x2(a, a2, 32) ;
    for(i=0 ; i<16 ; i++) printf("%5.1f", a2[i])   ; printf("\n") ;
    for(i=0 ; i<32 ; i++) printf("%5.1f", a[i+32]) ; printf("\n") ;
    for(i=0 ; i<32 ; i++) printf("%5.1f", a[i+ 0]) ; printf("\n\n") ;

    expand_rows_64x2(a, a2, 64) ;
    for(i=0 ; i<32 ; i++) printf("%5.1f", a2[i])   ; printf("\n") ;
    for(i=0 ; i<64 ; i+=2) printf("%5.1f", a[i+64]) ; printf("\n  ") ;
    for(i=1 ; i<64 ; i+=2) printf("%5.1f", a[i+64]) ; printf("\n") ;
    for(i=0 ; i<64 ; i+=2) printf("%5.1f", a[i+ 0]) ; printf("\n  ") ;
    for(i=1 ; i<64 ; i+=2) printf("%5.1f", a[i+ 0]) ; printf("\n\n") ;

    expand_rows_nx2(a, a2, 70, 70) ;
    for(i=0 ; i<35 ; i++) printf("%5.1f", a2[i])   ; printf("\n") ;
    for(i=0 ; i<70 ; i+=2) printf("%5.1f", a[i+70]) ; printf("\n  ") ;
    for(i=1 ; i<70 ; i+=2) printf("%5.1f", a[i+70]) ; printf("\n") ;
    for(i=0 ; i<70 ; i+=2) printf("%5.1f", a[i+ 0]) ; printf("\n  ") ;
    for(i=1 ; i<70 ; i+=2) printf("%5.1f", a[i+ 0]) ; printf("\n\n") ;
  }

  for(i=0 ; i<64*64 ; i++) a[i] = i ;

  if(strcmp(test, "shuffle") == 0) {
    printf("===== SHUFFLE test =====\n");
    for(i=0 ; i<128 ; i++) { e[i] = o[i] = -1.0f ; }
    SplitEvenOdd_8(a, e, o) ;
    for(i=0 ; i< 8 ; i++) printf("%5.1f", a[i  ]) ; printf("\n") ;
    for(i=0 ; i< 5 ; i++) printf("%10.1f",e[i  ]) ; printf("\n") ;
    for(i=0 ; i< 5 ; i++) printf("%10.1f",o[i  ]) ; printf("\n") ;
    for(i=0 ; i<128 ; i++) b[i] = -1.0 ;
    ShuffleEvenOdd_8(b, e, o) ;
    for(i=0 ; i< 9 ; i++) printf("%5.1f", b[i  ]) ; printf("\n\n") ;

    for(i=0 ; i<128 ; i++) { e[i] = o[i] = -1.0f ; }
    SplitEvenOdd_16(a, e, o) ;
    for(i=0 ; i<16 ; i++) printf("%5.1f", a[i  ]) ; printf("\n") ;
    for(i=0 ; i< 9 ; i++) printf("%10.1f",e[i  ]) ; printf("\n") ;
    for(i=0 ; i< 9 ; i++) printf("%10.1f",o[i  ]) ; printf("\n") ;
    for(i=0 ; i<128 ; i++) b[i] = -1.0 ;
    ShuffleEvenOdd_16(b, e, o) ;
    for(i=0 ; i<17 ; i++) printf("%5.1f", b[i  ]) ; printf("\n\n") ;

    for(i=0 ; i<128 ; i++) { e[i] = o[i] = -1.0f ; }
    SplitEvenOdd_32(a, e, o) ;
    for(i=0 ; i<32 ; i++) printf("%5.1f", a[i  ]) ; printf("\n") ;
    for(i=0 ; i<17 ; i++) printf("%10.1f",e[i  ]) ; printf("\n") ;
    for(i=0 ; i<17 ; i++) printf("%10.1f",o[i  ]) ; printf("\n") ;
    for(i=0 ; i<128 ; i++) b[i] = -1.0 ;
    ShuffleEvenOdd_32(b, e, o) ;
    for(i=0 ; i<33 ; i++) printf("%5.1f", b[i  ]) ; printf("\n\n") ;

    for(i=0 ; i<128 ; i++) { e[i] = o[i] = -1.0f ; }
    SplitEvenOdd_64(a, e, o) ;
    for(i=0 ; i<64 ; i+=2) printf("%5.1f", a[i  ]) ; printf("\n  ") ;
    for(i=1 ; i<64 ; i+=2) printf("%5.1f", a[i  ]) ; printf("\n") ;
    for(i=0 ; i<33 ; i++) printf("%5.1f",e[i  ]) ; printf("\n ") ;
    for(i=0 ; i<33 ; i++) printf("%5.1f",o[i  ]) ; printf("\n") ;
    for(i=0 ; i<128 ; i++) b[i] = -1.0 ;
    ShuffleEvenOdd_64(b, e, o) ;
    for(i=0 ; i<65 ; i+=2) printf("%5.1f", b[i  ]) ; printf("\n  ") ;
    for(i=1 ; i<65 ; i+=2) printf("%5.1f", b[i  ]) ; printf("\n\n") ;

    for(i=0 ; i<128 ; i++) { e[i] = o[i] = -1.0f ; }
    SplitEvenOdd_n(a, e, o, 67) ;
    for(i=0 ; i<67 ; i+=2) printf("%5.1f", a[i  ]) ; printf("\n  ") ;
    for(i=1 ; i<67 ; i+=2) printf("%5.1f", a[i  ]) ; printf("\n") ;
    for(i=0 ; i<35 ; i++) printf("%5.1f",e[i  ]) ; printf("\n ") ;
    for(i=0 ; i<34 ; i++) printf("%5.1f",o[i  ]) ; printf("\n") ;
    for(i=0 ; i<128 ; i++) b[i] = -1.0 ;
    ShuffleEvenOdd_n(b, e, o, 67) ;
    for(i=0 ; i<68 ; i+=2) printf("%5.1f", b[i  ]) ; printf("\n  ") ;
    for(i=1 ; i<68 ; i+=2) printf("%5.1f", b[i  ]) ; printf("\n") ;

    tmin = 1000000000;
    for(i=0 ; i<1000 ; i++){  // prime the pump
      SplitEvenOdd_n(a, b, b+64*32, 64*64) ;
      SplitEvenOdd_n(b, a, a+64*32, 64*64) ;
    }
    tavg = 0.0 ;
    for(i=0 ; i<NTIMES*1000 ; i++){
      t0 = elapsed_cycles() ;
      SplitEvenOdd_n(a, b, b+64*32, 64*64) ;
      SplitEvenOdd_n(b, a, a+64*32, 64*64) ;
      t1 = elapsed_cycles() ;
      t[i] = t1 - t0 ;
      tavg += t[i];
      tmin = (t1 - t0 < tmin) ? t1 - t0 : tmin ;
    }
    tavg /= 8192 ;
    printf("===== SplitEvenOdd_64x64 =====\n");
    for(i=0 ; i<NTIMES ; i++){
      printf("%7ld ",t[i]) ;
    }
    printf("cycles (min = %ld) (avg = %6.0f)\n", tmin, tavg) ;
    for(i=0 ; i<NTIMES ; i++){
      printf("%7.2f ",t[i]*nano) ;
    }
    printf("ns     (min = %7.2f) (avg = %7.2f)\n", tmin*nano, tavg*nano) ;
    for(i=0 ; i<NTIMES ; i++){
      printf("%7.2f ", t[i]*nano/8192) ;
    }
    printf("ns/pt  (min = %7.2f) (avg = %7.2f)\n\n", tmin*nano/8192, tavg*nano/8192) ;
  }

  if(strcmp(test, "average") == 0) {
    for(i=0 ; i<64*64 ; i++) a[i] = i ;
    printf("===== AVERAGE test =====\n");
    average_rows_8x2(a, a2, 8) ;
    for(i=0 ; i<8 ; i++) printf("%5.1f", a[i+8]) ; printf("\n") ;
    for(i=0 ; i<8 ; i++) printf("%5.1f", a[i  ]) ; printf("\n") ;
    for(i=0 ; i<4 ; i++) printf("%10.2f", a2[i]) ; printf("\n\n") ;

    average_rows_16x2(a, a2, 16) ;
    for(i=0 ; i<16 ; i++) printf("%8.2f", a[i+8]) ; printf("\n") ;
    for(i=0 ; i<16 ; i++) printf("%8.2f", a[i  ]) ; printf("\n") ;
    for(i=0 ; i<8  ; i++) printf("%16.2f", a2[i]) ; printf("\n\n") ;

    average_rows_32x2(a, a2, 32) ;
    for(i=0 ; i<16  ; i++) printf("%6.2f", a2[i]) ; printf("\n\n") ;

    average_rows_64x2(a, a2, 64) ;
    for(i=0 ; i<32  ; i++) printf("%6.2f", a2[i]) ; printf("\n\n") ;

    averages_64x8(a, a2, a4, a8, 64) ;
    for(i=0 ; i<16  ; i++) printf("%8.1f", a4[i+16]) ; printf("\n") ;
    for(i=0 ; i<16  ; i++) printf("%8.1f", a4[i   ]) ; printf("\n\n") ;
    for(i=0 ; i< 8  ; i++) printf("%8.1f", a8[i   ]) ; printf("\n\n") ;

    averages_64x64(a, a2, a4, a8, 64) ;
  //   for(i=0 ; i<32  ; i++) printf("%6.1f", a2[i+32*3]) ; printf("\n") ;
  //   for(i=0 ; i<32  ; i++) printf("%6.1f", a2[i   ]) ; printf("\n\n") ;
    for(i=0 ; i<8  ; i++) printf("%8.1f", a8[i+56]) ; printf("\n") ;
    for(i=0 ; i<8  ; i++) printf("%8.1f", a8[i   ]) ; printf("\n\n") ;
  }

  if(strcmp(test, "timing") == 0) {
    printf("===== TIMING test =====\n");
    tmin = 1000000000;
    for(i=0 ; i<NTIMES*1000 ; i++){
      t0 = elapsed_cycles() ;
      averages_64x64(a, a2, a4, a8, 64) ;
      t1 = elapsed_cycles() ;
      t[i] = t1 - t0 ;
      tmin = (t1 - t0 < tmin) ? t1 - t0 : tmin ;
    }
    printf("===== averages_64x64 =====\n");
    for(i=0 ; i<NTIMES ; i++){
      printf("%7ld ",t[i]) ;
    }
    printf("cycles (min = %ld)\n", tmin) ;
    for(i=0 ; i<NTIMES ; i++){
      printf("%7.1f ",t[i]*nano) ;
    }
    printf("ns     (min = %7.2f)\n", tmin*nano) ;
    for(i=0 ; i<NTIMES ; i++){
      printf("%7.2f ",t[i]*nano/4096) ;
    }
    printf("ns/pt  (min = %7.2f)\n\n", tmin*nano/4096) ;

    tmin = 1000000000;
    for(i=0 ; i<NTIMES*1000 ; i++){
      t0 = elapsed_cycles() ;
      expand_rows_64x64(b, a2, 64) ;
      t1 = elapsed_cycles() ;
      t[i] = t1 - t0 ;
      tmin = (t1 - t0 < tmin) ? t1 - t0 : tmin ;
    }
    printf("===== expand_rows_64x64 =====\n");
    for(i=0 ; i<NTIMES ; i++){
      printf("%7ld ",t[i]) ;
    }
    printf("cycles (min = %ld)\n", tmin) ;
    for(i=0 ; i<NTIMES ; i++){
      printf("%7.1f ",t[i]*nano) ;
    }
    printf("ns     (min = %7.2f)\n", tmin*nano) ;
    for(i=0 ; i<NTIMES ; i++){
      printf("%7.2f ",t[i]*nano/4096) ;
    }
    printf("ns/pt  (min = %7.2f)\n\n", tmin*nano/4096) ;
    averages_64x64(b, b2, b4, b8, 64) ;
    for(j=56 ; j>=0 ; j-=8){
      for(i=0 ; i<8 ; i++) printf("%8.2f",a8[i+j]);
      printf(" | ") ;
      for(i=0 ; i<8 ; i++) printf("%8.2f",b8[i+j]);
      printf("\n") ;
    }
    printf("\n") ;
  }

}
