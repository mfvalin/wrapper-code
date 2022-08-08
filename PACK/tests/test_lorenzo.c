
#include <stdio.h>
#include <stdint.h>
#include <misc_timers.h>

void LorenzoPredict2D(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);
void LorenzoUnpredict2D(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);

#define NI 64
#define LNI 131
#define LNIO 129
#define NJ 64
#define NTIMES 100000

int main(int argc, char **argv){
  int src[LNI*NJ], p[LNIO*NJ], dst[LNI*NJ] ;
  int i, j, ij, errors ;
  uint64_t t0, t1, t[NTIMES], freq, tmin ;
  double nano ;

  freq = cycles_counter_freq() ;
  nano = 1000000000 ;
  nano /= freq ;
#if defined(__AVX2__) && defined(__x86_64__) && defined(WITH_SIMD)
  printf("using SIMD version\n");
#endif
  printf("NI = %d, NJ = %d, LNI = %d, LNIO = %d\n", NI, NJ, LNI, LNIO) ;

  for(j=0 ; j<NJ ; j++){
    ij = j*LNI ;
    for(i=0 ; i<LNI ; i++){
      src[ij+i] = i + j - (NI+NJ-1)/2 ;
      dst[ij+i] = -999 ;
    }
  }
  for(j=0 ; j<NJ ; j++){
    ij = j*LNIO ;
    for(i=0 ; i<LNIO ; i++){
      p[ij+i] = -999 ;
    }
  }
  if(NI < 13){
    for(j=NJ-1 ; j>=0 ; j--){
      ij = j*LNI ;
      for(i=0 ; i<NI ; i++){
        printf("%4d",src[ij+i]) ;
      }
      printf("\n");
    }
    printf("\n");
  }
  tmin = 1000000000;
  for(i=0 ; i<NTIMES ; i++){
    t0 = elapsed_cycles() ;
    LorenzoPredict2D(src, p, NI, LNI, LNIO, NJ) ;
    t1 = elapsed_cycles() ;
    t[i] = t1 - t0 ;
    tmin = (t1 - t0 < tmin) ? t1 - t0 : tmin ;
  }
  printf("tmin LorenzoPredict2D   = %ld cycles, %8.1f ns, %8.2f ns/pt\n", tmin, tmin*nano, tmin*nano/NI/NJ) ;
  if(NI < 13){
    for(j=NJ-1 ; j>=0 ; j--){
      ij = j*LNIO ;
      for(i=0 ; i<NI ; i++){
        printf("%4d",p[ij+i]) ;
      }
      printf("\n");
    }
    printf("\n");
  }

  tmin = 1000000000;
  for(i=0 ; i<NTIMES ; i++){
    t0 = elapsed_cycles() ;
    LorenzoUnpredict2D(dst, p, NI, LNI, LNIO, NJ) ;
    t1 = elapsed_cycles() ;
    t[i] = t1 - t0 ;
    tmin = (t1 - t0 < tmin) ? t1 - t0 : tmin ;
  }
  printf("tmin LorenzoUnpredict2D = %ld cycles, %8.1f ns, %8.2f ns/pt\n", tmin, tmin*nano, tmin*nano/NI/NJ) ;
  if(NI < 13){
    for(j=NJ-1 ; j>=0 ; j--){
      ij = j*LNI ;
      for(i=0 ; i<NI ; i++) printf("%4d",dst[ij+i]) ; 
//       printf(" | ");
//       for(i=0 ; i<NI ; i++) printf("%4d",src[ij+i]) ;
      printf("\n");
    }
  }
  errors = 0 ;
  for(j=0 ; j<NJ ; j++){
    ij = j*LNI ;
    for(i=0 ; i<NI ; i++){
      if(dst[ij+i] != src[ij+i]) errors++ ;
    }
  }
  printf("\nerrors = %d\n",errors);
}
