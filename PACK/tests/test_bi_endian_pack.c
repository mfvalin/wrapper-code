#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>

#include <misc_timers.h>
#include <bi_endian_pack.h>

#if !defined(NPTS)
#define NPTS 4097
#endif

#define NTIMES 1000

#define TIME_CODE(tmin, tmax, tavg, niter, CODE) { \
    tmin = 1000000000.0 ; tmax = 0.0 ; tavg = 0.0 ;uint64_t t ; \
    for(j=0 ; j < niter ; j++) { t = elapsed_cycles() ; \
      CODE ; \
      t = elapsed_cycles() -t ; \
      tavg += t ; tmin = (t < tmin) ? t : tmin ; tmax = (t > tmax) ? t : tmax ; \
    } tavg /= niter ; \
  }


int main(int argc, char **argv){
  uint32_t unpacked[NPTS], packedle[NPTS], packedbe[NPTS], restored[NPTS] ;
  bitstream ple, pbe ;
  int i, j, nbits, errors ;
  uint32_t mask ;
  uint64_t t0, t1, tmin, tmax, tavg, freq ;
  double nano ;
  uint64_t t[NTIMES] ;

  freq = cycles_counter_freq() ;
  nano = 1000000000 ;
  nano /= freq ;
  for(i=0 ; i<NTIMES ; i++) t[i] = 0 ;

  for(i=0 ; i<NPTS ; i++) unpacked[i] = i + 16 ;
  printf("original  : ") ;
  for(i=0 ; i<8 ; i++) printf("%8.8x ", unpacked[i]); printf("\n") ;

  nbits = 12 ;
  LeStreamInit(&ple, packedle) ;
  LeStreamInsert(&ple, unpacked, nbits, -NPTS) ;
  printf("packedle %2d  : ", nbits) ; for(i=7 ; i>=0 ; i--) printf("%8.8x ", packedle[i]); printf("\n") ;
  LeStreamInit(&ple, packedle) ;
  for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
  LeStreamXtract(&ple, restored, nbits, NPTS) ;
  printf("restoredle %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;

  BeStreamInit(&pbe, packedbe) ;
  BeStreamInsert(&pbe, unpacked, nbits, -NPTS) ;
  printf("packedbe %2d  : ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", packedbe[i]); printf("\n") ;
  BeStreamInit(&pbe, packedbe) ;
  for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
  BeStreamXtract(&pbe, restored, nbits, NPTS) ;
  printf("restoredbe %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;
  printf("\n") ;

  for(nbits = 1 ; nbits <= 16 ; nbits += 1){
    TIME_LOOP(tmin, tmax, tavg, NTIMES, LeStreamInit(&ple, packedle) ; LeStreamInsert(&ple, unpacked, nbits, -NPTS) ) ;
    printf("nbits = %2d, ns(le) = %6.0f (%6.0f), %6.2f ns/pt", nbits, tmin*nano, tavg*nano, tavg*nano/NPTS);
  //     printf("packedle  : ") ;
  //     for(i=7 ; i>=0 ; i--) printf("%8.8x ", packedle[i]); printf("\n") ;

    TIME_LOOP(tmin, tmax, tavg, NTIMES, LeStreamInit(&pbe, packedbe) ; BeStreamInsert(&pbe, unpacked, nbits, -NPTS) ) ;
    printf(", ns(be) = %6.0f (%6.0f), %6.2f ns/pt", tmin*nano, tavg*nano, tavg*nano/NPTS);
  //     printf("packedbe  : ") ;
  //     for(i=0 ; i<8 ; i++) printf("%8.8x ", packedbe[i]); printf("\n") ;

    for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
    LeStreamInit(&ple, packedle) ;
    LeStreamXtract(&ple, restored, nbits, NPTS) ;
//     printf("restoredle: ") ;
//     for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;
    mask = RMask(nbits) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if((restored[i] & mask) != (unpacked[i] & mask) ) {
        if(errors < 1) printf("i = %4d, expected = %8.8x, got = %8.8x, raw = %8.8x\n",
                              i, unpacked[i] & mask, restored[i] & mask, unpacked[i]) ;
        errors++ ;
      }
    }
    printf(", errors (le) = %d / %d", errors, NPTS) ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, LeStreamInit(&ple, packedle) ; LeStreamXtract(&ple, restored, nbits, NPTS) ) ;
    printf(", ns(le) = %6.1f (%6.0f), %6.2f ns/pt", tmin*nano, tavg*nano, tavg*nano/NPTS);

    BeStreamInit(&pbe, packedbe) ;
    for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
    BeStreamXtract(&pbe, restored, nbits, NPTS) ;
//     printf("restoredbe: ") ;
//     for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;
    mask = RMask(nbits) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if((restored[i] & mask) != (unpacked[i] & mask) ) {
        if(errors < 1) printf("i = %4d, expected = %8.8x, got = %8.8x, raw = %8.8x\n",
                              i, unpacked[i] & mask, restored[i] & mask, unpacked[i]) ;
        errors++ ;
      }
    }
    printf(", errors (be) = %d / %d", errors, NPTS) ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, BeStreamInit(&pbe, packedbe) ; BeStreamXtract(&pbe, restored, nbits, NPTS) ) ;
    printf(", ns(be) = %6.0f (%6.0f), %6.2f ns/pt\n", tmin*nano, tavg*nano, tavg*nano/NPTS);
  }
}
