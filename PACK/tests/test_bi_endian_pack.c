#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>

#include <misc_timers.h>
#include <bi_endian_pack.h>

#if !defined(NPTS)
#define NPTS 4097
#endif

#define NTIMES 1000

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

  for(nbits = 1 ; nbits <= 32 ; nbits += 1){
    tmin = 1000000000 ;
    for(j=0 ; j < NTIMES ; j++) {
      t0 = elapsed_cycles() ;
      LeStreamInit(&ple, packedle) ;
      LeStreamInsert(&ple, unpacked, nbits, -NPTS) ;
      t1 = elapsed_cycles() ;
      t[j] = t1 - t0 ;
      tmin = (t[j] < tmin) ? t[j] : tmin ;
    }
    printf("nbits = %2d, ns(le) = %6.1f, %6.3f ns/pt] ", nbits, tmin*nano, tmin*nano/NPTS);
  //     printf("packedle  : ") ;
  //     for(i=7 ; i>=0 ; i--) printf("%8.8x ", packedle[i]); printf("\n") ;

    tmin = 1000000000 ;
    for(j=0 ; j < NTIMES ; j++) {
      t0 = elapsed_cycles() ;
      LeStreamInit(&pbe, packedbe) ;
      BeStreamInsert(&pbe, unpacked, nbits, -NPTS) ;
      t1 = elapsed_cycles() ;
      t[j] = t1 - t0 ;
      tmin = (t[j] < tmin) ? t[j] : tmin ;
    }
    printf(", ns(be) = %6.1f, %6.3f ns/pt]", tmin*nano, tmin*nano/NPTS);
  //     printf("packedbe  : ") ;
  //     for(i=0 ; i<8 ; i++) printf("%8.8x ", packedbe[i]); printf("\n") ;

    for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
    LeStreamInit(&ple, packedle) ;
    LeStreamXtract(&ple, restored, nbits, NPTS) ;
//     printf("restoredle: ") ;
//     for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;
    mask = MaskNbits(nbits) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if((restored[i] & mask) != (unpacked[i] & mask) ) {
        if(errors < 1) printf("i = %4d, expected = %8.8x, got = %8.8x, raw = %8.8x\n",
                              i, unpacked[i] & mask, restored[i] & mask, unpacked[i]) ;
        errors++ ;
      }
    }
    printf(", mask = %8.8x, errors (le) = %d / %d", mask, errors, NPTS) ;
    tmin = 1000000000 ;
    for(j=0 ; j < NTIMES ; j++) {
      t0 = elapsed_cycles() ;
      LeStreamInit(&ple, packedle) ;
      LeStreamXtract(&ple, restored, nbits, NPTS) ;
      t1 = elapsed_cycles() ;
      t[j] = t1 - t0 ;
      tmin = (t[j] < tmin) ? t[j] : tmin ;
    }
    printf(", ns(le) = %6.1f, %6.3f ns/pt]", tmin*nano, tmin*nano/NPTS);

    BeStreamInit(&pbe, packedbe) ;
    for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
    BeStreamXtract(&pbe, restored, nbits, NPTS) ;
//     printf("restoredbe: ") ;
//     for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;
    mask = MaskNbits(nbits) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if((restored[i] & mask) != (unpacked[i] & mask) ) {
        if(errors < 1) printf("i = %4d, expected = %8.8x, got = %8.8x, raw = %8.8x\n",
                              i, unpacked[i] & mask, restored[i] & mask, unpacked[i]) ;
        errors++ ;
      }
    }
    printf(", errors (be) = %d / %d", errors, NPTS) ;
    tmin = 1000000000 ;
    for(j=0 ; j < NTIMES ; j++) {
      t0 = elapsed_cycles() ;
      BeStreamInit(&pbe, packedbe) ;
      BeStreamXtract(&pbe, restored, nbits, NPTS) ;
      t1 = elapsed_cycles() ;
      t[j] = t1 - t0 ;
      tmin = (t[j] < tmin) ? t[j] : tmin ;
    }
    printf(", ns(be) = %6.1f, %6.3f ns/pt]\n", tmin*nano, tmin*nano/NPTS);
  }
}
