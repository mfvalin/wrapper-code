#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>

#include <misc_timers.h>
#include <bi_endian_pack.h>

#if !defined(NPTS)
#define NPTS 524289
#endif

#define NTIMES 100

int main(int argc, char **argv){
  uint32_t unpacked[NPTS], packedle[NPTS], packedbe[NPTS], restored[NPTS] ;
  int32_t unpacked_signed[NPTS], signed_restored[NPTS] ;
  bitstream ple, pbe ;
  int i, j, nbits, errors, errorsle, errorsbe, errorsles, errorsbes ;
  uint32_t mask ;
  uint64_t t0, t1, tmin, tmax, tavg, freq ;
  double nano ;
  uint64_t t[NTIMES] ;
  char buf[1024] ;
  size_t bufsiz = sizeof(buf) ;

  freq = cycles_counter_freq() ;
  nano = 1000000000 ;
  nano /= freq ;
  for(i=0 ; i<NTIMES ; i++) t[i] = 0 ;

  for(i=0 ; i<NPTS ; i++) unpacked[i] = i + 16 ;
  for(i=0 ; i<NPTS   ; i+=2) unpacked_signed[i] = -unpacked[i] ;
  for(i=1 ; i<NPTS-1 ; i+=2) unpacked_signed[i] =  unpacked[i] ;
  printf("original(u)  : ") ;
  for(i=0 ; i<8 ; i++) printf("%8.8x ", unpacked[i]); printf("\n") ;
  printf("original(s)  : ") ;
  for(i=0 ; i<8 ; i++) printf("%8.8x ", unpacked_signed[i]); printf("\n") ;
  printf("\n") ;

  nbits = 12 ;
  LeStreamInit(&ple, packedle, sizeof(packedle)) ;
  LeStreamInsert(&ple, unpacked, nbits, -NPTS) ;
  printf("packedle %2d  : ", nbits) ; for(i=7 ; i>=0 ; i--) printf("%8.8x ", packedle[i]); printf("\n") ;
  LeStreamInit(&ple, packedle, sizeof(packedle)) ;
  for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
  LeStreamXtract(&ple, restored, nbits, NPTS) ;
  printf("restoredle %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;

  BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ;
  BeStreamInsert(&pbe, unpacked, nbits, -NPTS) ;
  printf("packedbe %2d  : ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", packedbe[i]); printf("\n") ;
  BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ;
  for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
  BeStreamXtract(&pbe, restored, nbits, NPTS) ;
  printf("restoredbe %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;
  printf("\n") ;

  LeStreamInit(&ple, packedle, sizeof(packedle)) ;
  LeStreamInsert(&ple, (void *) unpacked_signed, nbits, -NPTS) ;
  printf("packedle %2d  : ", nbits) ; for(i=7 ; i>=0 ; i--) printf("%8.8x ", packedle[i]); printf("\n") ;
  LeStreamInit(&ple, packedle, sizeof(packedle)) ;
  for(i=0 ; i<NPTS ; i++) signed_restored[i] = 0xFFFFFFFFu ;
  LeStreamXtractSigned(&ple, signed_restored, nbits, NPTS) ;
  printf("restoredle %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", signed_restored[i]); printf("\n") ;

  BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ;
  BeStreamInsert(&pbe, (void *) unpacked_signed, nbits, -NPTS) ;
  printf("packedbe %2d  : ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", packedbe[i]); printf("\n") ;
  BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ;
  for(i=0 ; i<NPTS ; i++) signed_restored[i] = 0xFFFFFFFFu ;
  BeStreamXtractSigned(&pbe, signed_restored, nbits, NPTS) ;
  printf("restoredbe %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", signed_restored[i]); printf("\n") ;
  printf("\n") ;

  printf("%6d points,              insert                            extract (unsigned)                       extract (signed)\n", NPTS) ;
  for(nbits = 1 ; nbits <= 32 ; nbits += 1){
    mask = RMask(nbits) ;
    for(i=0 ; i<NPTS ; i++)    unpacked[i] = (i + 15) ;
    for(i=0 ; i<NPTS   ; i+=2) unpacked_signed[i] = -(((unpacked[i]) & mask) >> 1) ;
    for(i=1 ; i<NPTS-1 ; i+=2) unpacked_signed[i] =  (((unpacked[i]) & mask) >> 1) ;
    printf("nbits = %2d", nbits) ;

//  time little endian insertion
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, LeStreamInit(&ple, packedle, sizeof(packedle)) ; LeStreamInsert(&ple, unpacked, nbits, -NPTS) ) ;
    printf(", %6.2f ns/pt (le)", tavg*nano/NPTS);

//  time big endian insertion
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, LeStreamInit(&pbe, packedbe, sizeof(packedbe)) ; BeStreamInsert(&pbe, unpacked, nbits, -NPTS) ) ;
    printf(", %6.2f ns/pt (be)", tavg*nano/NPTS);

//  time little endian unsigned extraction
    for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
    LeStreamInit(&ple, packedle, sizeof(packedle)) ;
    LeStreamXtract(&ple, restored, nbits, NPTS) ;
    mask = RMask(nbits) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if((restored[i] & mask) != (unpacked[i] & mask) ) errors++ ;
    }
    errorsle = errors ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, LeStreamInit(&ple, packedle, sizeof(packedle)) ; LeStreamXtract(&ple, restored, nbits, NPTS) ) ;
    printf(", = %6.2f ns/pt (le)", tavg*nano/NPTS);

//  time big endian unsigned extraction
    for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
    BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ;
    BeStreamXtract(&pbe, restored, nbits, NPTS) ;
    mask = RMask(nbits) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if((restored[i] & mask) != (unpacked[i] & mask) ) errors++ ;
    }
    errorsbe = errors ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ; BeStreamXtract(&pbe, restored, nbits, NPTS) ) ;
    printf(", %6.2f ns/pt (be)", tavg*nano/NPTS);

//  time little endian signed extraction
    for(i=0 ; i<NPTS ; i++) signed_restored[i] = 0xFFFFFFFFu ;
    LeStreamInit(&ple, packedle, sizeof(packedle)) ;
    LeStreamInsert(&ple, (void *) unpacked_signed, nbits, -NPTS) ;
    LeStreamInit(&ple, packedle, sizeof(packedle)) ;
    LeStreamXtractSigned(&ple, signed_restored, nbits, NPTS) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if(unpacked_signed[i] != signed_restored[i]) errors++;
    }
    errorsles = errors ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, LeStreamInit(&ple, packedle, sizeof(packedle)) ; LeStreamXtractSigned(&ple, signed_restored, nbits, NPTS) ) ;
    printf(", %6.2f ns/pt (les)", tavg*nano/NPTS);

//  time big endian signed extraction
    for(i=0 ; i<NPTS ; i++) signed_restored[i] = 0xFFFFFFFFu ;
    BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ;
    BeStreamInsert(&pbe, (void *) unpacked_signed, nbits, -NPTS) ;
    BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ;
    BeStreamXtractSigned(&pbe, signed_restored, nbits, NPTS) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if(unpacked_signed[i] != signed_restored[i]) errors++;
    }
    errorsbes = errors ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, BeStreamInit(&pbe, packedbe, sizeof(packedbe)) ; BeStreamXtractSigned(&pbe, signed_restored, nbits, NPTS) ) ;
    printf(", %6.2f ns/pt (bes)", tavg*nano/NPTS);
//
    printf(" (%d/%d/%d/%d errors)", errorsle, errorsbe, errorsles, errorsbes);
//
//     printf(" %8.8x %8.8x", signed_restored[NPTS/2 + nbits], signed_restored[NPTS/2 + nbits + 1]);
    printf("\n");
  }
}
