//
// Copyright (C) 2022  Environnement Canada
//
// This is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// Author:
//     M. Valin,   Recherche en Prevision Numerique, 2022
//
#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>

#include <misc_timers.h>
#include <rmn/bi_endian_pack.h>

#if !defined(NPTS)
#define NPTS 524289
#endif

#define NTIMES 100

// syntax and functional test for the bi endian pack/unpack macros and functions
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
  uint32_t peek_u, peek_u2 ;
  int32_t peek_i, peek_i2 ;

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

// ========================== functional and syntax test ==========================
  nbits = 12 ;
  // insert in Little Endian mode (unsigned)
  LeStreamInit(&ple, packedle, sizeof(packedle), 0) ;  // force unknown stream mode
  printf("Stream mode = %s (Unknown expected)", StreamMode(ple)) ;
  LE64_STREAM_INSERT_BEGIN(ple) ;                      // this should force insert (write) mode
  printf(", then %s (Write expected)\n", StreamMode(ple)) ;
  LeStreamInsert(&ple, unpacked, nbits, -NPTS) ;
  printf("packedle %2d  : ", nbits) ; for(i=7 ; i>=0 ; i--) printf("%8.8x ", packedle[i]); printf("\n") ;
  // rewind and extract in Little Endian mode (unsigned)
  if(0) LeStreamInit(&ple, packedle, sizeof(packedle), BIT_XTRACT_MODE) ;  // not necessary, syntax check
  LeStreamRewind(&ple) ;                               // this should force extract (read) mode
  printf("Stream mode = %s (Read expected)\n", StreamMode(ple)) ;
  LE64_PEEK_NBITS(ple.accum, ple.xtract, peek_u, nbits) ;  // explicit peek macro (unsigned)
  LE64_STREAM_PEEK_NBITS(ple, peek_u2, nbits) ;            // stream peek macro (unsigned)
  printf("peekle(u) = %8.8x %8.8x\n", peek_u, peek_u2) ;
  for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
  LeStreamXtract(&ple, restored, nbits, NPTS) ;
  printf("restoredle %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;

  // insert in Big Endian mode (unsigned)
  BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_INSERT_MODE) ;
  printf("Stream mode = %s (Write expected)\n", StreamMode(pbe)) ;
  BE64_STREAM_INSERT_BEGIN(pbe) ;
  BeStreamInsert(&pbe, unpacked, nbits, -NPTS) ;
  printf("packedbe %2d  : ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", packedbe[i]); printf("\n") ;
  // rewind and extract in Big Endian mode (unsigned)
  if(0) BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_XTRACT_MODE) ;  // not necessary, syntax check
  BeStreamRewind(&pbe) ;                               // this should force extract (read) mode
  printf("Stream mode = %s (Read expected)\n", StreamMode(pbe)) ;
  BE64_PEEK_NBITS(pbe.accum, pbe.xtract, peek_u, nbits) ;  // explicit peek macro (unsigned)
  BE64_STREAM_PEEK_NBITS(pbe, peek_u2, nbits) ;            // stream peek macro (unsigned)
  printf("peekbe(u) = %8.8x %8.8x\n", peek_u, peek_u2) ;
  for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
  BeStreamXtract(&pbe, restored, nbits, NPTS) ;
  printf("restoredbe %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", restored[i]); printf("\n") ;
  printf("\n") ;

  // insert in Little Endian mode (signed)
  LeStreamInit(&ple, packedle, sizeof(packedle), BIT_INSERT_MODE) ;
  LeStreamInsert(&ple, (void *) unpacked_signed, nbits, -NPTS) ;
  printf("packedle %2d  : ", nbits) ; for(i=7 ; i>=0 ; i--) printf("%8.8x ", packedle[i]); printf("\n") ;
  // rewind and extract in Little Endian mode (signed)
  if(0) LeStreamInit(&ple, packedle, sizeof(packedle), BIT_XTRACT_MODE) ;  // not necessary, syntax check
  LeStreamRewind(&ple) ;                               // this should force extract (read) mode
  LE64_PEEK_NBITS((int64_t) ple.accum, ple.xtract, peek_i, nbits) ;  // explicit peek macro (signed)
  LE64_STREAM_PEEK_NBITS_S(ple, peek_i2, nbits) ;                    // stream peek macro (signed)
  printf("peekle(s) = %8.8x %8.8x\n", peek_i, peek_i2) ;
  for(i=0 ; i<NPTS ; i++) signed_restored[i] = 0xFFFFFFFFu ;
  LeStreamXtractSigned(&ple, signed_restored, nbits, NPTS) ;
  printf("restoredle %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", signed_restored[i]); printf("\n") ;

  // insert in Big Endian mode (signed)
  BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_INSERT_MODE) ;
  BeStreamInsert(&pbe, (void *) unpacked_signed, nbits, -NPTS) ;
  printf("packedbe %2d  : ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", packedbe[i]); printf("\n") ;
  // rewind and extract in Big Endian mode (signed)
  if(0) BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_XTRACT_MODE) ;  // not necessary, syntax check
  BeStreamRewind(&pbe) ;                               // this should force extract (read) mode
  BE64_PEEK_NBITS((int64_t) pbe.accum, pbe.xtract, peek_i, nbits) ;  // explicit peek macro (signed)
  BE64_STREAM_PEEK_NBITS_S(pbe, peek_i2, nbits) ;                    // stream peek macro (signed)
  printf("peekbe(s) = %8.8x %8.8x\n", peek_i, peek_i2) ;
  for(i=0 ; i<NPTS ; i++) signed_restored[i] = 0xFFFFFFFFu ;
  BeStreamXtractSigned(&pbe, signed_restored, nbits, NPTS) ;
  printf("restoredbe %2d: ", nbits) ; for(i=0 ; i<8 ; i++) printf("%8.8x ", signed_restored[i]); printf("\n") ;
  printf("\n") ;

// ==================================================== timing tests ====================================================
  printf("%6d items,               insert                            extract (unsigned)                       extract (signed)\n", NPTS) ;
  for(nbits = 1 ; nbits <= 32 ; nbits += 1){
    mask = RMask(nbits) ;
    for(i=0 ; i<NPTS ; i++)    unpacked[i] = (i + 15) ;
    for(i=0 ; i<NPTS   ; i+=2) unpacked_signed[i] = -(((unpacked[i]) & mask) >> 1) ;
    for(i=1 ; i<NPTS-1 ; i+=2) unpacked_signed[i] =  (((unpacked[i]) & mask) >> 1) ;
    printf("nbits = %2d", nbits) ;

//  time little endian insertion
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, \
    LeStreamInit(&ple, packedle, sizeof(packedle), BIT_INSERT_MODE) ; LeStreamInsert(&ple, unpacked, nbits, -NPTS) ) ;
    printf(", %6.2f ns/pt (le)", tavg*nano/NPTS);

//  time big endian insertion
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, \
    BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_INSERT_MODE) ; BeStreamInsert(&pbe, unpacked, nbits, -NPTS) ) ;
    printf(", %6.2f ns/pt (be)", tavg*nano/NPTS);

//  time little endian unsigned extraction + error check
    for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
    LeStreamInit(&ple, packedle, sizeof(packedle), BIT_XTRACT_MODE) ;
    LeStreamRewind(&ple) ;
    LeStreamXtract(&ple, restored, nbits, NPTS) ;
    mask = RMask(nbits) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if((restored[i] & mask) != (unpacked[i] & mask) ) errors++ ;
    }
    errorsle = errors ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, \
    LeStreamInit(&ple, packedle, sizeof(packedle), BIT_XTRACT_MODE) ; LeStreamXtract(&ple, restored, nbits, NPTS) ) ;
    printf(", = %6.2f ns/pt (le)", tavg*nano/NPTS);

//  time big endian unsigned extraction + error check
    for(i=0 ; i<NPTS ; i++) restored[i] = 0xFFFFFFFFu ;
    BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_XTRACT_MODE) ;
    BeStreamRewind(&pbe) ;
    BeStreamXtract(&pbe, restored, nbits, NPTS) ;
    mask = RMask(nbits) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if((restored[i] & mask) != (unpacked[i] & mask) ) errors++ ;
    }
    errorsbe = errors ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, \
    BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_XTRACT_MODE) ; BeStreamXtract(&pbe, restored, nbits, NPTS) ) ;
    printf(", %6.2f ns/pt (be)", tavg*nano/NPTS);

//  time little endian signed extraction + error check
    for(i=0 ; i<NPTS ; i++) signed_restored[i] = 0xFFFFFFFFu ;
    LeStreamInit(&ple, packedle, sizeof(packedle), BIT_INSERT_MODE) ;
    LeStreamInsert(&ple, (void *) unpacked_signed, nbits, -NPTS) ;
    LeStreamInit(&ple, packedle, sizeof(packedle), BIT_XTRACT_MODE) ;
    LeStreamXtractSigned(&ple, signed_restored, nbits, NPTS) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if(unpacked_signed[i] != signed_restored[i]) errors++;
    }
    errorsles = errors ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, \
    LeStreamInit(&ple, packedle, sizeof(packedle), BIT_XTRACT_MODE) ; LeStreamXtractSigned(&ple, signed_restored, nbits, NPTS) ) ;
    printf(", %6.2f ns/pt (les)", tavg*nano/NPTS);

//  time big endian signed extraction + error check
    for(i=0 ; i<NPTS ; i++) signed_restored[i] = 0xFFFFFFFFu ;
    BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_INSERT_MODE) ;
    BeStreamInsert(&pbe, (void *) unpacked_signed, nbits, -NPTS) ;
    BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_XTRACT_MODE) ;
    BeStreamXtractSigned(&pbe, signed_restored, nbits, NPTS) ;
    errors = 0 ;
    for(i=0 ; i<NPTS ; i++){
      if(unpacked_signed[i] != signed_restored[i]) errors++;
    }
    errorsbes = errors ;
    TIME_LOOP(tmin, tmax, tavg, NTIMES, NPTS, buf, bufsiz, \
    BeStreamInit(&pbe, packedbe, sizeof(packedbe), BIT_XTRACT_MODE) ; BeStreamXtractSigned(&pbe, signed_restored, nbits, NPTS) ) ;
    printf(", %6.2f ns/pt (bes)", tavg*nano/NPTS);
//
    printf(" (%d/%d/%d/%d errors)", errorsle, errorsbe, errorsles, errorsbes);
//
//     printf(" %8.8x %8.8x", signed_restored[NPTS/2 + nbits], signed_restored[NPTS/2 + nbits + 1]);
    printf("\n");
  }
}
