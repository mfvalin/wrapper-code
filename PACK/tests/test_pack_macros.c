#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>

#include <pack_macros.h>
// need RMASK operators
#include <rmn/misc_operators.h>

#define TIMING_PACK 67108864
#define MAX_PACK 4096
#define MAX_ELEMENTS 23

uint32_t my_packer(uint32_t *packed, uint32_t *data, uint32_t nbits, int ndata){
  uint64_t acc ;
  uint32_t nfree ;
  uint32_t *stream = packed ;
  int i ;
  PACK32_PUT_INIT(acc, nfree) ;
  for(i = 0 ; i < ndata ; i++){
    PACK32_PUT(stream, acc, nfree, data[i], nbits) ;
  }
  PACK32_PUT_FLUSH(stream, acc, nfree) ;
  return stream - packed ;
}

uint32_t my_unpacker(uint32_t *packed, uint32_t *data, uint32_t nbits, int ndata){
  uint64_t acc ;
  uint32_t nfree ;
  uint32_t *stream = packed ;
  int i ;
  PACK32_GET_INIT(stream, acc, nfree) ;
  for(i = 0 ; i < ndata ; i++){
    PACK32_GET(stream, acc, nfree, data[i], nbits) ;
  }
  return stream - packed ;
}

int main(int argc, char **argv){
  uint64_t acc, t64, mask64, lmask64 ;
  uint32_t token ;
  uint32_t nfree, navail, nbits, mask, lmask ;
  uint32_t *stream ;
  uint32_t buffer[MAX_PACK] ;
  int i, nstream, errors ;
  struct timeval t1, t2 ;
  uint32_t buf[TIMING_PACK] ;
  uint32_t tok[TIMING_PACK] ;
  uint32_t lenpak ;

  nbits = 28 ;
  mask   = RMASK32(nbits) ;
  lmask  = LMASK32(nbits) ;
  nbits = 36 ;
  mask64  = RMASK64(nbits) ;
  lmask64 = LMASK64(nbits) ;
  fprintf(stderr,"mask  = %8.8x, mask64  = %16.16lx\n", mask,  mask64) ;
  fprintf(stderr,"lmask = %8.8x, lmask64 = %16.16lx\n", lmask, lmask64) ;

  nbits = 28 ;
  mask   = RMASK32(nbits) ;

  // pack
  stream = &buffer[0] ;
  PACK32_PUT_INIT(acc, nfree) ;
  for(i = 0 ; i < MAX_ELEMENTS ; i++){
    token = i & mask ;
    PACK32_PUT(stream, acc, nfree, token, nbits) ;
  }
  fprintf(stderr,"last in stream = %8.8x, acc = %16.16lx, free = %u\n", stream[-1], acc, nfree) ;
  PACK32_PUT_FLUSH(stream, acc, nfree) ;
  nstream = stream - buffer ;
  fprintf(stderr,"acc = %lu, free = %u, stream elements = %d\n", acc, nfree, nstream) ;
//   for(i = 0 ; i < nstream ; i++) fprintf(stderr," %8.8x", buffer[i]) ;
//   fprintf(stderr,"\n\n") ;

  for(i = 0 ; i < MAX_PACK ; i++) buffer[i] = 0 ;
  stream = &buffer[0] ;
  PACK32_PUT_INIT(acc, nfree) ;
  for(i = 0 ; i < MAX_ELEMENTS ; i++){
    t64 = i & mask ;
    PACK32_PUT_L(stream, acc, nfree, t64, nbits) ;
  }
  fprintf(stderr,"last in stream = %8.8x, acc = %16.16lx, free = %u\n", stream[-1], acc, nfree) ;
  PACK32_PUT_FLUSH(stream, acc, nfree) ;
  nstream = stream - buffer ;
  fprintf(stderr,"acc = %lu, free = %u, stream elements = %d\n", acc, nfree, nstream) ;
//   for(i = 0 ; i < nstream ; i++) fprintf(stderr," %8.8x", buffer[i]) ;
//   fprintf(stderr,"\n\n") ;

  // unpack
  stream = &buffer[0] ;
  PACK32_GET_INIT(stream, acc, navail) ;
  for(i = 0 ; i < MAX_ELEMENTS ; i++){
    PACK32_GET(stream, acc, navail, token, nbits) ;
//     fprintf(stderr," %8.8x", token) ;
  }
//   fprintf(stderr,"\n") ;

  // unpack
  stream = &buffer[0] ;
  PACK32_GET_INIT(stream, acc, navail) ;
  for(i = 0 ; i < MAX_ELEMENTS ; i++){
    PACK32_GET_L(stream, acc, navail, t64, nbits) ;
//     fprintf(stderr," %8.8lx", t64) ;
  }
//   fprintf(stderr,"\n") ;

  for(nbits = 2 ; nbits <= 32 ; nbits+=2) {
    mask = RMASK32(nbits) ;
    for(i = 0 ; i < TIMING_PACK ; i++) tok[i] = i & mask ;
    gettimeofday(&t1, NULL) ;
    lenpak = my_packer(buf, tok, nbits, TIMING_PACK) ;
    gettimeofday(&t2, NULL) ;
    t64  = (t2.tv_usec - t1.tv_usec) ;
    t64 += (t2.tv_sec - t1.tv_sec) * 1000000 ;
    fprintf(stderr,"PACK:  nbits = %2d, lenpak = %9d, %6ld us (%7.2g ns/pt), %8.3g PTS/s ", 
            nbits, lenpak, t64, t64*1000./TIMING_PACK, TIMING_PACK*1.0E6/t64) ;
    gettimeofday(&t1, NULL) ;
    lenpak = my_unpacker(buf, tok, nbits, TIMING_PACK) ;
    gettimeofday(&t2, NULL) ;
    t64  = (t2.tv_usec - t1.tv_usec) ;
    t64 += (t2.tv_sec - t1.tv_sec) * 1000000 ;
    errors = 0 ;
    for(i = 0 ; i < TIMING_PACK ; i++) if(tok[i] != (i & mask)) errors++ ;
//     fprintf(stderr,"UNPAK: nbits = %2d, lenpak = %9d, time = %6ld us (%7.2g ns/pt), PTS/s = %8.3g, errors = %d\n", 
//             nbits, lenpak, t64, t64*1000./TIMING_PACK, TIMING_PACK*1.0E6/t64, errors) ;
    fprintf(stderr,"| UNPAK: %6ld us (%7.2g ns/pt), %8.3g PTS/s, errors = %d\n", 
            t64, t64*1000./TIMING_PACK, TIMING_PACK*1.0E6/t64, errors) ;
  }
  return 0 ;
}
