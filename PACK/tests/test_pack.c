#include <stdio.h>
#include <misc_helpers.h>

#if 0
void unpack_stream(void *p, uint32_t *u, int nbits, int n){
  stream32 ps32 ;
  int i ;
  stream32_init(&ps32, p) ;
  if(nbits > 16){   // items are processed 1 by 1
    for(i=0 ; i<n ; i++) u[i] = stream32_get(&ps32, nbits) ;
  }else{            // items are processed 2 by 2
    for(i=0 ; i<n-1 ; i+=2) {
      u[i] = stream32_get(&ps32, nbits) ;
      u[i+1] = stream32_get_fast(&ps32, nbits) ;  // second item is safe if 2 items <= 32 bibts
    }
    for( ; i<n ; i+=1) { u[i] = stream32_get(&ps32, nbits) ; }
  }
}

void unpack_stream1(void *p, uint32_t *u, int nbits, int n){
  stream32 ps32 ;
  int i ;
  stream32_init(&ps32, p) ;
  for(i=0 ; i<n ; i+=1) {
    stream32_get_check(&ps32) ;
    u[i]   = stream32_get_fast(&ps32, nbits) ;
  }
}

void unpack_stream2(void *p, uint32_t *u, int nbits, int n){
  stream32 ps32 ;
  int i ;
  stream32_init(&ps32, p) ;
  for(i=0 ; i<n-1 ; i+=2) {
    stream32_get_check(&ps32) ;
    u[i]   = stream32_get_fast(&ps32, nbits) ;
    u[i+1] = stream32_get_fast(&ps32, nbits) ;
  }
  stream32_get_check(&ps32) ;
  for( ; i<n ; i+=1) { u[i]   = stream32_get_fast(&ps32, nbits) ; }
}

void unpack_stream4(void *p, uint32_t *u, int nbits, int n){
  stream32 ps32 ;
  int i ;
  stream32_init(&ps32, p) ;
  for(i=0 ; i<n-3 ; i+=4) {
    stream32_get_check(&ps32) ;
    u[i]   = stream32_get_fast(&ps32, nbits) ;
    u[i+1] = stream32_get_fast(&ps32, nbits) ;
    u[i+2] = stream32_get_fast(&ps32, nbits) ;
    u[i+3] = stream32_get_fast(&ps32, nbits) ;
  }
  stream32_get_check(&ps32) ;
  for( ; i<n ; i+=1) { u[i]   = stream32_get_fast(&ps32, nbits) ; }
}

void pack_stream(void *p, uint32_t *u, int nbits, int n){
  stream32 ps32 ;
  int i ;
  stream32_init(&ps32, p) ;
  if(nbits > 16){                // items are processed 1 by 1
    for(i=0 ; i<n ; i++) {
      stream32_put(&ps32, u[i], nbits) ;
    }
  }else{                         // items are processed 2 by 2
    for(i=0 ; i<n-1 ; i+=2) {
      stream32_put(&ps32, u[i], nbits) ;
      stream32_put_fast(&ps32, u[i+1], nbits) ;  // second item is safe if 2 items <= 32 bibts
    }
    for( ; i<n ; i+=1) { stream32_put(&ps32, u[i], nbits) ; }
  }
  stream32_flush(&ps32) ;
}

void pack_stream1(void *p, uint32_t *u, int nbits, int n){
  stream32 ps32 ;
  int i ;
  stream32_init(&ps32, p) ;
  if(nbits > 16){
    for(i=0 ; i<n ; i+=1){
      stream32_put_check(&ps32) ;
      stream32_put_fast(&ps32, u[i], nbits) ;
    }
  }else{
    for(i=0 ; i<n-1 ; i+=2){
      stream32_put_check(&ps32) ;
      stream32_put_fast(&ps32, u[i], nbits) ;
      stream32_put_fast(&ps32, u[i+1], nbits) ;
    }
    stream32_put_check(&ps32) ;
    for( ; i<n ; i+=1) { stream32_put_fast(&ps32, u[i], nbits) ; }
  }
  stream32_flush(&ps32) ;
}

void pack_stream2(void *p, uint32_t *u, int nbits, int n){
  stream32 ps32 ;
  int i ;
  stream32_init(&ps32, p) ;
  for(i=0 ; i<n-1 ; i+=2){
    stream32_put_check(&ps32) ;
    stream32_put_fast(&ps32, u[i], nbits) ;
    stream32_put_fast(&ps32, u[i+1], nbits) ;
  }
  stream32_put_check(&ps32) ;
  for( ; i<n ; i+=1) { stream32_put_fast(&ps32, u[i], nbits) ; }
  stream32_flush(&ps32) ;
}

void pack_stream4(void *p, uint32_t *u, int nbits, int n){
  stream32 ps32 ;
  int i ;
  stream32_init(&ps32, p) ;
  for(i=0 ; i<n-3 ; i+=4){
    stream32_put_check(&ps32) ;
    stream32_put_fast(&ps32, u[i], nbits) ;
    stream32_put_fast(&ps32, u[i+1], nbits) ;
    stream32_put_fast(&ps32, u[i+2], nbits) ;
    stream32_put_fast(&ps32, u[i+3], nbits) ;
  }
  stream32_put_check(&ps32) ;
  for( ; i<n ; i+=1) { stream32_put_fast(&ps32, u[i], nbits) ; }
  stream32_flush(&ps32) ;
}
#endif

#define NITEMS 655300

// test for the stream packers, from 1 to 32 bits
int main(int argc, char **argv){
  uint32_t packed[NITEMS] ;
  uint32_t src[NITEMS] ;
  uint32_t dst[NITEMS] ;
  uint32_t i, nbits, errors ;
  uint32_t mask ;
  uint64_t t0, t1, t2 ;

  for(nbits = 1 ; nbits < 33 ; nbits++){
    mask = 0xFFFFFFFF ;
    mask >>= (32-nbits) ;
    for(i=0 ; i<NITEMS ; i++) src[i] = i & mask ;
    for(i=0 ; i<NITEMS ; i++) packed[i] = -1 ;

    t0 = elapsed_us() ;
    pack_stream(packed, src, nbits, NITEMS) ;
    t1 = elapsed_us() ;
    unpack_stream(packed, dst, nbits, NITEMS) ;
    t2 = elapsed_us() ;

    errors = 0 ;
    for(i=0 ; i<NITEMS ; i++) {
      if(src[i] != dst[i]) {
        errors++ ;
        if(errors<5)printf("error at %d, src = %8.8x, dst = %8.8x, mask = %8.8x\n", i, src[i], dst[i], mask);
//         break ;
      }
    }
    printf("nbits = %2d, npts = %6d, errors = %6d, t(us) = %8ld %8ld, pts/us = %3ld %3ld\n", 
           nbits, NITEMS, errors, t1-t0, t2-t1, NITEMS/(t1-t0), NITEMS/(t2-t1)) ;
  }
}
