#include <stdint.h>

// insert nbits from right into acc64 assuming that there is room space for up to 32 bits
#define in64(acc64,nbits,nleft,src) { acc64 = (acc64<<nbits) | (src) ; nleft -= nbits ; }

// insert nbits from right into acc64 (no assumptions about available space)
#define chk64in(acc64,nbits,nleft,src,dest) {uint64_t t=acc64>>(32-nleft); acc64<<=(nbits); acc64|=(src); if(nleft<32){nleft += 32; (dest)=t;}; nleft-=nbits;}

// make sure that there is room enough in acc64 for up to 32 bits, write 32 bits into where if need be
#define chk64(acc64,nleft,where) { if (nleft<0) { *++where = acc64 >> (-nleft) ; nleft += 32 ; } }

//align acc64 and flush to memory assuming at most 32 useful bits in acc64
#define flush64(acc64,nleft,dest) { if (nleft<32){ (dest)=acc64>>(32-nleft); nleft+=32; } ; if(nleft<64){ acc64<<=nleft; (dest)=acc64>>32 ;}  }

// get leftmost nbits from acc64 assuming that there are at least 32 useful bits available
#define out64(acc64,nbits,nleft,where) { (where) = (acc64 >> (64-nbits) ) ; acc64 <<= nbits ; nleft -= 32 ; }

// make sure that  there are at least 32 useful bits available in acc64, get 32 bits from FROM if need be
#define chk64out(acc64,nleft,FROM) { if (nleft<0) { acc64 >>= (32-nleft); acc64 += (FROM)  ; acc64 <<= (32-nleft) ; nleft += 32 ;} }

uint64_t rdtscp(uint32_t *cpu) {   // version "in order" avec "serialization"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi, c;
  __asm__ volatile ("rdtscp"
      : /* outputs */ "=a" (lo), "=d" (hi), "=c" (c)
      : /* no inputs */
      : );
//       : /* clobbers */ "%rcx");
  __asm__ volatile ("mfence");
  *cpu = c;
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  return time0++;
#endif
}

uint64_t rdtsc(void) {   // version rapide "out of order"
#if defined(__x86_64__) || defined( __i386__ )
  uint32_t lo, hi;
  __asm__ volatile ("rdtsc"
      : /* outputs */ "=a" (lo), "=d" (hi)
      : /* no inputs */
      : /* clobbers */ "%rcx");
  return (uint64_t)lo | (((uint64_t)hi) << 32);
#else
  return time0++;
#endif
}

void stream_unpack_16(uint16_t *where, uint32_t *from, uint32_t nbits, uint32_t n)
{
  uint64_t acc64;
  int32_t nleft=64;
  int i;

  acc64 = *from++ ;                     // prime the pump with 64 bits 
  acc64 = (acc64 << 32) | *from++ ;

  for(i=0 ; i<n-1 ; i+=2){
    out64(acc64,nbits,nleft,*where++) ;  // nbits assumed <= 16
    out64(acc64,nbits,nleft,*where++) ;
    chk64out(acc64,nleft,*from) ;
  }
  if(n & 1) {
    *where = acc64 >> 64-nbits ;
  }
  
}

int stream_pack_16_2(uint16_t *src, uint32_t *dest, uint16_t *src_2, uint32_t *dest_2, uint32_t nbits, uint32_t n)
{
  uint64_t acc64=0;
  uint64_t acc64_2=0;
  int32_t nleft=64;
  int32_t nleft_2=64;
  int i;
  uint32_t *dest_start = dest;

  for(i=0 ; i<n-1 ; i+=2){
    chk64in(acc64  , nbits, nleft  , *src++  , *dest++  ) ;  // nbits assumed <= 16
    chk64in(acc64_2, nbits, nleft_2, *src_2++, *dest_2++) ;  // nbits assumed <= 16
    in64(acc64  , nbits, nleft  , *src++  ) ;
    in64(acc64_2, nbits, nleft_2, *src_2++) ;
  }
  while(i++ < n) {
    chk64in(acc64  , nbits, nleft  , *src++  , *dest++  ) ;
    chk64in(acc64_2, nbits, nleft_2, *src_2++, *dest_2++) ;
  }
  flush64(acc64  , nleft  , *dest++  ) ;
  flush64(acc64_2, nleft_2, *dest_2++) ;
  return dest - dest_start;
}

int stream_pack_16(uint16_t *src, uint32_t *dest, uint32_t nbits, uint32_t n)
{
  uint64_t acc64=0;
  int32_t nleft=64;
  int i;
  uint32_t *dest_start = dest;
  uint32_t nb2 = nbits + nbits;

  if(nbits <=1) {
    for(i=0 ; i<n-3 ; i+=4){
      chk64in(acc64, nbits, nleft, *src++, *dest++) ;  // nbits assumed <= 16
      in64(acc64, nbits, nleft, *src++) ;
      in64(acc64, nbits, nleft, *src++) ;
      in64(acc64, nbits, nleft, *src++) ;
    }
  }else{
    for(i=0 ; i<n-1 ; i+=2){
      chk64in(acc64, nbits, nleft, *src++, *dest++) ;  // nbits assumed <= 16
      in64(acc64, nbits, nleft, *src++) ;
    }
  }
  while(i++ < n) {
    chk64in(acc64, nbits, nleft, *src++, *dest++) ;
  }
  flush64(acc64, nleft, *dest++) ;
  return dest - dest_start;
}

int stream_pack_16x4(uint16_t *src, uint32_t *packed, uint32_t nbits, uint32_t n)
{
  uint64_t acc64[4], src64[4];
  int32_t nleft=64;
  int i, j;
  uint32_t *packed_start = packed;
  uint32_t nb2 = nbits + nbits;

  for(i=0 ; i<n-7 ; i+=8){
    for(j=0 ; j<4 ; j++) { src64[j] = src[j  ] ; }
    for(j=0 ; j<4 ; j++) { acc64[j] <<= nbits ; acc64[j] |= src64[j] ; } ;
    for(j=0 ; j<4 ; j++) { src64[j] = src[j+4] ; }
    for(j=0 ; j<4 ; j++) { acc64[j] <<= nbits ; acc64[j] |= src64[j] ; } ;
    nleft -= nb2;
    src += 8;
    if(nleft < 32) {
      for(j=0 ; j<4 ; j++) packed[j] = acc64[j] >> (32-nleft);
      packed += 4;
      nleft += 32;
    }
  }
  if(nleft < 64){
    for(j=0 ; j<4 ; j++) { acc64[j] <<= nleft ; packed[j] = acc64[j] >>32; }
  }
  nleft = 64;
  while(i++ < n) {
    chk64in(acc64[0], nbits, nleft, *src++, *packed++) ;
  }
  flush64(acc64[0], nleft, *packed++) ;
  return packed - packed_start;
}

#if defined(SELF_TEST)
uint64_t nsec(uint64_t ticks){
  double t0;
  t0 = ticks ; t0 /= 3.7 ; return( t0 + .5 ) ;
}
#include <stdio.h>
main(){
  uint32_t cpu;
  uint64_t t1, t2;
  double t0;
  uint16_t unpacked[1800];
  uint32_t packed[1800];
  uint16_t unpacked2[1800];
  uint32_t packed2[1800];
  int i, n, j;


  for(j=4 ; j<=16 ; j++) {
    for(i=0 ; i<1800 ; i++) { unpacked[i] = i ; packed[i] = -1 ;}
    for(i=0 ; i<1800 ; i++) { unpacked2[i] = i ; packed2[i] = -1 ;}
    t1 = rdtscp(&cpu);
    n = stream_pack_16x4(unpacked,packed,j,1800);
//     n = stream_pack_16_2(unpacked,packed,unpacked2,packed2,j,1800);
  //   sleep(1);
    t2 = rdtscp(&cpu) - t1;
    printf("j = %d,ticks = %ld, n = %d\n",j,nsec(t2),n);
  }
  for(i=0 ; i<10 ; i++) printf("%8.8x ",packed[i]); printf("\n");
}
#endif
