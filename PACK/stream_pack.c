#include <stdint.h>

// insert nbits from right into acc64 assuming that there is room enough for up to 32 bits
#define in64(acc64,nbits,nleft,what) { acc64 = (acc64<<nbits) | (what) ; nleft -= nbits ; }

// make sure that there is room enough in acc64 for up to 32 bits, write 32 bits into where if need be
#define chk64(acc64,nleft,where) { if (nleft<0) { *++where = acc64 >> (-nleft) ; nleft += 32 ; } }

//align acc64 and flush to memory assuming at most 32 useful bits in acc64
#define flush64(acc64,nleft,where) { if (nleft<32) { *++where = acc64 << nleft ; nleft += 32 ; } }

// insert nbits from right into acc64 (no asumptions)
#define insert64(acc64,nbits,nleft,what,where) { in64(acc64,nbits,nleft,what) ; chk64(acc64,nleft,where) ; }

// get leftmost nbits from acc64 assuming that there are at least 32 useful bits available
#define out64(acc64,nbits,nleft,where) { (where) = (acc64 >> (64-nbits) ) ; acc64 <<= nbits ; nleft -= 32 ; }

// make sure that  there are at least 32 useful bits available in acc64, get 32 bits from FROM if need be
#define chk64out(acc64,nleft,FROM) { if (nleft<0) { acc64 >>= (32-nleft); acc64 += (FROM)  ; acc64 <<= (32-nleft) ; nleft += 32 ;} }

void stream_unpack_16(uint16_t *where, uint32_t *from, uint32_t nbits, uint32_t n)
{
  uint64_t acc64;
  int32_t nleft=32;
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

int stream_pack_16(uint16_t *what, uint32_t *where, uint32_t nbits, uint32_t n)
{
  uint64_t acc64=0;
  int32_t nleft=32;
  int i;
  uint32_t *where_start = where;
  uint32_t nb2 = nbits + nbits;

  for(i=0 ; i<n-1 ; i+=2){
    in64(acc64, nbits, nleft, *what++) ;  // nbits assumed <= 16
    in64(acc64, nbits, nleft, *what++) ;
    chk64(acc64, nleft, where) ;
  }
  if(n & 1) {
    in64(acc64,nbits,nleft,*what) ;
    chk64(acc64,nleft,where) ;
  }
  flush64(acc64, nleft, where) ;
  return where_start - where ;
}
