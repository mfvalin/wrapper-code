#include <stdio.h>
#include <stdint.h>

#include <rmn/misc_encode.h>

int main(int argc, char **argv){
  int i, nbits ;
  uint32_t src[64], block[64] ;
  for(i=0 ; i<64 ; i++) src[i] = 2*i + 125 ;
  printf("src[63] = %d, nbits = %d\n", src[63], nbits) ;
  nbits = stream_get_block_8x8(src, 8, block) ;
  printf("src[63] = %d, nbits = %d\n", src[63], nbits) ;
}

