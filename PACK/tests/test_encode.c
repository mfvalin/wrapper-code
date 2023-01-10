#include <stdio.h>
#include <stdint.h>

#include <rmn/misc_encode.h>

int main(int argc, char **argv){
  int i, nbits ;
  uint32_t src[64], block[64], pop[34], gain[3] ;

  for(i=0 ; i<64 ; i++) src[i] = 5*i ;
  nbits = stream_get_block_8x8(src, 8, block, pop, gain) ;
  printf("src[63] = %d, nbits = %d\n", src[63], nbits) ;
  for(i=0 ; i<16 ; i++) printf("%3d", pop[i]); printf("\n") ;
}

