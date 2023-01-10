/* Hopefully useful routines for C and FORTRAN
 * Copyright (C) 2020-2023  Recherche en Prevision Numerique
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * Author : M. Valin (RPN-SI)
 */

#include <stdint.h>
#include <rmn/misc_encode.h>
#include <immintrin.h>

uint32_t stream_get_block_8x8(uint32_t * restrict src, int lni, uint32_t * restrict block, uint32_t * restrict pop, uint32_t * restrict gain ){
  int i, j, nbits, nbits0 ;
  uint32_t max = 0 ;
  uint32_t *block0 = block ;

  for(j=0 ; j<8 ; j++){
    for(i=0 ; i<8 ; i++) block[i] = src[i] ;
    src += lni ;
    block += 8 ;
  }
  for(i=0 ; i<64 ; i++) max = (block0[i] > max) ? block0[i] : max ; 
  nbits = 32 - _lzcnt_u32(max) ;
  for(i=0 ; i<33 ; i++) pop[i] = 0 ;
  for(i=0 ; i<64 ; i++) pop[32-_lzcnt_u32(block0[i])]++ ;
  for(i=1 ; i<33 ; i++) pop[i] = pop[i] + pop[i-1] ;
  nbits0 = (nbits+1) / 2 ;
  return nbits ;
}

void stream_encode_init(bitstream *bstream, void *buffer, size_t bufsize){
  BeStreamInit(bstream, buffer, bufsize, BIT_INSERT_MODE) ;
}

// encode a block of nx * ny unsigned integers, representable using nbits bits into stream
uint32_t stream_encode_ublock(uint32_t *src, int nx, int ny, int nbits, bitstream *bstream){
  uint64_t accum   = bstream->accum ;          // 64 bit unsigned accumulator
  uint32_t *stream = bstream->stream ;         // pointer to 32 bit unsigned int packed stream
  uint32_t insert  = bstream->insert ;         // number of unused bits in a
  int nbits0 = (nbits + 1) >> 1 ;              // length of "short" tokens
  uint32_t mask0 = (~0) << nbits0 ;            // mask to detect "short" tokens
  uint32_t flag = 1 << nbits ;                 // flag for full length tokens
  uint32_t header, header_length, total_length ;
  uint32_t token ;
  int token_length ;

  header = nbits ; header_length = 5 ;
  if(nx == 8 & ny == 8){                // full size block
    header <<= 1 ; header |= 1 ; header_length += 1 ;
  }else{
    header <<= 4 ; header |= nx ;
    header <<= 3 ; header |= ny ;
    header_length += 7 ;
  }
  header <<= 1 ; header |= 1 ; header_length += 1 ;  // add short/long encoding flag
  BE64_PUT_NBITS(accum, insert, header, nbits, stream) ;
  total_length = header_length ;

  while(ny-- >0){
    while(nx-- > 0){
      token = *src ; src++ ;
      if(token & mask0){           // not a "short" token"
        token |= flag ; token_length = nbits + 1 ;
      }else{
        token_length = nbits0 ;    // "short" token"
      }
      total_length += token_length ;
      BE64_PUT_NBITS(accum, insert, token, token_length, stream) ;
    }
  }

  bstream->accum  = accum ;       // update accumulator
  bstream->stream = stream ;      // update pointer into packed stream
  bstream->insert = insert ;      // update free bit count
  return total_length ;
}

// decode a block of nx * ny unsigned integers, representable using nbits bits into stream
uint32_t stream_decode_ublock(uint32_t *dst, bitstream *bstream){
  uint64_t accum   = bstream->accum ;          // 64 bit unsigned accumulator
  uint32_t *stream = bstream->stream ;         // pointer to 32 bit unsigned int packed stream
  uint32_t xtract  = bstream->xtract ;         // number of available bits in a
  uint32_t nbits, nbits0 ;
  uint32_t tmp, nx, ny ;

  BE64_GET_NBITS(accum, xtract, nbits, 5, stream) ;
  nbits0 = (nbits + 1) >> 1 ;              // length of "short" tokens
  BE64_GET_NBITS(accum, xtract, tmp, 1, stream) ;
  if(tmp == 1){         // full size block
    nx = 8 ; ny = 8 ;
  }else{
    BE64_GET_NBITS(accum, xtract, nx, 3, stream) ;
    BE64_GET_NBITS(accum, xtract, ny, 3, stream) ;
  }
  BE64_GET_NBITS(accum, xtract, tmp, 1, stream) ;  // get short/long encoding flag
  if(tmp == 1){
    while(ny-- > 0){
      while(nx-- > 0){
        BE64_GET_NBITS(accum, xtract, tmp, 1, stream) ;
        if(tmp == 1){                               // "long" token
          BE64_GET_NBITS(accum, xtract, dst[0], nbits, stream) ;
        }else{                                      // "short" token
          BE64_GET_NBITS(accum, xtract, dst[0], nbits0, stream) ;
        }
      }
    }
  }else{
    nx = 0 ; ny = 0 ;    // error if short/long encoding flag is 0
  }

  bstream->accum  = accum ;       // update accumulator
  bstream->stream = stream ;      // update pointer into packed stream
  bstream->xtract = xtract ;      // update free bit count

  return (nx << 3) | ny ;         // return block dimensions
}

void stream_encode_decode_test(){
  bitstream  thebits ;
}
