// Hopefully useful code for C
// Copyright (C) 2022  Recherche en Prevision Numerique
//
// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// set of functions to handle insertion/extraction of 1-32 bit tokens into/from a stream of 32 bit unsigned integers
// should a token be langer thatn 32 bits, it must be split into smaller tokens before insertion/extraction
//
#if ! defined(MISC_PACK)
#define MISC_PACK

#include <stdint.h>
#include <stddef.h>

// the basic struct to handle insertion/extraction into/from a packed stream
typedef struct{
  uint64_t    a ;  // 64 bit unsigned accumulator
  uint32_t   *s ;  // pointer to 32 bit unsigned int packed stream
  uint32_t used ;  // number of bits used in a (used by unpacker)
  uint32_t free ;  // number of unused bits in a (used by packer)
} stream32 ;

// initialize the pack/unpack stream control struct
static inline void stream32_init(stream32 *p, void *stream){
  p->free = 64 ;       // 64 unused bit positions in accumulator
  p->used = 0 ;        // no used bit
  p->a    = 0 ;        // accumulator is empty
  p->s    = stream ;   // point to beginning of stream
}

// store whatever is left in the accumulator into the packed stream
static inline void stream32_flush(stream32 *p){
  if(p->free < 32){
    *p->s = p->a >> (32 -p->free) ;
    p->s++ ;
    p->free += 32 ;
  }
  if(p->free < 64){
    p->a <<= (p->free - 32) ;
    *p->s = p->a ;
  }
}

// insert 1 item (<= 32 bits) assuming that enough space is available
static inline void stream32_put_fast(stream32 *p, uint32_t item, int nbits){
  p->a <<= nbits ;
  p->a |= item ;
  p->free = p->free - nbits ;
}

// check that at least 32 bits are available to insert an item
static inline void stream32_put_check(stream32 *p){
  if(p->free >= 32) return ;
  *p->s = (p->a >> (32 - p->free) ) ;   // store upper 32 bits
  p->s++ ;                              // bump pointer topacked stream
  p->free = p->free + 32 ;              // adjust free bit count
}

// insert 1 item (<= 32 bits) with check for space availibility in accumulator
static inline void stream32_put(stream32 *p, uint32_t item, int nbits){
  int full = (p->free < 32) ;                    // less than 32 positions available ?
  *p->s = (p->a >> (32 - p->free) ) ;            // store upper 32 bits, in case it is "full"
  p->a <<= nbits ;                               // insert item into accumulator
  p->a |= item ;
  p->free = p->free - nbits + (full ? 32 : 0 ) ; // adjust count, add 32 if "full"
  p->s += (full ? 1 : 0 ) ;                      // bump store stream pointer if "full"
}

// pack n items, nbits (<=32) bits long into a stream
static void pack_stream(void *p, uint32_t *u, int nbits, int n){
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

// extract 1 item (<= 32 bits) assuming that enough data is available in accumulator
static inline uint32_t stream32_get_fast(stream32 *p, int nbits){
  uint32_t item = p->a >> (64 - nbits) ;
  p->a <<= nbits ;
  p->used = p->used - nbits ;
  return item ;
}

// check that at least 32 bits are available to extract an item
static inline void stream32_get_check(stream32 *p){
  uint64_t temp ;
  if(p->used >= 32) return ;
  temp = *p->s ;                        // get new data from stream
  p->a = p->a | ( temp << (32 - p->used)) ;
  p->s++ ;                              // bump pointer topacked stream
  p->used = p->used + 32 ;              // adjust free bit count
}

// extract 1 item (<= 32 bits) with check for data availibility in accumulator
static inline uint32_t stream32_get(stream32 *p, int nbits){
  uint32_t item ;
  int empty = (p->used < 32) ;                    // less than 32 bits available ?
  uint64_t temp = *p->s ;                         // fetch new data just in case
  p->a = empty ? (p->a | (temp << (32 - p->used))) : p->a ;
  item = (p->a >> (64 - nbits)) ;                 // extract item
  p->a <<= nbits ;                                // flush item from accumulator
  p->used = p->used - nbits + (empty ? 32 : 0 ) ; // adjust count, add 32 if "empty"
  p->s += (empty ? 1 : 0 ) ;                      // bump fetch stream pointer if "empty"
  return item ;
}

// extract n items nbits(<= 32) bits long from a stream
static void unpack_stream(void *p, uint32_t *u, int nbits, int n){
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

#endif
