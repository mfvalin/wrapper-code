/*
 * Copyright (C) 2021  Environnement Canada
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * Author:
 *     M. Valin,   Recherche en Prevision Numerique, 2021
 */
#if ! defined(PACK_STREAM_FUNCTIONS)
#define PACK_STREAM_FUNCTIONS

#include <stdint.h>

// inline functions to initialize, insert into, extract from a bit stream
// said packed stream is a sequence of 32 bit unsigned integers
// the pstream structure MUST be initialized with pstream_init before any other call

typedef struct{   // packing stream (the packed stream is a sequence of 32 bit unsigned integers)
  uint64_t a ;    // 64bit accumulator
  uint32_t *s ;   // start of stream buffer
  uint32_t *p ;   // current insertion/extraction pointer into stream buffer
  int32_t nf ;    // number of bits free for insertion
  int32_t na ;    // number of bits available for extraction
} pstream ;

void inline pstream_init(pstream *ps, void *buffer){  // initialize a stream for read or write
  ps->a = 0 ;                     // not really necessary
  ps->s = (uint32_t *) buffer ;   // saved address of user buffer
  ps->p = (uint32_t *) buffer ;   // current address into user buffer
  ps->nf = 64 ;                   // 64 bits available for insertion (put)
  ps->na = 0 ;                    // 0 bits available for extraction (get)
}

void inline pstream_put_32(pstream *ps, uint32_t token, int nbits){   // insert nbits (<= 32) bits into stream
  ps->a <<= nbits ;     // make space for token
  ps->a  |= token ;     // add token
  ps->nf -= nbits ;     // bump down free bits counter
  if(ps->nf < 32) {     // make sure that at least 32 bits can always be inserted without problem
    *(ps->p) = ps->a >> (32 - (ps->nf)) ;   // write 32 leftmost used bits from accumulator
    ps->p = ps->p + 1 ;                     // bump write pointer
    ps->nf = ps->nf + 32 ;                  // bump free space counter
  }
}

void inline pstream_put_64(pstream *ps, uint64_t token, int nbits){   // insert nbits (<= 64) bits into stream
  uint32_t token32 ;
  if(nbits > 32) {                               // more than 32 bits to write
    token32 = token >> 32 ;
    pstream_put_32(ps, token32 , nbits - 32) ;   // insert upper nbits - 32 bits
    token32 = token & 0xFFFFFFFFL ;              // lower 32 bits of token
    pstream_put_32(ps, token32 , 32) ;           // insert lower 32 bits
  }else{                                         // 32 bits or less to write
    token32 = token & 0xFFFFFFFFL ;              // lower 32 bits of token
    pstream_put_32(ps, token32 , nbits) ;
  }
}

void inline pstream_flush(pstream *ps){    // write whatever is left in the accumulator into the stream buffer
  uint32_t token32 ;                       // forget whatever was left in the accumulator for a read stream
  if(ps->nf < 64) {                        // NO-OP if read stream
    token32 = ps->a << (ps->nf - 32) ;     // align letfovers 
    *(ps->p) = token32 ;                   // store leftovers
    ps->p = ps->p + 1 ;                    // bump store pointer
  }
  ps->a  = 0 ;                             // accumulator is zeroed
  ps->na = 0 ;                             // no available bits
  ps->nf = 64 ;                            // 64 free bits (accumulator is empty)
}

void inline pstream_rewind(pstream *ps){   // bring a read or write stream back to the beginning
  pstream_flush(ps) ;     // in case it is a write stream
  ps->na = 0 ;            // reset read counter to empty
  ps->a  = 0;
  ps->p = ps->s ;         // set current pointer back to beginning of buffer
}

uint32_t  inline pstream_get_32(pstream *ps, int nbits){   // extract nbits (<= 32) bits from stream
  uint32_t t32 ;
  uint64_t t64 ;
  if(ps->na < 32){                     // always keep at least 32 bits available in read buffer
    ps->a = ps->a >> (32 - ps->na) ;   // align usefule bits at halfway position
    t64 = *(ps->p) ;                   // get 32 bits from stream
    ps->p = ps->p + 1 ;                // bump read pointer
    ps->na = ps->na + 32 ;             // bump bits available counter
    ps->a = ps->a | t64 ;              // insert fresh 32 bits
    ps->a = ps->a << (64 - ps->na) ;   // realign left
  }
  ps->na = ps->na - nbits ;            // bump down bits available counter
  t32 = (ps->a) >> (64 - nbits) ;      // extract leftmost nbits bits
  ps->a = (ps->a) << nbits ;           // keep bits to be extracted in the most significant position
  return t32 ;
}

uint64_t  inline pstream_get_64(pstream *ps, int nbits){   // extract nbits (> 32) bits from stream
  uint32_t t32 ;
  uint64_t t64 ;
  if(nbits <= 32){                            // <= 32 bits to extract
    t64 = pstream_get_32(ps, nbits) ;
  }else{                                      // more than 32 bits to extract
    t64 = pstream_get_32(ps, nbits - 32) ;    // upper nbits - 32 bits
    t32 = pstream_get_32(ps, 32) ;            // lower 32 bits
    t64 = (t64 << 32) | t32 ;                 // merge both parts
  }
  return t64 ;
}

#endif
