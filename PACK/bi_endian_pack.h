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
//     M. Valin,   Recherche en Prevision Numerique, 2022/09
//
#include <stdint.h>

#if ! defined(STATIC)
#define STATIC static
#endif

typedef struct{
  uint64_t  accum ;   // 64 bit unsigned bit accumulator
  uint32_t *start ;   // pointer to start of stream data storage
  uint32_t *stream ;  // pointer into packed stream
  uint32_t  insert ;  // insertion point into accumulator (number of bits already in accumulator)
  uint32_t  xtract ;  // extraction point from accumulator (number of bits available in accumulator)
} bitstream ;
//
// little endian style (right to left) bit stream packing
//
#define LE64_INSERT_BEGIN(accum, insert) \
        { accum = 0 ; insert = 0 ; }
#define LE64_INSERT_NBITS(accum, insert, w32, nbits) \
        { uint32_t mask = ~0 ; mask >>= (32-(nbits)) ; uint64_t w64 = (w32) & mask ; accum |= (w64 << insert) ; insert += (nbits) ; }
#define LE64_INSERT_CHECK(accum, insert, stream) \
        { if(insert > 32) { *stream = accum ; stream++ ; insert -= 32 ; accum >>= 32 ; } ; }
#define LE64_INSERT_FINAL(accum, insert, stream) \
        { LE64_INSERT_CHECK(accum, insert, stream) ; { if(insert > 0) { *stream = accum ; stream++ ;} ; } }
#define LE64_PUT_NBITS(accum, insert, w32, nbits) \
        { LE64_INSERT_CHECK(accum, insert, stream) ; LE64_INSERT_NBITS(accum, insert, w32, nbits) ; }

#define LE64_XTRACT_BEGIN(accum, xtract, stream) \
        { accum = 0 ; xtract = 0 ; }
#define LE64_XTRACT_NBITS(accum, xtract, w32, nbits) \
        { uint32_t mask = ~0 ; mask >>= (32-(nbits)) ; w32 = accum & mask ; accum >>= nbits ; xtract -= (nbits) ;}
#define LE64_XTRACT_CHECK(accum, xtract, stream) \
        { if(xtract < 32) { uint64_t w64 = *(stream) ; accum |= (w64 << xtract) ; (stream)++ ; xtract += 32 ; } ; }
#define LE64_XTRACT_FINAL(accum, xtract) \
        { accum = 0 ; xtract = 0 ; }
#define LE64_GET_NBITS(accum, xtract, w32, nbits) \
        { LE64_XTRACT_CHECK(accum, xtract, stream) ; LE64_XTRACT_NBITS(accum, xtract, w32, nbits) ; }
//
// big endian style (left to right) bit stream packing
//
#define BE64_INSERT_BEGIN(accum, insert) \
        { accum = 0 ; insert = 0 ; }
#define BE64_INSERT_NBITS(accum, insert, w32, nbits) \
        {  uint32_t mask = ~0 ; mask >>= (32-(nbits)) ; accum <<= (nbits) ; insert += (nbits) ; accum |= ((w32) & mask) ; }
#define BE64_INSERT_CHECK(accum, insert, stream) \
        { if(insert > 32) { insert -= 32 ; *(stream) = accum >> insert ; (stream)++ ; } ; }
#define BE64_INSERT_FINAL(accum, insert, stream) \
        { BE64_INSERT_CHECK(accum, insert, stream) ; if(insert > 0) { *stream = accum << (32 - insert) ; stream++ ; } }
#define BE64_PUT_NBITS(accum, insert, w32, nbits) \
        { BE64_INSERT_CHECK(accum, insert, stream) ; BE64_INSERT_NBITS(accum, insert, w32, nbits) ; }

#define BE64_XTRACT_BEGIN(accum, xtract, stream) \
        { accum = 0 ; xtract = 0 ; }
#define BE64_XTRACT_NBITS(accum, xtract, w32, nbits) \
        { w32 = accum >> (64 - (nbits)) ; accum <<= (nbits) ; xtract -= (nbits) ; }
#define BE64_XTRACT_CHECK(accum, xtract, stream) \
        { if(xtract < 32) { accum >>= (32-xtract) ; accum |= *stream ; accum <<= (32-xtract) ; xtract += 32 ; (stream)++ ; } ; }
#define BE64_XTRACT_FINAL(accum, xtract) \
        { accum = 0 ; xtract = 0 ; }
#define BE64_GET_NBITS(accum, xtract, w32, nbits) \
        { BE64_XTRACT_CHECK(accum, xtract, stream) ; BE64_XTRACT_NBITS(accum, xtract, w32, nbits) ; }

STATIC inline void  LeStreamInit(bitstream *p, uint32_t *buffer){
  p->accum  = 0 ;         // accumulator is empty
  p->insert = 0 ;         // insertion point at Least Significant Bit
  p->xtract = 0 ;         // extraction point at Least Significant Bit
  p->start  = buffer ;    // stream storage
  p->stream = buffer ;    // stream is empty and starts at 
}

STATIC inline void  BeStreamInit(bitstream *p, uint32_t *buffer){
  p->accum  = 0 ;         // accumulator is empty
  p->insert = 0 ;         // no data has been inserted
  p->xtract = 0 ;         // no data available for extraction
  p->start  = buffer ;    // stream storage
  p->stream = buffer ;    // stream is empty and starts at 
}

void  LeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw);
void  LeStreamInsertM(bitstream *p, uint32_t *w32, int *nbits, int *n);

void  BeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw);
void  BeStreamInsertM(bitstream *p, uint32_t *w32, int *nbits, int *n);

void  LeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n);
void  LeStreamXtractM(bitstream *p, uint32_t *w32, int *nbits, int *n);

void  BeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n);
void  BeStreamXtractM(bitstream *p, uint32_t *w32, int *nbits, int *n);

STATIC inline uint32_t MaskNbits1(uint32_t nbits){
  uint32_t mask = ~0 ;
  return ~(mask << nbits) ;
}

STATIC inline uint32_t MaskNbits2(uint32_t nbits){
  uint32_t mask = ~0 ;
  return  ( mask >> (32 - nbits)) ;
}

STATIC inline void  BeStreamReset(bitstream *p, uint32_t *buffer){
  p->accum = 0 ;
  p->insert = 0 ;
  p->stream = buffer ;
  p->start  = buffer ;
}
