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

#if ! defined(MAKE_SIGNED_32)

#if ! defined(STATIC)
#define STATIC static
#define STATIC_DEFINED_HERE
#endif

// bit stream descriptor. ONLY ONE of insert / extract should be non zero
// in insertion mode, extract should be 0
// in extraction mode, insert should be 0
typedef struct{
  uint64_t  accum ;   // 64 bit unsigned bit accumulator
  uint32_t *start ;   // pointer to start of stream data storage
  uint32_t *stream ;  // pointer into packed stream (both insert and extract mode)
  uint32_t  insert ;  // # of bits used in accumulator (0 <= insert <= 64)
  uint32_t  xtract ;  // # of bits extractable from accumulator (0 <= xtract <= 64)
} bitstream ;

// for this macro to produce meaningful results, w32 MUST BE int32_t (32 bit signed int)
#define MAKE_SIGNED_32(w32, nbits) { w32 <<= (32 - (nbits)) ; w32 >>= (32 - (nbits)) ; }
// for this macro to produce meaningful results, w64 MUST BE int64_t (64 bit signed int)
#define MAKE_SIGNED_64(w64, nbits) { w64 <<= (64 - (nbits)) ; w64 >>= (64 - (nbits)) ; }

//
// macro arguments description
// accum  [INOUT] : 64 bit accumulator
// insert [INOUT] : # of bits used in accumulator (0 <= insert <= 64)
// xtract [INOUT] : # of bits extractable from accumulator (0 <= xtract <= 64)
// stream [INOUT] : pointer to next position in packed stream
// w32    [IN]    : 32 bit integer containing data to be inserted (expression allowed)
//        [OUT]   : 32 bit integer receiving extracted data (MUST be a variable)
// nbits  [IN]    : number of bits to insert / extract in w32 (<= 32 bits)
//
// N.B. : if w32 and accum are signed variables, extraction will produce a "signed" result
//        if w32 and accum are unsigned variables, extraction will produce an "unsigned" result
//
// little endian style (right to left) bit stream packing
// insertion at top, most significant part
// the useful bits are at the bottom (least significant part) of accum
//
// initialize stream for insertion
#define LE64_INSERT_BEGIN(accum, insert) \
        { accum = 0 ; insert = 0 ; }
// insert the lower nbits bits from w32 into accum, update insert, accum
#define LE64_INSERT_NBITS(accum, insert, w32, nbits) \
        { uint32_t mask = ~0 ; mask >>= (32-(nbits)) ; uint64_t w64 = (w32) & mask ; accum |= (w64 << insert) ; insert += (nbits) ; }
// check that 32 bits can be safely inserted into accum
// if not possible, store lower 32 bits of accum into stream, update accum, insert, stream
#define LE64_INSERT_CHECK(accum, insert, stream) \
        { if(insert > 32) { *stream = accum ; stream++ ; insert -= 32 ; accum >>= 32 ; } ; }
// store any residual data from accum into stream, update accum, insert, stream
#define LE64_INSERT_FINAL(accum, insert, stream) \
        { LE64_INSERT_CHECK(accum, insert, stream) ; { if(insert > 0) { *stream = accum ; stream++ ;} ; } }
// combined INSERT_CHECK and INSERT_NBITS, update accum, insert, stream
#define LE64_PUT_NBITS(accum, insert, w32, nbits, stream) \
        { LE64_INSERT_CHECK(accum, insert, stream) ; LE64_INSERT_NBITS(accum, insert, w32, nbits) ; }

// N.B. : if w32 and accum are signed variables, the extract will produce a "signed" result
//        if w32 and accum are unsigned variables, the extract will produce an "unsigned" result
// initialize stream for extraction
#define LE64_XTRACT_BEGIN(accum, xtract, stream) \
        { accum = 0 ; xtract = 0 ; }
// extract nbits bits into w32 from accum, update xtract, accum
#define LE64_XTRACT_NBITS(accum, xtract, w32, nbits) \
        { w32 = (accum << (64-nbits)) >> (64-nbits) ; accum >>= nbits ; xtract -= (nbits) ;}
// check that 32 bits can be safely extracted from accum
// if not possible, get extra 32 bits into accum from stresm, update accum, xtract, stream
#define LE64_XTRACT_CHECK(accum, xtract, stream) \
        { if(xtract < 32) { uint64_t w64 = *(stream) ; accum |= (w64 << xtract) ; (stream)++ ; xtract += 32 ; } ; }
// finalize extraction, update accum, xtract
#define LE64_XTRACT_FINAL(accum, xtract) \
        { accum = 0 ; xtract = 0 ; }
// combined XTRACT_CHECK and XTRACT_NBITS, update accum, xtract, stream
#define LE64_GET_NBITS(accum, xtract, w32, nbits, stream) \
        { LE64_XTRACT_CHECK(accum, xtract, stream) ; LE64_XTRACT_NBITS(accum, xtract, w32, nbits) ; }
//
// big endian style (left to right) bit stream packing
// insertion at bottom, least significant part
//
// initialize stream for insertion
#define BE64_INSERT_BEGIN(accum, insert) \
        { accum = 0 ; insert = 0 ; }
// insert the lower nbits bits from w32 into accum, update insert, accum
#define BE64_INSERT_NBITS(accum, insert, w32, nbits) \
        {  uint32_t mask = ~0 ; mask >>= (32-(nbits)) ; accum <<= (nbits) ; insert += (nbits) ; accum |= ((w32) & mask) ; }
// check that 32 bits can be safely inserted into accum
// if not possible, store lower 32 bits of accum into stream, update accum, insert, stream
#define BE64_INSERT_CHECK(accum, insert, stream) \
        { if(insert > 32) { insert -= 32 ; *(stream) = accum >> insert ; (stream)++ ; } ; }
// store any residual data from accum into stream, update accum, insert, stream
#define BE64_INSERT_FINAL(accum, insert, stream) \
        { BE64_INSERT_CHECK(accum, insert, stream) ; if(insert > 0) { *stream = accum << (32 - insert) ; stream++ ; } }
// combined INSERT_CHECK and INSERT_NBITS, update accum, insert, stream
#define BE64_PUT_NBITS(accum, insert, w32, nbits, stream) \
        { BE64_INSERT_CHECK(accum, insert, stream) ; BE64_INSERT_NBITS(accum, insert, w32, nbits) ; }

// N.B. : if w32 and accum are signed variables, the extract will produce a "signed" result
//        if w32 and accum are unsigned variables, the extract will produce an "unsigned" result
// initialize stream for extraction
#define BE64_XTRACT_BEGIN(accum, xtract, stream) \
        { accum = 0 ; xtract = 0 ; }
// extract nbits bits into w32 from accum, update xtract, accum
#define BE64_XTRACT_NBITS(accum, xtract, w32, nbits) \
        { w32 = accum >> (64 - (nbits)) ; accum <<= (nbits) ; xtract -= (nbits) ; }
// check that 32 bits can be safely extracted from accum
// if not possible, get extra 32 bits into accum from stresm, update accum, xtract, stream
#define BE64_XTRACT_CHECK(accum, xtract, stream) \
        { if(xtract < 32) { accum >>= (32-xtract) ; accum |= *stream ; accum <<= (32-xtract) ; xtract += 32 ; (stream)++ ; } ; }
// finalize extraction, update accum, xtract
#define BE64_XTRACT_FINAL(accum, xtract) \
        { accum = 0 ; xtract = 0 ; }
// combined XTRACT_CHECK and XTRACT_NBITS, update accum, xtract, stream
#define BE64_GET_NBITS(accum, xtract, w32, nbits, stream) \
        { BE64_XTRACT_CHECK(accum, xtract, stream) ; BE64_XTRACT_NBITS(accum, xtract, w32, nbits) ; }

// initialize a LittleEndian stream
STATIC inline void  LeStreamInit(bitstream *p, uint32_t *buffer){
  p->accum  = 0 ;         // accumulator is empty
  p->insert = 0 ;         // insertion point at Least Significant Bit
  p->xtract = 0 ;         // extraction point at Least Significant Bit
  p->start  = buffer ;    // stream storage
  p->stream = buffer ;    // stream is empty and starts at 
}

// initialize a BigEndian stream
STATIC inline void  BeStreamInit(bitstream *p, uint32_t *buffer){
  p->accum  = 0 ;         // accumulator is empty
  p->insert = 0 ;         // no data has been inserted
  p->xtract = 0 ;         // no data available for extraction
  p->start  = buffer ;    // stream storage
  p->stream = buffer ;    // stream is empty and starts at 
}

// insert multiple values (unsigned)
void  LeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw);
// insert multiple values from list (unsigned)
void  LeStreamInsertM(bitstream *p, uint32_t *w32, int *nbits, int *n);

// insert multiple values (unsigned)
void  BeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw);
// insert multiple values from list (unsigned)
void  BeStreamInsertM(bitstream *p, uint32_t *w32, int *nbits, int *n);

// extract multiple values (unsigned)
void  LeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n);
// extract multiple values (signed)
void  LeStreamXtractSigned(bitstream *p, int32_t *w32, int nbits, int n);
// extract multiple values from list (unsigned)
void  LeStreamXtractM(bitstream *p, uint32_t *w32, int *nbits, int *n);

// extract multiple values (unsigned)
void  BeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n);
// extract multiple values (signed)
void  BeStreamXtractSigned(bitstream *p, int32_t *w32, int nbits, int n);
// extract multiple values from list (unsigned)
void  BeStreamXtractM(bitstream *p, uint32_t *w32, int *nbits, int *n);

STATIC inline uint32_t RMask(uint32_t nbits){
  uint32_t mask = ~0 ;
  return  ( mask >> (32 - nbits)) ;
}

STATIC inline uint32_t LMask(uint32_t nbits){
  uint32_t mask = ~0 ;
  return  ( mask << (32 - nbits)) ;
}

// flush stream being written into
STATIC inline void  LeStreamFlush(bitstream *p){
  if(p->insert > 0) LE64_INSERT_FINAL(p->accum, p->insert, p->stream) ;
  p->accum = 0 ;
  p->insert = 0 ;
  p->xtract = 0 ;
}

// flush stream being written into
STATIC inline void  BeStreamFlush(bitstream *p){
  if(p->insert > 0) BE64_INSERT_FINAL(p->accum, p->insert, p->stream) ;
  p->accum = 0 ;
  p->insert = 0 ;
}

// rewind stream to read it from the beginning
STATIC inline void  LeStreamRewind(bitstream *p){
  if(p->insert > 0) LeStreamFlush(p) ;   // something left to write
  p->accum = 0 ;
  p->insert = 0 ;
  p->xtract = 0 ;
  p->stream = p->start ;
}

// rewind stream to read it from the beginning
STATIC inline void  BeStreamRewind(bitstream *p){
  if(p->insert > 0) BeStreamFlush(p) ;   // something left to write
  p->accum = 0 ;
  p->insert = 0 ;
  p->xtract = 0 ;
  p->stream = p->start ;
}

// STATIC inline void  LeStreamReset(bitstream *p, uint32_t *buffer){
//   p->accum = 0 ;
//   p->insert = 0 ;
//   p->xtract = 0 ;
//   p->stream = buffer ;
//   p->start  = buffer ;
// }
// 
// STATIC inline void  BeStreamReset(bitstream *p, uint32_t *buffer){
//   p->accum = 0 ;
//   p->insert = 0 ;
//   p->xtract = 0 ;
//   p->stream = buffer ;
//   p->start  = buffer ;
// }

#if defined(STATIC_DEFINED_HERE)
#undef STATIC
#undef STATIC_DEFINED_HERE
#endif

#endif
