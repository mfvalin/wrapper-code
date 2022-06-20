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
//     M. Valin,   Recherche en Prevision Numerique, june 2022
// 

#if ! defined(PACK_MACROS_INCLUDE)
#define PACK_MACROS_INCLUDE

#include <stdint.h>

#if ! defined(RMASK)
// 32 and 64 bit left aligned masks
#define LMASK32(nbits)  ((~0)   << (32-nbits))
#define LMASK64(nbits)  ((~0l)  << (64-nbits))
// 32 and 64 bit left aligned masks
#define RMASK32(nbits)  (~((~0)  << nbits))
#define RMASK64(nbits)  (~((~0l) << nbits))
#endif

// macros to insert into, extract from a 32 bit token stream

// =============================== token insertion ===============================

// acc   : 64 bit accumulator (uint64_t)
// sp    : pointer into uint32_t array
// count : numbre of available/free bits in accumulator
// token : 32/64 bit unsigned integer (uint32_t/uint64_t)
// nbits : token length (the token is ASSUMED to only have 1s in the rightmost nbits bits)

// initialize bit count and accumulator to store into a fresh stream of 32 bit UNSIGNED elements
// an empty accumulator with 64 bits available for insertion
#define PACK32_IN_INIT(acc, count) {acc = 0 ; count = 64 ; }

// insert "short" variable length tokens into a stream of UNSIGNED 32 bit elements
// nbits <= 32, acc MUST be uint64_t, token MUST NOT have extra bits, sp MUST be a pointer to uint32_t

// UNCHECKED insertion, assuming that at least nbits bits are available for insertion in accumulator
#define PACK32_IN_FAST(acc,count,token,nbits) { acc <<= (nbits) ; acc |= (token) ; count -= (nbits) ; }

// CHECK to make sure that at least 32 bits are available for insertion
#define PACK32_IN_CHECK(sp,acc,count) if(count<32) { *sp = acc >> (32-count) ; sp ++ ; count += 32 ; }

// safe insertion, at least 32 bits are available before insertion
#define PACK32_IN(sp,acc,count,token,nbits) PACK32_IN_CHECK(sp,acc,count) ; PACK32_IN_FAST(acc,count,token,nbits) ;

// insert "long" variable length tokens into a stream of UNSIGNED 32 bit elements
// nbits <= 64, acc MUST be uint64_t, token MUST NOT have extra bits, sp MUST be a pointer to uint32_t
// insertion is done in 2 steps, nbits - 32 upper bits first, then the lower 32 bits
#define PACK32_IN_L(sp,acc,count,token,nbits) \
  if((nbits) > 32) { PACK32_IN(sp,acc,count,(token>>32),((nbits)-32)) ; PACK32_IN(sp,acc,count,(token)&0xFFFFFFFFl,32) ; } \
  else { PACK32_IN(sp,acc,count,(token)&0xFFFFFFFFl,nbits) }

// store the remaining partial token from the accumulator into the packed stream (align it left first)
#define PACK32_FLUSH(sp,acc,count) if(count < 32) { PACK32_IN_CHECK(sp,acc,count) ; } \
        { uint32_t t32 = (acc << (count-32) ) ; *sp = t32 ; acc = 0 ; count = 64 ; sp++ ; }

// =============================== token extraction ===============================

// initialize bit count and accumulator to extract from a stream of 32 bit UNSIGNED elements
// 32 bits of packed data will be available for extraction, left aligned in accumulator
#define PACK32_OUT_INIT(sp,acc,count) { acc = *sp ; sp++ ; count = 32 ; acc <<= 32 ;}

// UNCHECKED extraction
#define PACK32_OUT_FAST(acc,count,token,nbits) { token = (acc >> (64-nbits)) ; acc <<= (nbits) ; count -= (nbits) ; }

// CHECK to make sure that at least 32 bits are available for extraction
#define PACK32_OUT_CHECK(sp,acc,count) { if(count<32) { acc >>= (32-count) ; acc |= *sp ; acc <<= (32-count) ; sp++ ; count += 32 ;} }

// extract variable length tokens (<= 32 bits) from a 32 bit token stream
// nbits <= 32, acc MUST be uint64_t, sp MUST be a pointer to uint32_t, token may be 32 or 64 bits
#define PACK32_OUT(sp,acc,count,token,nbits) PACK32_OUT_CHECK(sp,acc,count) ; PACK32_OUT_FAST(acc,count,token,nbits)

// extract variable length tokens (> 32 bits) from a 32 bit token stream, token expected to be 64 bits
// extraction is performed in 2 step, the upper nbits -32 bits first, then the lower 32 bits
#define PACK32_OUT_L(sp,acc,count,token,nbits) \
  if((nbits) > 32) { uint32_t t ; PACK32_OUT(sp,acc,count,token,((nbits)-32)) ; token <<= 32 ; PACK32_OUT(sp,acc,count,t,32) ; token |= t ; } \
  else { PACK32_OUT(sp,acc,count,token,nbits) }

#endif
