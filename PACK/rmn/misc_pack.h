/*
 * Hopefully useful code for C or Fortran
 * Copyright (C) 2022  Recherche en Prevision Numerique
 *
 * This code is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * set of functions to handle insertion/extraction of 1-32 bit tokens into/from a stream of 32 bit unsigned integers
 * should a token be langer thatn 32 bits, it must be split into smaller tokens before insertion/extraction
 */
#if defined(IN_FORTRAN_CODE) || defined(__GFORTRAN__)

type, bind(C) :: PackHeader      ! structure describing quantization information for a set of values
  integer(C_INT) :: offset = 0   ! integer offset reflecting minimum value of packed field
  integer(C_INT) :: exp = 0      ! largest exponent (min, max, range) (with IEEE bias removed)
  integer(C_INT) :: nbits = 0    ! number of useful bits in quantized token
  real(C_FLOAT)  :: quantum      ! quantization interval. quantized value = (value/quantum) - offset
end type

interface
! void float_quantize_prep(int nbits, QuantizeHeader *p, float maxval, float minval, float quantum);
  subroutine float_quantize_prep(nbits, p, maxval, minval, quantum) bind(C,name='float_quantize_prep')
    import :: C_INT, C_FLOAT, PackHeader
    implicit none
    integer(C_INT), intent(IN), value :: nbits
    real(C_FLOAT), intent(IN), value :: maxval, minval, quantum
    type(PackHeader), intent(INOUT) :: p
  end subroutine
! void float_quantize(void *iz, float *z, int ni, int lni, int lniz, int nj, QuantizeHeader *p );
  subroutine float_quantize(iz, z, ni, lni, lniz, nj, p ) bind(C,name='float_quantize')
    import :: C_INT, C_FLOAT, PackHeader
    implicit none
    integer(C_INT), intent(IN), value :: ni, lni, nj, lniz
    integer(C_INT), dimension(*), intent(INOUT) :: iz
    real(C_FLOAT), dimension(lni,nj), intent(IN) :: z
    type(PackHeader), intent(IN) :: p
  end subroutine
! void float_unquantize(void *iz, float *z, int ni, int lni, int lniz, int nj, QuantizeHeader *p);
  subroutine float_unquantize(iz, z, ni, lni, lniz, nj, p) bind(C,name='float_unquantize')
    import :: C_INT, C_FLOAT, PackHeader
    implicit none
    integer(C_INT), intent(IN), value :: ni, lni, nj, lniz
    integer(C_INT), dimension(*), intent(IN) :: iz
    real(C_FLOAT), dimension(lni,nj), intent(INOUT) :: z
    type(PackHeader), intent(IN) :: p
  end subroutine
! float quantum_adjust(float quantum);
  function quantum_adjust(quantum) result(adjusted) bind(C,name='quantum_adjust')
    import :: C_FLOAT
    implicit none
    real(C_FLOAT), intent(IN), value :: quantum
    real(C_FLOAT) :: adjusted
  end function
! uint32_t float_quantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, IntPair *t);
  function float_quantize_simple(z, q, ni, lniz, lniq, nj, quantum, t) result(nbits) bind(C,name='float_quantize_simple')
    import :: C_INT, C_FLOAT
    implicit none
    real(C_FLOAT), dimension(*), intent(IN) :: z
    integer(C_INT), dimension(*), intent(OUT) :: q
    integer(C_INT), intent(IN), value :: ni, lniq, nj, lniz
    real(C_FLOAT), intent(IN), value :: quantum
    integer(C_INT), dimension(2), intent(OUT) :: t
    integer(C_INT) :: nbits
  end function
! void float_unquantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, FloatPair *t);
  subroutine float_unquantize_simple(z, q, ni, lniz, lniq, nj, quantum, t) bind(C,name='float_unquantize_simple')
    import :: C_INT, C_FLOAT
    implicit none
    real(C_FLOAT), dimension(*), intent(OUT) :: z
    integer(C_INT), dimension(*), intent(IN) :: q
    integer(C_INT), intent(IN), value :: ni, lniq, nj, lniz
    real(C_FLOAT), intent(IN), value :: quantum
    real(C_FLOAT), dimension(2), intent(OUT) :: t
  end subroutine
end interface

! int float_info(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs, float *spval, uint32_t spmask);
interface float_info  ! generic interface for both missing and no missing cases
! int float_info_no_missing(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs);
  function float_info_no_missing(zz, ni, lni, nj, maxval, minval, minabs) result(n) bind(C,name='float_info_no_missing')
    import :: C_INT, C_FLOAT
    implicit none
#define IgnoreTypeKindRank zz
#define ExtraAttributes 
#include <rmn/IgnoreTypeKindRank.hf>
    integer(C_INT), intent(IN), value :: ni, lni, nj
!   real(C_FLOAT), dimension(lni,nj)), intent(IN) :: zz
    real(C_FLOAT), intent(OUT) :: maxval, minval, minabs
    integer(C_INT) :: n
  end function
! int float_info_missing(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs, float *spval, uint32_t spmask);
  function float_info_missing(zz, ni, lni, nj, maxval, minval, minabs, spval, spmask) result(n) bind(C,name='float_info_missing')
    import :: C_INT, C_FLOAT
    implicit none
#define IgnoreTypeKindRank zz
#define ExtraAttributes 
#include <rmn/IgnoreTypeKindRank.hf>
    integer(C_INT), intent(IN), value :: ni, lni, nj, spmask
!   real(C_FLOAT), dimension(lni,nj), intent(IN) :: zz
    real(C_FLOAT), intent(OUT) :: maxval, minval, minabs, spval
    integer(C_INT) :: n
  end function
end interface

#else

#if ! defined(MISC_PACK)
#define MISC_PACK

#include <stdint.h>
#include <stddef.h>

#include <rmn/misc_types.h>

// determine how many bits are needed to represent value x
static inline uint32_t NeedBits(int32_t nx){
  uint32_t n ;
  uint32_t x ;
  uint32_t xtra = (nx < 0) ? 1 : 0 ;  // one extra bit is necessary if nx is negative
  x = (nx < 0) ? -nx : nx ;           // absolute value of nx
#if defined(__x86_64__)
  __asm__ __volatile__ ("lzcnt{l %1, %0| %0, %1}" : "=r"(n) : "r"(x) : "cc");
  n -= xtra ;
  return 32 - n ;
#elif defined(__aarch64__)
  __asm__ __volatile__ ("clz %w[out], %w[in]" : [out]"=r"(n) : [in]"r"(x) );
  n -= xtra ;
  return 32 - n ;
#else
  uint32_t m = 0xFFFF ;
  n = 0 ;
  n  += ((x > m) ? 16 : 0) ;   // something in the upper 16 bits
  x >>= ((x > m) ? 16 : 0) ;   // lower 16 bits
  m >>= 8 ;
  n  += ((x > m) ? 8 : 0) ;    // something in bits 8-15
  x >>= ((x > m) ? 8 : 0) ;    // lower 8 bits
  m >>= 4 ;
  n  += ((x > m) ? 4 : 0) ;    // something in bits 4-7
  x >>= ((x > m) ? 4 : 0) ;    // lower 4 bits
  m >>= 2 ;
  n  += ((x > m) ? 2 : 0) ;    // something in bits 2-3
  x >>= ((x > m) ? 2 : 0) ;    // lower 2 bits
  m >>= 1 ;
  n  += ((x > m) ? 1 : 0) ;    // something in bit 1
  x >>= ((x > m) ? 1 : 0) ;    // lower 1 bit
  m >>= 1 ;
  n  += ((x > m) ? 1 : 0) ;    // something in bit 0
  return n + xtra ;
#endif
}

// linear quantization header
typedef struct {
  int o ;      // integer offset reflecting minimum value of packed field
  int e ;      // largest exponent (min, max, range) (with IEEE bias removed)
  int nbits ;  // number of useful bits in quantized token
  float q ;    // quantization interval. quantized value = (value/quantum) - offset
} QuantizeHeader;

// linear quantization functions
void float_quantize_prep(int nbits, QuantizeHeader *p, float maxval, float minval, float quantum);
void float_quantize(void *iz, float *z, int ni, int lni, int lniz, int nj, QuantizeHeader *p );
uint32_t float_quantize_simple_1D(float * restrict z, int32_t * restrict q, int ni, float quantum, IntPair *t);
uint32_t float_quantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, IntPair *t);

void float_unquantize(void *iz, float *z, int ni, int lni, int lniz, int nj, QuantizeHeader *p);
void float_unquantize_simple_1D(float * restrict z, int32_t * restrict q, int ni, float quantum, FloatPair *t);
void float_unquantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, FloatPair *t);

float quantum_adjust(float quantum);
int float_info_no_missing(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs);
int float_info_missing(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs, float *spval, uint32_t spmask);
int float_info(float *zz, int ni, int lni, int nj, float *maxval, float *minval, float *minabs, float *spval, uint32_t spmask);

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

#endif
