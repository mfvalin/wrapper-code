/*
Hopefully useful code for C and Fortran
Copyright (C) 2022  Recherche en Prevision Numerique

This code is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation,
version 2.1 of the License.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.
*/
#if ! defined(IN_FORTRAN_CODE) && ! defined(__GFORTRAN__)

#if ! defined(MISC_OPERATORS)

#include <stdint.h>

#if ! defined(STATIC)
#define STATIC static
#define STATIC_DEFINED_HERE
#endif

#define MISC_OPERATORS

// divide a signed integer by 2 with rounding toward +-infinity
#define IDIV2R(x) (( (x) + 1 + ((x)>>31) ) >> 1 )

// divide a signed integer by 4 with rounding toward +-infinity
#define IDIV4R(x) (( (x) + 2 + ((x)>>31) ) >> 2 )

// divide a signed integer by 8 with rounding toward +-infinity
#define IDIV8R(x) (( (x) + 4 + ((x)>>31) ) >> 3 )

#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define ABS(val)  ((val) < 0) ? (-(val)) : (val)

#define MINMAX(min,max,src,n)  \
{ int i=1 ; min = max = src[0] ; \
  while(i++ < n) { min = MIN(min,src[i]) ; max = MAX(max,src[i]); } \
}

// number of bits needed to represent range (positive number)
#define NEEDBITS(range,needed) { uint64_t rng = (range) ; needed = 1; while (rng >>= 1) needed++ ; }

// 32 and 64 bit left aligned masks
#define LMASK32(nbits)  ((~0)   << (32-nbits))
#define LMASK64(nbits)  ((~0l)  << (64-nbits))

// 32 and 64 bit right aligned masks
#define RMASK32(nbits)  (~((~0)  << nbits))
#define RMASK64(nbits)  (~((~0l) << nbits))

// leading zeros count (32 bit value)
STATIC inline uint32_t lzcnt_32(uint32_t what){
  uint32_t cnt ;
#if defined(__x86_64__)
  __asm__ __volatile__ ("lzcnt{l %1, %0| %0, %1}" : "=r"(cnt) : "r"(what) : "cc" ) ;
#elif defined(__aarch64__)
   __asm__ __volatile__ ("clz %w[out], %w[in]" : [out]"=r"(cnt) : [in]"r"(what) ) ;
#endif
  return cnt ;
}

// leading zeros count (64 bit value)
STATIC inline uint32_t lzcnt_64(uint64_t what){
  uint64_t cnt ;
#if defined(__x86_64__)
  __asm__ __volatile__ ("lzcnt{ %1, %0| %0, %1}" : "=r"(cnt) : "r"(what) : "cc" ) ;
#elif defined(__aarch64__)
   __asm__ __volatile__ ("clz %[out], %[in]" : [out]"=r"(cnt) : [in]"r"(what) ) ;
#endif
  return cnt ;
}

STATIC inline uint32_t iabs_32(int32_t what){
  int32_t sign = (what >> 31) ;
  return (what ^ sign) - sign ;
}

STATIC inline int32_t isign_32(int32_t what){
  return (what >> 31) ;
}

// convert to sign and magnitude form, sign is Least Significant Bit
STATIC inline uint32_t isignmag_32(int32_t what){
  return (iabs_32(what) << 1) - isign_32(what) ;
}
STATIC inline int32_t visignmag_32(int32_t * restrict src, int32_t * restrict dst, int ni){
  int i, max=0 ;
  for(i=0 ; i<ni ; i++){
    dst[i] = isignmag_32(src[i]) ;
    max = (dst[i] > max) ? dst[i] : max ;
  }
  return max ;
}
STATIC inline int32_t visignmag_32_inplace(int32_t * restrict src, int ni){
  int i, max=0 ;
  for(i=0 ; i<ni ; i++){
    src[i] = isignmag_32(src[i]) ;
    max = (src[i] > max) ? src[i] : max ;
  }
  return max ;
}

// number of bits needed to represent a 32 bit signed number
// uses lzcnt_32 function, that uses the lzcnt instruction
STATIC inline uint32_t BitsNeeded_32(int32_t what){
  union {
    int32_t  i ;
    uint32_t u ;
  }iu ;
  uint32_t nbits ;
  if(what >= 0) return 32 - lzcnt_32(what) ;
  iu.i = what ; // - 1 ;          // what < 0
  nbits = 33 - lzcnt_32(~iu.u) ;
  return (nbits > 32) ? 32 : nbits ;
}
int32_t vBitsNeeded_32(int32_t * restrict what, int32_t * restrict bits, int n);

// number of bits needed to represent a 64 bit signed number
// uses lzcnt_64 function, that uses the lzcnt instruction
STATIC inline uint32_t BitsNeeded_64(int64_t what){
  union {
    int64_t  i ;
    uint64_t u ;
  }iu ;
  uint32_t nbits ;
  if(what >= 0) return 64 - lzcnt_64(what) ;
  iu.i = what - 1 ;          // what < 0
  nbits = 65 - lzcnt_64(~iu.u) ;
  return (nbits > 64) ? 64 : nbits ;
}
int32_t vBitsNeeded_64(int64_t * restrict what, int32_t * restrict bits, int n);

// number of bits needed to represent a 32 bit signed number
// sleigh of hand using the IEEE double exponent to determine number of bits
STATIC inline uint32_t BitsNeeded32(int32_t what){
  int it ;
  union{
    double   f;
    uint64_t u;
  } t ;
  if(what == 0) return 0 ;

  t.f = what ;
  it = ((t.u >> 52) & 0x7FF) - 1023 + 1 ; // exponent - bias + 1
  it += (t.u >> 63) ;
  return (it > 32) ? 32 : it ;
}

// number of bits needed to represent a 24 bit signed number (inaccurate above 24 bits)
// sleigh of hand using the IEEE float exponent to determine number of bits
STATIC inline uint32_t BitsNeeded24(int32_t what){
  int it ;
  union{
    float    f;
    uint32_t u;
  } t ;
  if(what == 0) return 0 ;

  t.f = what ;
  it = ((t.u >> 23) & 0xFF) - 127 + 1 ; // exponent - bias + 1
  it += (t.u >> 31) ;                   // add 1 if number is negative
  return (it > 32) ? 32 : it ;
}

// add to number of bits needed distribution
// what : integer array of dimension n
// pop  ; integer array of dimension 33 containing the number of bits needed distribution
// n    : number of signed integers in what
STATIC inline void BitPop(int32_t *what, uint32_t *pop, int n){
  int i;
  for(i=0 ; i<n ; i++) {
    int nbits = BitsNeeded_32(what[i]) ;  // number of bits needed for this signed integer
    pop[nbits]++ ;                     // bump count for this number of bits
  }
}

// nearest integer, .5 goes toward +infinity, -.5 goes toward -infinity
// should be equivalent to Fortran NINT intrinsic
STATIC inline int Nint(float what){
  union{
    float f;
    int i;
  } a, b ;
  int i ;
  a.f = .5 ;
  b.f = what ;
  a.i = a.i | (b.i & 0x80000000) ;  // transfer sign of what to .5
  i = a.f + b.f ;                   // C conversion
  return i ;
}

#if defined(STATIC_DEFINED_HERE)
#undef STATIC
#undef STATIC_DEFINED_HERE
#endif

#endif

#else
  interface BitsNeeded  ! generic interface
    function BitsNeeded_32(what) result(nbits) bind(C,name='BitsNeeded_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function BitsNeeded_32
    function vBitsNeeded_32_0(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN)  :: what
      integer(C_INT32_T), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_32_0
    function vBitsNeeded_32_1(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN)  :: what
      integer(C_INT32_T), dimension(*), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_32_1
    function vBitsNeeded_32_2(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(1,*), intent(IN)  :: what
      integer(C_INT32_T), dimension(1,*), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_32_2
    function vBitsNeeded_32_3(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(1,1,*), intent(IN)  :: what
      integer(C_INT32_T), dimension(1,1,*), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_32_3
    function BitsNeeded_64(what) result(nbits) bind(C,name='BitsNeeded_64')
      import C_INT32_T, C_INT64_T
      implicit none
      integer(C_INT64_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function BitsNeeded_64
    function vBitsNeeded_64(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_64')
      import C_INT32_T, C_INT64_T
      implicit none
      integer(C_INT64_T), dimension(*), intent(IN)  :: what
      integer(C_INT32_T), dimension(*), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_64
  end interface
  interface lzcnt
    function lzcnt_32(what) result(nbits) bind(C,name='lzcnt_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function lzcnt_32
    function lzcnt_64(what) result(nbits) bind(C,name='lzcnt_32')
      import C_INT32_T, C_INT64_T
      implicit none
      integer(C_INT64_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function lzcnt_64
  end interface
  interface isignmag
    function isignmag_32(what) result(r) bind(C,name='isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: r
    end function isignmag_32
    function visignmag_32_0(src, dest, n) result(r) bind(C,name='visignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN)  :: src
      integer(C_INT32_T), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function visignmag_32_0
    function visignmag_32_1(src, dest, n) result(r) bind(C,name='visignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN)  :: src
      integer(C_INT32_T), dimension(*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function visignmag_32_1
    function visignmag_32_2(src, dest, n) result(r) bind(C,name='visignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(1,*), intent(IN)  :: src
      integer(C_INT32_T), dimension(1,*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function visignmag_32_2
  end interface
  interface
    function BitsNeeded24(what) result(nbits) bind(C,name='BitsNeeded24')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function BitsNeeded24
    subroutine BitPop(what, pop, n) bind(C, name='BitPop')
      import :: C_INT32_T
      implicit none
      integer(C_INT32_t), dimension(n), intent(IN) :: what
      integer(C_INT32_t), dimension(34), intent(INOUT) :: pop
      integer(C_INT32_t), intent(IN), value :: n
    end subroutine BitPop
    function BitEntropy(bitstream, npts, nbits, rshift) result(entropy) bind(C,name='BitEntropy')
      import :: C_INT32_T, C_FLOAT
      implicit none
      integer(C_INT32_t), intent(IN), value :: npts, nbits, rshift
      integer(C_INT32_t), dimension(*) :: bitstream
      real(C_FLOAT) :: entropy
    end function BitEntropy
  end interface
#endif
