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
#if defined(__x86_64__x)
  // X86 family of processors
  __asm__ __volatile__ ("lzcnt{l %1, %0| %0, %1}" : "=r"(cnt) : "r"(what) : "cc" ) ;
#elif defined(__aarch64__)
  // ARM family of processors
   __asm__ __volatile__ ("clz %w[out], %w[in]" : [out]"=r"(cnt) : [in]"r"(what) ) ;
#else
   // generic code
   cnt = 32;
   if(what >> 16) { cnt-= 16 ; what >>= 16 ; } ; // bits (16-31) not 0
   if(what >>  8) { cnt-=  8 ; what >>=  8 ; } ; // bits ( 8-15) not 0
   if(what >>  4) { cnt-=  4 ; what >>=  4 ; } ; // bits ( 4- 7) not 0
   if(what >>  2) { cnt-=  2 ; what >>=  2 ; } ; // bits ( 2- 3) not 0
   if(what >>  1) { cnt-=  1 ; what >>=  1 ; } ; // bits ( 1- 1) not 0
   if(what) cnt-- ;                            // bit 0 not 0 ;
#endif
  return cnt ;
}

// leading ones count (32 bit value)
STATIC inline uint32_t lnzcnt_32(uint32_t what){
  return lzcnt_32(~what) ;
}

// leading zeros count (64 bit value)
STATIC inline uint32_t lzcnt_64(uint64_t what){
  uint64_t cnt ;
#if defined(__x86_64__x)
  // X86 family of processors
  __asm__ __volatile__ ("lzcnt{ %1, %0| %0, %1}" : "=r"(cnt) : "r"(what) : "cc" ) ;
#elif defined(__aarch64__)
  // ARM family of processors
   __asm__ __volatile__ ("clz %[out], %[in]" : [out]"=r"(cnt) : [in]"r"(what) ) ;
#else
   // generic code
   cnt = 64;
   if(what >> 32) { cnt-= 32 ; what >>= 32 ; } ; // bits (32-63) not 0
   if(what >> 16) { cnt-= 16 ; what >>= 16 ; } ; // bits (16-31) not 0
   if(what >>  8) { cnt-=  8 ; what >>=  8 ; } ; // bits ( 8-15) not 0
   if(what >>  4) { cnt-=  4 ; what >>=  4 ; } ; // bits ( 4- 7) not 0
   if(what >>  2) { cnt-=  2 ; what >>=  2 ; } ; // bits ( 2- 3) not 0
   if(what >>  1) { cnt-=  1 ; what >>=  1 ; } ; // bits ( 1- 1) not 0
   if(what) cnt-- ;                            // bit 0 not 0 ;
#endif
  return cnt ;
}

// leading ones count (64 bit value)
STATIC inline uint32_t lnzcnt_64(uint64_t what){
  return lzcnt_64(~what) ;
}

// absolute value of signed 32 bit integer (no if)
STATIC inline uint32_t iabs_32(int32_t what){
  int32_t sign = (what >> 31) ;
  return (what ^ sign) - sign ;
}

// sign of signed 32 bit integer (return 0 or -1)
STATIC inline int32_t isign_32(int32_t what){
  return (what >> 31) ;
}

// convert to sign and magnitude form, sign is Least Significant Bit
STATIC inline uint32_t to_isignmag_32(int32_t what){
  return (what << 1) ^ (what >> 31) ;
//   return (iabs_32(what) << 1) - isign_32(what) ;
}

// convert from sign and magnitude form, sign is Least Significant Bit
STATIC inline int32_t from_isignmag_32(uint32_t what){
  int32_t sign = -(what & 1) ;
  return ((what >> 1) ^ sign) ;
}

STATIC inline int32_t v_to_isignmag_32(int32_t * restrict src, uint32_t * restrict dst, int ni){
  int i ;
  uint32_t max=0 ;
  for(i=0 ; i<ni ; i++){
    dst[i] = to_isignmag_32(src[i]) ;
    max = (dst[i] > max) ? dst[i] : max ;
  }
  return max ;
}

STATIC inline int32_t v_from_isignmag_32(uint32_t * restrict src, int32_t * restrict dst, int ni){
  int32_t i, max=0 ;
  for(i=0 ; i<ni ; i++){
    dst[i] = from_isignmag_32(src[i]) ;
    max = (dst[i] > max) ? dst[i] : max ;
  }
  return max ;
}

// number of bits needed to represent a 32 bit unsigned number
// uses lzcnt_32 function, that uses the lzcnt instruction
STATIC inline uint32_t BitsNeeded_u32(uint32_t what){
  return 32 - lzcnt_32(what) ;
}

// number of bits needed to represent a 32 bit signed number
// uses lzcnt_32 function, that uses the lzcnt instruction
STATIC inline uint32_t BitsNeeded_32(int32_t what){
  union {
    int32_t  i ;
    uint32_t u ;
  }iu ;
  uint32_t nbits = 33 - lzcnt_32(what) ; // there must be a 0 bit at the front
  if(what >= 0) goto end ;
  iu.i = what ;                          // what < 0
  nbits = 33 - lzcnt_32(~iu.u) ;         // one's complement, then count leading zeros
end:
  return (nbits > 32) ? 32 : nbits ;     // max is 32 bits
}

// vector versions of above
int32_t vBitsNeeded_32(int32_t * restrict what, int32_t * restrict bits, int n);
int32_t vBitsNeeded_u32(uint32_t * restrict what, int32_t * restrict bits, int n);

// number of bits needed to represent a 64 bit unsigned number
// uses lzcnt_64 function, that uses the lzcnt instruction
STATIC inline uint32_t BitsNeeded_u64(uint64_t what){
  return 64 - lzcnt_64(what) ;
}

// number of bits needed to represent a 64 bit signed number
// uses lzcnt_64 function, that uses the lzcnt instruction
STATIC inline uint32_t BitsNeeded_64(int64_t what){
  union {
    int64_t  i ;
    uint64_t u ;
  }iu ;
  uint32_t nbits = 65 - lzcnt_64(what) ; // there must be a 0 bit at the front
  if(what >= 0) goto end ;
  iu.i = what ;                          // what < 0
  nbits = 65 - lzcnt_64(~iu.u) ;         // one's complement, then count leading zeros
end:
  return (nbits > 64) ? 64 : nbits ;     // max is 64 bits
}

// vector versions of above
int32_t vBitsNeeded_64(int64_t * restrict what, int32_t * restrict bits, int n);
int32_t vBitsNeeded_u64(uint64_t * restrict what, int32_t * restrict bits, int n);

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
  interface BitsNeeded_u  ! generic interface
    function BitsNeeded_u32(what) result(nbits) bind(C,name='BitsNeeded_u32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function BitsNeeded_u32
    function vBitsNeeded_u32_0(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_u32') ! rank 0
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN)  :: what
      integer(C_INT32_T), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_u32_0
    function vBitsNeeded_u32_1(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_u32') ! rank 1
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN)  :: what
      integer(C_INT32_T), dimension(*), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_u32_1
  end interface

  interface BitsNeeded  ! generic interface
    function BitsNeeded_32(what) result(nbits) bind(C,name='BitsNeeded_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function BitsNeeded_32
    function vBitsNeeded_32_0(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_32') ! rank 0
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN)  :: what
      integer(C_INT32_T), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_32_0
    function vBitsNeeded_32_1(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_32') ! rank 1
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN)  :: what
      integer(C_INT32_T), dimension(*), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_32_1
    function vBitsNeeded_32_2(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_32') ! rank 2
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(1,*), intent(IN)  :: what
      integer(C_INT32_T), dimension(1,*), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_32_2
    function vBitsNeeded_32_3(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_32') ! rank 3
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
    function vBitsNeeded_64_1(what, bits, n) result(nbits) bind(C,name='vBitsNeeded_64') ! rank 1
      import C_INT32_T, C_INT64_T
      implicit none
      integer(C_INT64_T), dimension(*), intent(IN)  :: what
      integer(C_INT32_T), dimension(*), intent(OUT) :: bits
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: nbits
    end function vBitsNeeded_64_1
  end interface

  interface lzcnt  ! generic interface
    function lzcnt_32(what) result(nbits) bind(C,name='lzcnt_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function lzcnt_32
    function lzcnt_64(what) result(nbits) bind(C,name='lzcnt_64')
      import C_INT32_T, C_INT64_T
      implicit none
      integer(C_INT64_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function lzcnt_64
  end interface

  interface lnzcnt  ! generic interface
    function lnzcnt_32(what) result(nbits) bind(C,name='lnzcnt_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function lnzcnt_32
    function lnzcnt_64(what) result(nbits) bind(C,name='lnzcnt_64')
      import C_INT32_T, C_INT64_T
      implicit none
      integer(C_INT64_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function lnzcnt_64
  end interface

  interface from_isignmag  ! generic interface
    function from_isignmag_32(what) result(r) bind(C,name='from_isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: r
    end function
    function v_from_isignmag_32_0(src, dest, n) result(r) bind(C,name='v_from_isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN)  :: src
      integer(C_INT32_T), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function v_from_isignmag_32_1(src, dest, n) result(r) bind(C,name='v_from_isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN)  :: src
      integer(C_INT32_T), dimension(*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function v_from_isignmag_32_2(src, dest, n) result(r) bind(C,name='v_from_isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(1,*), intent(IN)  :: src
      integer(C_INT32_T), dimension(1,*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
  end interface

  interface to_isignmag  ! generic interface
    function to_isignmag_32(what) result(r) bind(C,name='to_isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: r
    end function
    function v_to_isignmag_32_0(src, dest, n) result(r) bind(C,name='v_to_isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN)  :: src
      integer(C_INT32_T), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function v_to_isignmag_32_1(src, dest, n) result(r) bind(C,name='v_to_isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN)  :: src
      integer(C_INT32_T), dimension(*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function v_to_isignmag_32_2(src, dest, n) result(r) bind(C,name='v_to_isignmag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(1,*), intent(IN)  :: src
      integer(C_INT32_T), dimension(1,*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
  end interface

  interface
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
