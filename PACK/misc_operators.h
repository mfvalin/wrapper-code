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

// C interfaces and declarations

#if ! defined(MISC_OPERATORS)

#include <stdint.h>

#if ! defined(STATIC)
#define STATIC static
#define STATIC_DEFINED_HERE
#endif

#if defined(__x86_64__) && defined(__AVX2__)
#include <immintrin.h>
#endif

#define MISC_OPERATORS

#if defined(__x86_64__) && defined(__AVX2__)
static inline __m128i _mm_lower128(__m256i v256) { return _mm256_extracti128_si256(v256, 0) ; }
static inline __m128i _mm_upper128(__m256i v256) { return _mm256_extracti128_si256(v256, 1) ; }
#endif

// convert a float to a rounded integer
#define FLOAT_TO_INT(x) ( (x)<0 ? (int)((x)-0.5) : (int)((x)+0.5) )
static inline int32_t float_to_int(float x) {
  return FLOAT_TO_INT(x) ;
}

// divide a signed integer by 2 with rounding toward +-infinity
#define IDIV2R(x) (( (x) + 1 + ((x)>>31) ) >> 1 )
#define IDIV2R_256(v) _mm256_srai_epi32(_mm256_add_epi32(_mm256_sub_epi32(v, _mm256_cmpeq_epi32(v, v)), _mm256_srai_epi32(v, 31)), 1)
#define IDIV2R_128(v) _mm_srai_epi32(_mm_add_epi32(_mm_sub_epi32(v, _mm_cmpeq_epi32(v, v)), _mm_srai_epi32(v, 31)), 1)

#if defined(__x86_64__) && defined(__AVX2__)
static __m256i _mm256_idiv2r_epi32(__m256i v){
  return IDIV2R_256(v) ;
}
static __m128i _mm_idiv2r_epi32(__m128i v){
  return IDIV2R_128(v) ;
}
#endif

// divide a signed integer by 4 with rounding toward +-infinity
#define IDIV4R(x) (( (x) + 2 + ((x)>>31) ) >> 2 )
#define IDIV4R_256(v) _mm256_srai_epi32(_mm256_sub_epi32(_mm256_add_epi32(v, _mm256_srai_epi32(v, 31)), _mm256_slli_epi32(_mm256_cmpeq_epi32(v, v), 1)), 2)
#define IDIV4R_128(v) _mm_srai_epi32(_mm_sub_epi32(_mm_add_epi32(v, _mm_srai_epi32(v, 31)), _mm_slli_epi32(_mm_cmpeq_epi32(v, v), 1)), 2)

#if defined(__x86_64__) && defined(__AVX2__)
static inline __m256i _mm256_idiv4r_epi32(__m256i v){
  return IDIV4R_256(v) ;
}
static __m128i _mm_idiv4r_epi32(__m128i v){
  return IDIV4R_128(v) ;
}
#endif

// divide a signed integer by 8 with rounding toward +-infinity
#define IDIV8R(x) (( (x) + 4 + ((x)>>31) ) >> 3 )
#define IDIV8R_256(v) _mm256_srai_epi32(_mm256_sub_epi32(_mm256_add_epi32(v, _mm256_srai_epi32(v, 31)), _mm256_slli_epi32(_mm256_cmpeq_epi32(v, v), 2)), 3)
#define IDIV8R_128(v) _mm_srai_epi32(_mm_sub_epi32(_mm_add_epi32(v, _mm_srai_epi32(v, 31)), _mm_slli_epi32(_mm_cmpeq_epi32(v, v), 2)), 3)

#if defined(__x86_64__) && defined(__AVX2__)
static __m256i _mm256_idiv8r_epi32(__m256i v){
  return IDIV8R_256(v) ;
}
static __m128i _mm_idiv8r_epi32(__m128i v){
  return IDIV8R_128(v) ;
}
#endif

// divide a signed integer by 2 truncating toward zero
#define IDIV2T(x) ((x) + (int32_t) ((uint32_t) (x) >> 31 ) ) >> 1
#define IDIV2T_256(v) _mm256_srai_epi32(_mm256_add_epi32(v, _mm256_srli_epi32(v, 31)), 1)
#define IDIV2T_128(v) _mm_srai_epi32(_mm_add_epi32(v, _mm_srli_epi32(v, 31)), 1)

#if defined(__x86_64__) && defined(__AVX2__)
static __m256i _mm256_idiv2t_epi32(__m256i v){
  return IDIV2T_256(v) ;
}
static __m128i _mm_idiv2t_epi32(__m128i v){
  return IDIV2T_128(v) ;
}
#endif

// divide a signed integer by 4 truncating toward zero
#define IDIV4T(x) ((x) + (int32_t) ((uint32_t) ( (x) >> 1 ) >> 30 ) ) >> 2
#define IDIV4T_256(v) _mm256_srai_epi32(_mm256_add_epi32(v, _mm256_srli_epi32(_mm256_srai_epi32(v, 1), 30)), 2)
#define IDIV4T_128(v) _mm_srai_epi32(_mm_add_epi32(v, _mm_srli_epi32(_mm_srai_epi32(v, 1), 30)), 2)

#if defined(__x86_64__) && defined(__AVX2__)
static __m256i _mm256_idiv4t_epi32(__m256i v){
  return IDIV4T_256(v) ;
}
static __m128i _mm_idiv4t_epi32(__m128i v){
  return IDIV4T_128(v) ;
}
#endif

// divide a signed integer by 8 truncating toward zero
#define IDIV8T(x) ((x) + (int32_t) ((uint32_t) ( (x) >> 1 ) >> 29 ) ) >> 3
#define IDIV8T_256(v) _mm256_srai_epi32(_mm256_add_epi32(v, _mm256_srli_epi32(_mm256_srai_epi32(v, 1), 29)), 3)
#define IDIV8T_128(v) _mm_srai_epi32(_mm_add_epi32(v, _mm_srli_epi32(_mm_srai_epi32(v, 1), 29)), 3)

#if defined(__x86_64__) && defined(__AVX2__)
static __m256i _mm256_idiv8t_epi32(__m256i v){
  return IDIV8T_256(v) ;
}
static __m128i _mm_idiv8t_epi32(__m128i v){
  return IDIV8T_128(v) ;
}
#endif

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
#define LMASK32(nbits)  ((nbits) ? ((~0 )  << (32-nbits)) : 0 )
#define LMASK64(nbits)  ((nbits) ? ((~0l)  << (64-nbits)) : 0 )

// 32 and 64 bit left aligned masks (nbits == 0) NOT supported
#define LMASK32Z(nbits)  ((~0 )  << (32-nbits))
#define LMASK64Z(nbits)  ((~0l)  << (64-nbits))

// 32 and 64 bit right aligned masks
#define RMASK32(nbits)  (((nbits) == 32) ? (~0 ) : (~((~0)  << nbits)))
#define RMASK64(nbits)  (((nbits) == 64) ? (~0l) : (~((~0l) << nbits)))

// 31 and 63 bit right aligned masks (same as above but faster, assuming 32/64 bit case not needed)
#define RMASK31(nbits)  (~((~0)  << nbits))
#define RMASK63(nbits)  (~((~0l) << nbits))

// population count
STATIC inline uint32_t popcnt_32(uint32_t what){
  uint32_t cnt ;
#if defined(__x86_64__)
  // X86 family of processors
  __asm__ __volatile__ ("popcnt{l %1, %0| %0, %1}" : "=r"(cnt) : "r"(what) : "cc" ) ;
#else
  cnt = 0 ;
  while(what & 1){
    cnt++ ;
    what >> = 1 ;
  }
#endif
  return cnt ;
}

STATIC inline uint32_t popcnt_64(uint64_t what){
  uint64_t cnt ;
#if defined(__x86_64__)
  // X86 family of processors
  __asm__ __volatile__ ("popcnt{ %1, %0| %0, %1}" : "=r"(cnt) : "r"(what) : "cc" ) ;
#else
  cnt = 0 ;
  while(what & 1){
    cnt++ ;
    what >> = 1 ;
  }
#endif
  return cnt ;
}

// leading zeros count (32 bit value)
STATIC inline uint32_t lzcnt_32(uint32_t what){
  uint32_t cnt ;
#if defined(__x86_64__)
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
#if defined(__x86_64__)
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

#define TO_ZIGZAG(x)   ( ((x) << 1) ^ ((x) >> 31) )
#define FROM_ZIGZAG(x) ( ((x) >> 1) ^ (-((x) & 1)) )

// convert to sign and magnitude form, sign is Least Significant Bit
STATIC inline uint32_t to_zigzag_32(int32_t what){
  return (what << 1) ^ (what >> 31) ;
//   return (iabs_32(what) << 1) - isign_32(what) ;
}

// convert from sign and magnitude form, sign is Least Significant Bit
STATIC inline int32_t from_zigzag_32(uint32_t what){
  int32_t sign = -(what & 1) ;
  return ((what >> 1) ^ sign) ;
}

STATIC inline int32_t v_to_zigzag_32(int32_t * restrict src, uint32_t * restrict dst, int ni){
  int i ;
  uint32_t max=0 ;
  for(i=0 ; i<ni ; i++){
    dst[i] = to_zigzag_32(src[i]) ;
    max = (dst[i] > max) ? dst[i] : max ;
  }
  return max ;
}

STATIC inline int32_t v_from_zigzag_32(uint32_t * restrict src, int32_t * restrict dst, int ni){
  int32_t i, max=0 ;
  for(i=0 ; i<ni ; i++){
    dst[i] = from_zigzag_32(src[i]) ;
    max = (dst[i] > max) ? dst[i] : max ;
  }
  return max ;
}

#define NBMASK 0xaaaaaaaau /* negabinary<-> 2's complement binary conversion mask */

#define TO_NEGABINARY(x)   ( ((uint32_t)(x) + NBMASK) ^ NBMASK  )
#define FROM_NEGABINARY(x) ( (int32_t)(((x) ^ NBMASK) - NBMASK) )

// signed integer (2's complement) to negabinary (base -2) conversion
STATIC inline uint32_t int_to_negabinary(int32_t x)
{
  return ((uint32_t)x + NBMASK) ^ NBMASK;
}

STATIC inline void v_int_to_negabinary(int32_t x, int32_t n)
{
  int i ;
  for(i=0 ; i<n ; i++) x = int_to_negabinary(x) ;
}

// negabinary (base -2) to signed integer (2's complement) conversion
STATIC inline int32_t negabinary_to_int(uint32_t x)
{
  return (int32_t)((x ^ NBMASK) - NBMASK);
}

STATIC inline void v_negabinary_to_int(uint32_t x, int32_t n)
{
  int i ;
  for(i=0 ; i<n ; i++) x = negabinary_to_int(x) ;
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

// entropy calculator functions
float BitEntropy(int32_t *bitstream, int npts, int nbits, int rshift);
void BitEntropy4(float entropy[4], uint32_t *bitstream, int npts, int nbits, int rshift);

// add to number of bits needed distribution
// what : signed integer array of dimension n
// pop  ; integer array of dimension 33 containing the number of bits needed distribution
// n    : number of signed integers in what
STATIC inline void BitPop(int32_t *what, uint32_t *pop, int n){
  int i;
  for(i=0 ; i<n ; i++) {
    int nbits = BitsNeeded_32(what[i]) ;  // number of bits needed for this signed integer
    pop[nbits]++ ;                     // bump count for this number of bits
  }
}

// add to number of bits needed distribution
// what : unsigned integer array of dimension n
// pop  ; integer array of dimension 33 containing the number of bits needed distribution
// n    : number of signed integers in what
STATIC inline void BitPopU(uint32_t *what, uint32_t *pop, int n){
  int i;
  for(i=0 ; i<n ; i++) {
    int nbits = BitsNeeded_u32(what[i]) ;  // number of bits needed for this unsigned integer
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

!  Fortran interfaces and declarations

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

  interface from_zigzag  ! generic interface
    function from_zigzag_32(what) result(r) bind(C,name='from_zigzag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: r
    end function
    function v_from_zigzag_32_0(src, dest, n) result(r) bind(C,name='v_from_zigzag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN)  :: src
      integer(C_INT32_T), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function v_from_zigzag_32_1(src, dest, n) result(r) bind(C,name='v_from_zigzag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN)  :: src
      integer(C_INT32_T), dimension(*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function v_from_zigzag_32_2(src, dest, n) result(r) bind(C,name='v_from_zigzag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(1,*), intent(IN)  :: src
      integer(C_INT32_T), dimension(1,*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
  end interface

  interface to_negabinary
    function to_negabinary_32(what) result(r) bind(C,name='int_to_negabinary')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: r
    end function
    function to_negabinary_32_0(src, n) result(r) bind(C,name='v_int_to_negabinary')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(INOUT)  :: src
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function to_negabinary_32_1(src, n) result(r) bind(C,name='v_int_to_negabinary')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(INOUT)  :: src
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
  end interface

  interface from_negabinary
    function from_negabinary_32(what) result(r) bind(C,name='negabinary_to_int')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: r
    end function
    function from_negabinary_32_0(src, n) result(r) bind(C,name='v_negabinary_to_int')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(INOUT)  :: src
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function from_negabinary_32_1(src, n) result(r) bind(C,name='v_negabinary_to_int')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(INOUT)  :: src
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
  end interface

  interface to_zigzag  ! generic interface
    function to_zigzag_32(what) result(r) bind(C,name='to_zigzag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: r
    end function
    function v_to_zigzag_32_0(src, dest, n) result(r) bind(C,name='v_to_zigzag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN)  :: src
      integer(C_INT32_T), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function v_to_zigzag_32_1(src, dest, n) result(r) bind(C,name='v_to_zigzag_32')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN)  :: src
      integer(C_INT32_T), dimension(*), intent(OUT) :: dest
      integer(C_INT32_T), intent(IN), value :: n
      integer(C_INT32_T) :: r
    end function
    function v_to_zigzag_32_2(src, dest, n) result(r) bind(C,name='v_to_zigzag_32')
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
      integer(C_INT32_t), dimension(32), intent(INOUT) :: pop
      integer(C_INT32_t), intent(IN), value :: n
    end subroutine BitPop

    subroutine BitPopU(what, pop, n) bind(C, name='BitPopU')
      import :: C_INT32_T
      implicit none
      integer(C_INT32_t), dimension(n), intent(IN) :: what
      integer(C_INT32_t), dimension(32), intent(INOUT) :: pop
      integer(C_INT32_t), intent(IN), value :: n
    end subroutine BitPopU

    function BitEntropy(bitstream, npts, nbits, rshift) result(entropy) bind(C,name='BitEntropy')
      import :: C_INT32_T, C_FLOAT
      implicit none
      integer(C_INT32_t), intent(IN), value :: npts, nbits, rshift
      integer(C_INT32_t), dimension(*) :: bitstream
      real(C_FLOAT) :: entropy
    end function BitEntropy
  end interface
#endif
