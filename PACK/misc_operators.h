/*
Hopefully useful code for C
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
#endif

#define MISC_OPERATORS

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

// number of bits needed to represent a 32 bit signed number
STATIC inline uint32_t BitsNeeded(int32_t what){
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

STATIC inline uint32_t BitsNeeded24(int32_t what){
  int it ;
  union{
    float    f;
    uint32_t u;
  } t ;
  if(what == 0) return 0 ;

  t.f = what ;
  it = ((t.u >> 23) & 0xFF) - 127 + 1 ; // exponent - bias + 1
  it += (t.u >> 31) ;
  return (it > 32) ? 32 : it ;
}

// add to number of bits needed distribution
// what : integer array of dimension n
// pop  ; integer array of dimension 33 containing the number of bits needed distribution
// n    : number of signed integers in what
STATIC inline void BitPop(int32_t *what, uint32_t *pop, int n){
  int i;
  for(i=0 ; i<n ; i++) {
    int nbits = BitsNeeded(what[i]) ;  // number of bits needed for this signed integer
    pop[nbits]++ ;                     // bump count for this number of bits
//     pop[0]++ ;
  }
}

#endif

#else
  interface
    function BitsNeeded(what) result(nbits) bind(C,name='BitsNeeded')
      import C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: what
      integer(C_INT32_T) :: nbits
    end function BitsNeeded
    function BitsNeeded24(what) result(nbits) bind(C,name='BitsNeeded')
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
