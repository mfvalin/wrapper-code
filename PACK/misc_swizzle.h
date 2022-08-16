/*
 * Hopefully useful code for C
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
 */
#if defined(IN_FORTRAN_CODE)

interface
! void SplitEvenOdd_8(void *pa, void *pe, void *po)
  subroutine SplitEvenOdd_8(a, e, o) bind(C, name='SplitEvenOdd_8')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void SplitEvenOdd_16(void *pa, void *pe, void *po)
  subroutine SplitEvenOdd_16(a, e, o) bind(C, name='SplitEvenOdd_16')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void SplitEvenOdd_32(void *pa, void *pe, void *po)
  subroutine SplitEvenOdd_32(a, e, o) bind(C, name='SplitEvenOdd_32')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void SplitEvenOdd_64(void *pa, void *pe, void *po)
  subroutine SplitEvenOdd_64(a, e, o) bind(C, name='SplitEvenOdd_64')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void SplitEvenOdd_n(void *pa, void *pe, void *po, int n)
  subroutine SplitEvenOdd_n(a, e, o) bind(C, name='SplitEvenOdd_n')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void ShuffleEvenOdd_8(void *pa, void *pe, void *po)
  subroutine shuffle_eo_8(a, e, o) bind(C, name='ShuffleEvenOdd_8')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void ShuffleEvenOdd_16(void *pa, void *pe, void *po)
  subroutine shuffle_eo_16(a, e, o) bind(C, name='ShuffleEvenOdd_16')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void ShuffleEvenOdd_32(void *pa, void *pe, void *po)
  subroutine shuffle_eo_32(a, e, o) bind(C, name='ShuffleEvenOdd_32')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void ShuffleEvenOdd_64(void *pa, void *pe, void *po)
  subroutine shuffle_eo_64(a, e, o) bind(C, name='ShuffleEvenOdd_64')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine

! void ShuffleEvenOdd_n(void *pa, void *pe, void *po, int n)
  subroutine shuffle_eo_n(a, e, o) bind(C, name='ShuffleEvenOdd_n')
    implicit none
#define IgnoreTypeKindRank a, e, o
#define ExtraAttributes
#include <IgnoreTypeKindRank.hf>
  end subroutine
end interface

#else

// inline C code
// split into even and odd indices
// shuffle even and odd indices back

#include <stdint.h>
#if defined(__x86_64__) && defined(__AVX__) && defined(WITH_SIMD)
#include <immintrin.h>
#endif

static int b_const = 0b01010101 ;

#if ! defined(STATIC)
#define STATIC static
#endif

static int index[8] = {0, 2, 4, 6, 1, 3, 5, 7 } ;

// split a into even and odd terms (8 items)
STATIC inline void SplitEvenOdd_8(void *pa, void *pe, void *po){
  int i, j ;
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 v ;
  __m256i vx ;
  __m128 ve, vo ;
  vx = (__m256i) _mm256_loadu_ps((float *) index) ;
  v = _mm256_loadu_ps(pa) ;
  v = _mm256_permutevar8x32_ps(v, vx) ;
  _mm_storeu_ps(e, _mm256_extractf128_ps(v, 0) ) ;
  _mm_storeu_ps(o, _mm256_extractf128_ps(v, 1) ) ;
// printf("+++\n");
#else
  for(i=0, j=0 ; j<4 ; j++, i+=2){
    e[j] = a[i] ;
    o[j] = a[i+1] ;
  }
#endif
}

// split a into even and odd terms (16 items)
STATIC inline void SplitEvenOdd_16(void *pa, void *pe, void *po){
  int i, j ;
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  SplitEvenOdd_8(a  , e  , o  ) ;
  SplitEvenOdd_8(a+8, e+4, o+4) ;
#else
  for(i=0, j=0 ; j<8 ; j++, i+=2){
    e[j] = a[i] ;
    o[j] = a[i+1] ;
  }
#endif
}

// split a into even and odd terms (32 items)
STATIC inline void SplitEvenOdd_32(void *pa, void *pe, void *po){
  int i, j ;
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  SplitEvenOdd_16(a   , e  , o  ) ;
  SplitEvenOdd_16(a+16, e+8, o+8) ;
#else
  for(i=0, j=0 ; j<16 ; j++, i+=2){
    e[j] = a[i] ;
    o[j] = a[i+1] ;
  }
#endif
}

// split a into even and odd terms (64 items)
STATIC inline void SplitEvenOdd_64(void *pa, void *pe, void *po){
  int i, j ;
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  SplitEvenOdd_32(a   , e   , o   ) ;
  SplitEvenOdd_32(a+32, e+16, o+16) ;
#else
  for(i=0, j=0 ; j<32 ; j++, i+=2){
    e[j] = a[i] ;
    o[j] = a[i+1] ;
  }
#endif
}

// split a into even and odd terms (n items)
// if n is odd, n-1 is used
STATIC inline void SplitEvenOdd_n(void *pa, void *pe, void *po, int n){
  int i, j ;
  int nodd = n >> 1 ;
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;

  while(n >= 64){
    SplitEvenOdd_64(a , e , o ) ;
    a += 64 ; e += 32 ; o += 32 ; n -= 64 ;
  }
  if(n >= 32){
    SplitEvenOdd_32(a , e , o ) ;
    a += 32 ; e += 16 ; o += 16 ; n -= 32 ;
  }
  if(n >= 16){
    SplitEvenOdd_16(a , e , o ) ;
    a += 16 ; e += 8 ; o += 8 ; n -= 16 ;
  }
  if(n >= 8){
    SplitEvenOdd_8(a , e , o ) ;
    a += 8 ; e += 4 ; o += 4 ; n -= 8 ;
  }
  if(n < 8 && n > 0){
    for(i=0, j=0 ; i<(n-1) ; j++, i+=2){
      e[j] = a[i] ;
      o[j] = a[i+1] ;
    }
    if(i<n) e[j] = a[i] ;  // n was odd ?
  }
  return ;
}

// interleave even and odd terms into a (8 items)
STATIC inline void ShuffleEvenOdd_8(void *pa, void *pe, void *po){
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX__) && defined(WITH_SIMD)
  __m128 ve, vo ;
  ve = _mm_loadu_ps(e) ;
  vo = _mm_loadu_ps(o) ;
  _mm_storeu_ps(a  , _mm_unpacklo_ps(ve, vo) ) ;
  _mm_storeu_ps(a+4, _mm_unpackhi_ps(ve, vo) ) ;
#else
  int i, j ;
  for(i=0, j=0 ; j<4 ; j++, i+=2){
    a[i]   = e[j] ;
    a[i+1] = o[j] ;
  }
#endif
}

// interleave even and odd terms into a (16 items)
STATIC inline void ShuffleEvenOdd_16(void *pa, void *pe, void *po){
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX__) && defined(WITH_SIMD)
  __m256 ve, vo, v0, v1 ;
  ve = _mm256_loadu_ps(e) ;
  vo = _mm256_loadu_ps(o) ;
  v0 = _mm256_unpacklo_ps(ve, vo) ;
  v1 = _mm256_unpackhi_ps(ve, vo) ;
  _mm256_storeu_ps(a  , _mm256_permute2f128_ps(v0, v1, 32) ) ;
  _mm256_storeu_ps(a+8, _mm256_permute2f128_ps(v0, v1, 49) ) ;
#else
  int i, j ;
  for(i=0, j=0 ; j<8 ; j++, i+=2){
    a[i]   = e[j] ;
    a[i+1] = o[j] ;
  }
#endif
}

// interleave even and odd terms into a (32 items)
STATIC inline void ShuffleEvenOdd_32(void *pa, void *pe, void *po){
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX__) && defined(WITH_SIMD)
  __m256 ve0, ve1, vo0, vo1, v0, v1, v2, v3 ;
  ve0 = _mm256_loadu_ps(e) ;
  vo0 = _mm256_loadu_ps(o) ;
  ve1 = _mm256_loadu_ps(e+8) ;
  vo1 = _mm256_loadu_ps(o+8) ;
  v0 = _mm256_unpacklo_ps(ve0, vo0) ;
  v1 = _mm256_unpackhi_ps(ve0, vo0) ;
  v2 = _mm256_unpacklo_ps(ve1, vo1) ;
  v3 = _mm256_unpackhi_ps(ve1, vo1) ;
  _mm256_storeu_ps(a  , _mm256_permute2f128_ps(v0, v1, 32) ) ;
  _mm256_storeu_ps(a+8, _mm256_permute2f128_ps(v0, v1, 49) ) ;
  _mm256_storeu_ps(a+16, _mm256_permute2f128_ps(v2, v3, 32) ) ;
  _mm256_storeu_ps(a+24, _mm256_permute2f128_ps(v2, v3, 49) ) ;
#else
  int i, j ;
  for(i=0, j=0 ; j<16 ; j++, i+=2){
    a[i]   = e[j] ;
    a[i+1] = o[j] ;
  }
#endif
}

// interleave even and odd terms into a (64 items)
STATIC inline void ShuffleEvenOdd_64(void *pa, void *pe, void *po){
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX__) && defined(WITH_SIMD)
  ShuffleEvenOdd_32(a   , e   , o   ) ;   // first 32
  ShuffleEvenOdd_32(a+32, e+16, o+16) ;   // next 32
#else
  int i, j ;
  for(i=0, j=0 ; j<32 ; j++, i+=2){
    a[i]   = e[j] ;
    a[i+1] = o[j] ;
  }
#endif
}

// interleave even and odd terms into a (n items) (
// if n is odd, n-1 is used
STATIC inline void ShuffleEvenOdd_n(void *pa, void *pe, void *po, int n){
  int i, j ;
  int nodd = n >> 1;
  float *a = (float *) pa ;
  float *e = (float *) pe ;
  float *o = (float *) po ;
#if defined(__x86_64__) && defined(__AVX__) && defined(WITH_SIMD)
  while(n >= 64){
    ShuffleEvenOdd_64(a   , e   , o   ) ;   // batch of 64
    n -= 64 ; a += 64 ; e += 32 ; o += 32 ;
  }
  if(n >= 32){
    ShuffleEvenOdd_32(a, e, o) ;   // next 32
    n -= 32 ; a += 32 ; e += 16 ; o += 16 ;
  }
  if(n >= 16){
    ShuffleEvenOdd_16(a, e, o) ;   // next 32
    n -= 16 ; a += 16 ; e += 8 ; o += 8 ;
  }
  if(n >= 8){
    ShuffleEvenOdd_8(a, e, o) ;   // next 32
    n -= 8 ; a += 8 ; e += 4 ; o += 4 ;
  }
  for(i=0, j=0 ; j<nodd ; j++, i+=2){
    a[i]   = e[j] ;
    a[i+1] = o[j] ;
  }
#else
  for(i=0, j=0 ; j<nodd ; j++, i+=2){
    a[i]   = e[j] ;
    a[i+1] = o[j] ;
  }
#endif
  if(n&1) a[n-1] = e[nodd] ;   // odd number of values, one extra even value
}

#endif
