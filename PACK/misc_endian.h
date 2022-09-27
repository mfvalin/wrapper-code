/* Library of useful routines for C and FORTRAN programming
 * Copyright (C) 2019  Environnement Canada
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

// functions to perform endianness swaps

#if ! defined(FAST_ENDIAN)
#define FAST_ENDIAN

#if defined(PROTOTYPES_ONLY)

void Swap_8_in_16(void *s, void *d, int n);   // bytes in halfwords
void Swap_8_in_32(void *s, void *d, int n);   // bytes in words
void Swap_16_in_32(void *s, void *d, int n);  // halfwords in words
void Swap_8_in_64(void *s, void *d, int n);   // bytes in doublewords
void Swap_16_in_64(void *s, void *d, int n);  // halfwords in doublewords
void Swap_32_in_64(void *s, void *d, int n);  // words in doublewords

#else

#include <stdint.h>
#if !defined(NO_SIMD)
  #if defined(__x86_64__)
    #include <immintrin.h>
    #if !defined(WITH_SIMD)
      #define WITH_SIMD
    #endif
  #endif
#endif

// STATIC may be defined as extern, to generate real entry points
#if ! defined(STATIC)
#define STATIC static
#endif

#if defined(__AVX2__) && defined(WITH_SIMD)
// tables for SIMD endian swaps (byte per byte shuffle table) (16 bytes = 128 bits)
// endian swap of bytes in 16 bit tokens (s -> d), index list for 256 bit SIMD shuffle operation
static uint8_t swapindex_8_16[]  = {  1,  0,  3,  2,  5,  4,  7,  6,  9,  8, 11, 10, 13, 12, 15, 14};
// endian swap of bytes in 32 bit tokens (s -> d), index list for 256 bit SIMD shuffle operation
static uint8_t swapindex_8_32[]  = {  3,  2,  1,  0,  7,  6,  5,  4, 11, 10,  9,  8, 15, 14, 13, 12};
// endian swap of 16 bit halfwords in 32 bit tokens (s -> d), index list for 256 bit SIMD shuffle operation
static uint8_t swapindex_16_32[] = {  2,  3,  0,  1,  6,  7,  4,  5, 10, 11,  8,  9, 14, 15, 12, 13};
// endian swap of bytes in 64 bit tokens (s -> d), index list for 256 bit SIMD shuffle operation
static uint8_t swapindex_8_64[]  = {  7,  6,  5,  4,  3,  2,  1,  0, 15, 14, 13, 12, 11, 10,  9,  8};
// endian swap of 16 bit halfwords in 64 bit tokens (s -> d), index list for 256 bit SIMD shuffle operation
static uint8_t swapindex_16_64[] = {  6,  7,  4,  5,  2,  3,  0,  1, 14, 15, 12, 13, 10, 11,  8,  9};
// endian swap of 32 bit words in 64 bit tokens (s -> d), index list for 256 bit SIMD shuffle operation
#endif

// endian swap of 8 bit tokens in 16 bit tokens (s -> d)
// s [IN]  : pointer to source of 16 bit elements (halfwords)
// d [OUT] : pointer to destination 16 bit elements (halfwords)
// n [IN]  : number of 16 bit elements (halfwords) to swap
// s and d may be the same address (swap INPLACE)
STATIC void inline Swap_8_in_16(void *s, void *d, int n){
#if defined(__AVX2__) && defined(WITH_SIMD)
  __m256i ix, vs0, vs1;
  __m128i i0, v1;
  int n2;
  uint16_t *s0 = (uint16_t *) s;
  uint16_t *d0 = (uint16_t *) d;
#endif
  int i;
  uint16_t t;
  uint16_t *s1 = (uint16_t *) s;
  uint16_t *d1 = (uint16_t *) d;

  i = 0;
#if defined(__AVX2__) && defined(WITH_SIMD)
  s0 = s1;
  d0 = d1;
  n2 = (n >> 5);        // number of 32 token chunks
  s1 = s0 + (n2 << 4);  // "halfway" point
  d1 = d0 + (n2 << 4);  // "halfway" point
  ix = _mm256_setzero_si256();                         // unnecessary code to avoid compiler warning
  i0 = _mm_loadu_si128((__m128i const *) swapindex_8_16);
  ix = _mm256_inserti128_si256(ix, i0, 0);             // same shuffle pattern on low and high part
  ix = _mm256_inserti128_si256(ix, i0, 1);
  for( ; i < n - 31 ; i += 32){                        // endian swap of 16 32 bit tokens
    vs0 = _mm256_loadu_si256((__m256i const *) s0);    // load 8 32 bit tokens
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs0 = _mm256_shuffle_epi8 (vs0, ix);               // shuffle (vector bswap equivalent)
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d0, vs0);          // store 8 32 bit tokens
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s0 += 16;
    s1 += 16;
    d0 += 16;
    d1 += 16;
  }
  for( ; i < n - 15 ; i += 16){
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s1 += 16;
    d1 += 16;
  }
  for( ; i < n - 7 ; i += 8){
    v1 = _mm_loadu_si128((__m128i const *) s1);    // load 4 32 bit tokens
    v1 = _mm_shuffle_epi8 (v1, i0);                // shuffle (vector bswap equivalent)
    _mm_storeu_si128((__m128i *) d1, v1);          // store 4 32 bit tokens
    s1 += 8;
    d1 += 8;
  }
  s1 = (uint16_t *) s;
  d1 = (uint16_t *) d;
#endif
  for( ; i<n ; i++){ t = s1[i] ; d1[i] = (t >> 8) | (t<<8); }
}

// endian swap of 8 bit tokens in 32 bit tokens (s -> d)
// s [IN]  : pointer to source of 32 bit elements (words)
// d [OUT] : pointer to destination 32 bit elements (words)
// n [IN]  : number of 32 bit elements (words) to swap
// s and d may be the same address (swap INPLACE)
STATIC void inline Swap_8_in_32(void *s, void *d, int n){
#if defined(__AVX2__) && defined(WITH_SIMD)
  __m256i ix, vs0, vs1;
  __m128i i0, v1;
  int n2;
  uint32_t *s0, *d0;
#endif
  int i;
  uint32_t t;
  uint32_t *s1 = (uint32_t *) s;
  uint32_t *d1 = (uint32_t *) d;

  i = 0;
#if defined(__AVX2__) && defined(WITH_SIMD)
  s0 = s1;
  d0 = d1;
  n2 = (n >> 4);        // number of 16 token chunks
  s1 = s0 + (n2 << 3);  // "halfway" point
  d1 = d0 + (n2 << 3);  // "halfway" point
  ix = _mm256_setzero_si256();                         // unnecessary code to avoid compiler warning
  i0 = _mm_loadu_si128((__m128i const *) swapindex_8_32);
  ix = _mm256_inserti128_si256(ix, i0, 0);             // same shuffle pattern on low and high part
  ix = _mm256_inserti128_si256(ix, i0, 1);
  for( ; i < n - 15 ; i += 16){                        // endian swap of 16 32 bit tokens
    vs0 = _mm256_loadu_si256((__m256i const *) s0);    // load 8 32 bit tokens
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs0 = _mm256_shuffle_epi8 (vs0, ix);               // shuffle (vector bswap equivalent)
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d0, vs0);          // store 8 32 bit tokens
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s0 += 8;
    s1 += 8;
    d0 += 8;
    d1 += 8;
  }
  for( ; i < n - 7 ; i += 8){
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s1 += 8;
    d1 += 8;
  }
  for( ; i < n - 3 ; i += 4){
    v1 = _mm_loadu_si128((__m128i const *) s1);    // load 4 32 bit tokens
    v1 = _mm_shuffle_epi8 (v1, i0);                // shuffle (vector bswap equivalent)
    _mm_storeu_si128((__m128i *) d1, v1);          // store 4 32 bit tokens
    s1 += 4;
    d1 += 4;
  }
#endif
  for( ; i < n ; i++){    // loop over remainder
    t = *s1++;
    t = (t >> 24) | (t << 24) | ((t >> 8) & 0xFF00) | ((t << 8) & 0xFF0000);    // bswap
    *d1++ = t;
  }
}

// endian swap of 16 bit tokens in 32 bit tokens (s -> d)
// s [IN]  : pointer to source of 32 bit elements (words)
// d [OUT] : pointer to destination 32 bit elements (words)
// n [IN]  : number of 32 bit elements (words) to swap
// s and d may be the same address (swap INPLACE)
STATIC void inline Swap_16_in_32(void *s, void *d, int n){
#if defined(__AVX2__) && defined(WITH_SIMD)
  __m256i ix, vs0, vs1;
  __m128i i0, v1;
  int n2;
  uint32_t *s0, *d0;
#endif
  int i;
  uint32_t t;
  uint32_t *s1 = (uint32_t *) s;
  uint32_t *d1 = (uint32_t *) d;

  i = 0;
#if defined(__AVX2__) && defined(WITH_SIMD)
  s0 = s1;
  d0 = d1;
  n2 = (n >> 4);        // number of 16 token chunks
  s1 = s0 + (n2 << 3);  // "halfway" point
  d1 = d0 + (n2 << 3);  // "halfway" point
  ix = _mm256_setzero_si256();                         // unnecessary code to avoid compiler warning
  i0 = _mm_loadu_si128((__m128i const *) swapindex_16_32);
  ix = _mm256_inserti128_si256(ix, i0, 0);             // same shuffle pattern on low and high part
  ix = _mm256_inserti128_si256(ix, i0, 1);
  for( ; i < n - 15 ; i += 16){                        // endian swap of 16 32 bit tokens
    vs0 = _mm256_loadu_si256((__m256i const *) s0);    // load 8 32 bit tokens
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs0 = _mm256_shuffle_epi8 (vs0, ix);               // shuffle (vector bswap equivalent)
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d0, vs0);          // store 8 32 bit tokens
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s0 += 8;
    s1 += 8;
    d0 += 8;
    d1 += 8;
  }
  for( ; i < n - 7 ; i += 8){
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s1 += 8;
    d1 += 8;
  }
  for( ; i < n - 3 ; i += 4){
    v1 = _mm_loadu_si128((__m128i const *) s1);    // load 4 32 bit tokens
    v1 = _mm_shuffle_epi8 (v1, i0);                // shuffle (vector bswap equivalent)
    _mm_storeu_si128((__m128i *) d1, v1);          // store 4 32 bit tokens
    s1 += 4;
    d1 += 4;
  }
#endif
  for( ; i < n ; i++){    // loop over remainder
    t = *s1++;
    t = (t >> 16) | (t << 16) ;    // swap halfwords
    *d1++ = t;
  }
}

// endian swap of 8 bit tokens in 64 bit tokens (s -> d)
// s [IN]  : pointer to source of 64 bit elements (doublewords)
// d [OUT] : pointer to destination 64 bit elements (doublewords)
// n [IN]  : number of 64 bit elements (doublewords) to swap
// s and d may be the same address (swap INPLACE)
STATIC void inline Swap_8_in_64(void *s, void *d, int n){
#if defined(__AVX2__) && defined(WITH_SIMD)
  __m256i ix, vs0, vs1;
  __m128i i0, v1;
  int n2;
  uint32_t *s0, *d0;
#endif
  int i;
  uint32_t t1, t2;
  uint32_t *s1 = (uint32_t *) s;
  uint32_t *d1 = (uint32_t *) d;

  i = 0;
  n = n * 2;   // translate into number of 32 bit tokens
#if defined(__AVX2__) && defined(WITH_SIMD)
  s0 = s1;
  d0 = d1;
  n2 = (n >> 4);        // number of 16 token chunks
  s1 = s0 + (n2 << 3);  // "halfway" point
  d1 = d0 + (n2 << 3);  // "halfway" point
  ix = _mm256_setzero_si256();                         // unnecessary code to avoid compiler warning
  i0 = _mm_loadu_si128((__m128i const *) swapindex_8_64);
  ix = _mm256_inserti128_si256(ix, i0, 0);             // same shuffle pattern on low and high part
  ix = _mm256_inserti128_si256(ix, i0, 1);
  for( ; i < n - 15 ; i += 16){                        // endian swap of 16 32 bit tokens
    vs0 = _mm256_loadu_si256((__m256i const *) s0);    // load 8 32 bit tokens
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs0 = _mm256_shuffle_epi8 (vs0, ix);               // shuffle (vector bswap equivalent)
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d0, vs0);          // store 8 32 bit tokens
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s0 += 8;
    s1 += 8;
    d0 += 8;
    d1 += 8;
  }
  for( ; i < n - 7 ; i += 8){
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s1 += 8;
    d1 += 8;
  }
  for( ; i < n - 3 ; i += 4){
    v1 = _mm_loadu_si128((__m128i const *) s1);    // load 4 32 bit tokens
    v1 = _mm_shuffle_epi8 (v1, i0);                // shuffle (vector bswap equivalent)
    _mm_storeu_si128((__m128i *) d1, v1);          // store 4 32 bit tokens
    s1 += 4;
    d1 += 4;
  }
#endif
  for( ; i < n ; i=i+2){    // loop over remainder
    t1 = *s1++; t2 = *s1++;
    t1 = (t1 >> 24) | (t1 << 24) | ((t1 >> 8) & 0xFF00) | ((t1 << 8) & 0xFF0000);    // bswap
    t2 = (t2 >> 24) | (t2 << 24) | ((t2 >> 8) & 0xFF00) | ((t2 << 8) & 0xFF0000);    // bswap
    *d1++ = t2;  *d1++ = t1;
  }
}

// endian swap of 16 bit tokens in 64 bit tokens (s -> d)
// s [IN]  : pointer to source of 64 bit elements (doublewords)
// d [OUT] : pointer to destination 64 bit elements (doublewords)
// n [IN]  : number of 64 bit elements (doublewords) to swap
// s and d may be the same address (swap INPLACE)
STATIC void inline Swap_16_in_64(void *s, void *d, int n){
#if defined(__AVX2__) && defined(WITH_SIMD)
  __m256i ix, vs0, vs1;
  __m128i i0, v1;
  int n2;
  uint32_t *s0, *d0;
#endif
  int i;
  uint32_t t1, t2;
  uint32_t *s1 = (uint32_t *) s;
  uint32_t *d1 = (uint32_t *) d;

  i = 0;
  n = n * 2;   // translate into number of 32 bit tokens
#if defined(__AVX2__) && defined(WITH_SIMD)
  s0 = s1;
  d0 = d1;
  n2 = (n >> 4);        // number of 16 token chunks
  s1 = s0 + (n2 << 3);  // "halfway" point
  d1 = d0 + (n2 << 3);  // "halfway" point
  ix = _mm256_setzero_si256();                         // unnecessary code to avoid compiler warning
  i0 = _mm_loadu_si128((__m128i const *) swapindex_16_64);
  ix = _mm256_inserti128_si256(ix, i0, 0);             // same shuffle pattern on low and high part
  ix = _mm256_inserti128_si256(ix, i0, 1);
  for( ; i < n - 15 ; i += 16){                        // endian swap of 16 32 bit tokens
    vs0 = _mm256_loadu_si256((__m256i const *) s0);    // load 8 32 bit tokens
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs0 = _mm256_shuffle_epi8 (vs0, ix);               // shuffle (vector bswap equivalent)
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d0, vs0);          // store 8 32 bit tokens
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s0 += 8;
    s1 += 8;
    d0 += 8;
    d1 += 8;
  }
  for( ; i < n - 7 ; i += 8){
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 8 32 bit tokens
    vs1 = _mm256_shuffle_epi8 (vs1, ix);               // shuffle (vector bswap equivalent)
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 8 32 bit tokens
    s1 += 8;
    d1 += 8;
  }
  for( ; i < n - 3 ; i += 4){
    v1 = _mm_loadu_si128((__m128i const *) s1);    // load 4 32 bit tokens
    v1 = _mm_shuffle_epi8 (v1, i0);                // shuffle (vector bswap equivalent)
    _mm_storeu_si128((__m128i *) d1, v1);          // store 4 32 bit tokens
    s1 += 4;
    d1 += 4;
  }
#endif
  for( ; i < n ; i=i+2){    // loop over remainder
    t1 = *s1++; t2 = *s1++;
    t1 = (t1 >> 16) | (t1 << 16) ;    // swap halfwords
    t2 = (t2 >> 16) | (t2 << 16) ;    // swap halfwords
    *d1++ = t2;  *d1++ = t1;
  }
}

// endian swap of 32 bit tokens in 64 bit tokens (s -> d)
// no swap table is needed, an immediate argument is used
// s [IN]  : pointer to source of 64 bit elements (doublewords)
// d [OUT] : pointer to destination 64 bit elements (doublewords)
// n [IN]  : number of 64 bit elements (doublewords) to swap
// s and d may be the same address (swap INPLACE)
STATIC void Swap_32_in_64(void *s, void *d, int n){
#if defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vs0, vs1;
  int n2;
  uint64_t *s0, *d0;
#endif
  int i ;
  uint64_t t;
  uint64_t *s1 = (uint64_t *) s;
  uint64_t *d1 = (uint64_t *) d;

  i = 0;
#if defined(__AVX2__) && defined(WITH_SIMD)
  s0 = s1;
  d0 = d1;
  n2 = (n >> 3);        // number of full 8 token chunks
  s1 = s0 + (n2 << 2);  // "halfway" point
  d1 = d0 + (n2 << 2);  // "halfway" point
  for( ; i < n - 7 ; i +=8){                           // endian swap of 8 64 bit tokens
    vs0 = _mm256_loadu_si256((__m256i const *) s0);    // load 4 64 bit tokens
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 4 64 bit tokens
    vs0 = _mm256_shuffle_epi32(vs0 , 0xB1);            // 2 3 0 1  (0b10110001)
    vs1 = _mm256_shuffle_epi32(vs1 , 0xB1);            // 2 3 0 1  (0b10110001)
    _mm256_storeu_si256((__m256i *) d0, vs0);          // store 4 64 bit tokens
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 4 64 bit tokens
    s0 += 4;
    s1 += 4;
    d0 += 4;
    d1 += 4;
  }
  for( ; i < n - 3 ; i +=4){
    vs1 = _mm256_loadu_si256((__m256i const *) s1);    // load 4 64 bit tokens
    vs1 = _mm256_shuffle_epi32(vs1 , 0xB1);
    _mm256_storeu_si256((__m256i *) d1, vs1);          // store 4 64 bit tokens
    s1 += 4;
    d1 += 4;
  }
#endif
  for( ; i < n ; i++){    // loop over leftovers (everything if not AVX2)
    t = *s1++;
    t = (t >> 32) | (t << 32);
    *d1++ = t;
  }
}

// PROTOTYPES_ONLY
#endif

// FAST_ENDIAN
#endif
