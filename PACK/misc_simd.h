/*
 * extra SIMD functions for __x86_64__ (AVX and SSE2)
 * - logical not operator
 *
 * functions are implemented using the Intel® intrinsics
 */
#if ! defined(X86_SIMD_EXTRAS)
#define X86_SIMD_EXTRAS

#if ! defined(STATIC)
#define STATIC static
#define STATIC_DEFINED_HERE
#endif

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
#include <immintrin.h>
STATIC __m256i inline _mm256_not_si256(const __m256i v){
  return _mm256_xor_si256( v, _mm256_cmpeq_epi32(v,v) ) ;
}
#endif

#if defined(__x86_64__) && defined(__SSE2__) && defined(WITH_SIMD)
#include <emmintrin.h>
STATIC __m128i inline _mm_not_si128(const __m128i v){
  return _mm_xor_si128( v, _mm_cmpeq_epi32(v,v) ) ;
}
#endif

#if defined(STATIC_DEFINED_HERE)
#undef STATIC
#undef STATIC_DEFINED_HERE
#endif

#endif
