
#include <stdint.h>
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
#include <immintrin.h>
#endif

// static int b_conts = 0b01010101 ;

#define STATIC extern

// divide a signed integer by 2 with rounding toward +-infinity
#define IDIV2R(x) (( (x) + 1 + ((x)>>31) ) >> 1 )

// divide a signed integer by 4 with rounding toward +-infinity
#define IDIV4R(x) (( (x) + 2 + ((x)>>31) ) >> 2 )

// divide a signed integer by 8 with rounding toward +-infinity
#define IDIV8R(x) (( (x) + 4 + ((x)>>31) ) >> 3 )

// 2x2 averaging of 8 points on 2 lines (averaging along line and across lines)
// src1[8]   : first line
// src2[8]   : second line
// avg[4]    : result
STATIC inline void average_2x2_8_I32(int32_t * restrict src1, int32_t * restrict src2, int32_t * restrict avg){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i vs0, vs1, va1, va2, vb1, vb2, vk2 ;
  va1 = _mm_loadu_si128( (__m128i *)  src1   ) ;
  va2 = _mm_loadu_si128( (__m128i *) (src1+4)) ;
  vk2 = _mm_cmpeq_epi32(va1, va1) ;           // -1
  vk2 = _mm_add_epi32(vk2, vk2) ;             // -2
  vb1 = _mm_loadu_si128( (__m128i *)  src2   ) ;
  vb2 = _mm_loadu_si128( (__m128i *) (src2+4)) ;
  va1 = _mm_add_epi32(va1, vb1) ;             // add rows terms(0-3)
  va2 = _mm_add_epi32(va2, vb2) ;             // add rows terms(4-7)
  vs1 = _mm_hadd_epi32(va1,va2) ;            // add pairs of terms
  vs0 = _mm_sub_epi32(vs1, vk2) ;             // add 2 (subtract -2)
  vs1 = _mm_srai_epi32(vs1, 31) ;             // -1 or 0 according to sign
  vs0 = _mm_add_epi32(vs1, vs0) ;             // add to vs0
  vs0 = _mm_srai_epi32(vs0, 2) ;              // finalize divide by 4 with rounding
  _mm_storeu_si128( (__m128i *) avg, vs0) ;   // store result
#else
  int i, ii ;
  for(i=0, ii=0 ; ii<4 ; ii++, i+=2){
    avg[ii]  = src1[i] + src1[i+1] ;
    avg[ii] += src2[i] + src2[i+1] ;
    avg[ii]  = IDIV4R(avg[ii]) ;
  }
#endif
}

// 2x2 averaging of 16 points on 2 lines (averaging along line and across lines)
// src1[16]  : first line
// src2[16]  : second line
// avg[8]    : result
STATIC inline void average_2x2_16_I32(int32_t * restrict src1, int32_t * restrict src2, int32_t * restrict avg){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256i vs0, vs1, va1, va2, vb1, vb2, vk2 ;
  va1 = _mm256_loadu_si256( (__m256i *)  src1   ) ;
  vk2 = _mm256_cmpeq_epi32(va1, va1) ;           // -1
  vk2 = _mm256_add_epi32(vk2, vk2) ;             // -2
  va2 = _mm256_loadu_si256( (__m256i *) (src1+8)) ;
  vb1 = _mm256_loadu_si256( (__m256i *)  src2   ) ;
  vb2 = _mm256_loadu_si256( (__m256i *) (src2+8)) ;
  va1 = _mm256_add_epi32(va1, vb1) ;             // add rows terms(0-7)
  va2 = _mm256_add_epi32(va2, vb2) ;             // add rows terma(8-15)
  vs1 = _mm256_hadd_epi32(va1, va2) ;            // add pairs of terms
  vs1 = _mm256_permute4x64_epi64(vs1, 0b11011000) ;
  vs0 = _mm256_sub_epi32(vs1, vk2) ;             // add 2 (subtract -2)
  vs1 = _mm256_srai_epi32(vs1, 31) ;             // -1 or 0 according to sign
  vs0 = _mm256_add_epi32(vs1, vs0) ;             // add to vs0
  vs0 = _mm256_srai_epi32(vs0, 2) ;              // finalize divide by 4 with rounding
  _mm256_storeu_si256( (__m256i *) avg, vs0) ;   // store result
#else
  int i, ii ;
  for(i=0, ii=0 ; ii<8 ; ii++, i+=2){
    avg[ii]  = src1[i] + src1[i+1] ;
    avg[ii] += src2[i] + src2[i+1] ;
    avg[ii]  = IDIV4R(avg[ii]) ;
  }
#endif
}

// 2x2 averaging of n points on 2 lines (averaging along line and across lines)
// src1[n]      : first line
// src2[n]      : second line
// avg[(n+1)/2] : result
// n            : number of points (may be odd)
STATIC inline void average_2x2_I32(int32_t * restrict src1, int32_t * restrict src2, int32_t * restrict avg, uint32_t n){
  int i, ii ;
  int n2 = n>>1 ;
  int n1 = (n2 & 7) ? (n2 & 7) : 8 ;

  if(n1 <= 4) average_2x2_8_I32(src1, src2, avg) ;
  else        average_2x2_16_I32(src1, src2, avg) ;
  for(i=n1+n1, ii=n1 ; ii<n2 ; ii+=8, i+=16){
    average_2x2_16_I32(src1+i, src2+i, avg+ii) ;
  }
  if(n & 1) {    // odd number of points in lines, 2 points are averaged instead of 4
    avg[ii] = src1[i] + src2[i] ;
    avg[ii] = IDIV2R(avg[ii]) ;
  }
}

// average one row
// src1[n]      : row to be averaged
// avg[(n+1)/2] : result
// n            : number of points (may be odd)
// NOTE : average_2x2_I32 is used, simulating 2 identical rows
STATIC inline void average_2x1_I32(int32_t * restrict src, int32_t * restrict avg, uint32_t n){
  int i, ii ;
  int n2 = n>>1 ;
  average_2x2_I32(src, src, avg, n) ;
}

// 2x2 averaging of a 2 dimensional array
// src[lni       *        nj] : source array
// avg[((ni+1)/2)*((nj+1)/2)] : result
// lni  : storage dimension of rows in src
// ni   : number of points to be averaged in rows (may be odd)
// nj   : number of rows to be averaged (may be odd)
void average_2x2_2D_I32(int32_t * restrict src, int32_t * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj){
  int i, j ;
  int ni2 = (ni+1)/2 ;
  int lni2 = lni+lni ;
  for(j=0 ; j<nj/2 ; j++){
    average_2x2_I32(src, src+lni, avg, ni) ;
    src += lni2 ;
    avg += ni2 ;
  }
  if(nj & 1) average_2x2_I32(src, src, avg, ni) ;
}

// 2x2 average from src(ni,nj) -> avg((ni+1)/2,(nj+1)/2)
// residual -> res(ni,nj)
void average_2x2_2D_res_I32(int32_t * restrict src, int32_t * restrict avg, int32_t * restrict res, int ni, int lni, int nj){
  int i, ii, j ;
  for (j=0 ; j<nj-1 ; j+=2){
    for(i=0, ii=0 ; i<ni-1 ; i+=2, ii++){
      avg[ii] = (src[i] + src[i+1] + src[i+lni] +src[i+lni+1]) / 4 ;
      res[i]       = src[i]       - avg[ii] ;
      res[i+1]     = src[i+1]     - avg[ii] ;
      res[i+lni]   = src[i+lni]   - avg[ii] ;
      res[i+lni+1] = src[i+lni+1] - avg[ii] ;
    }
    if(i < ni) {   // odd number of points in row
      avg[ii] = (src[i] + src[i+lni]) / 2 ;
      res[i]       = src[i]       - avg[ii] ;
      res[i+lni]   = src[i+lni]   - avg[ii] ;
    }
    src += lni ;
    res += ni ;
    avg += (ni+1)/2 ;
  }
  if(j < nj){   // odd number of rows
    for(i=0, ii=0 ; i<ni-1 ; i+=2, ii++){
      avg[ii] = (src[i] + src[i+1]) / 2 ;
      res[i]       = src[i]       - avg[ii] ;
      res[i+1]     = src[i+1]     - avg[ii] ;
    }
    if(i < ni) {   // odd number of points in last row
      avg[ii] = src[i] ;
      res[i]  = 0 ;
    }
  }
}

// restore along i
STATIC void expand_2x2_row_along_i(int32_t * restrict dst, int32_t * restrict avg, int n){
  int i, ii ;
  int n2 = n>>1 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
#endif

  dst[0] = 5*avg[0] - avg[1] ; dst[0] = IDIV4R(dst[0]) ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  int nm2 = n2-1, i2 = 1 ;
  int n7 = (nm2 & 7) ? (nm2 & 7) : 8 ;
  __m256i va1, va2, va3, vb1, vb2, vb3, v31, v13, q31, q13, vlo, vhi, vk2, v07, v8f ;
  for(i=0 ; i<nm2 ; i+= n7, n7 = 8) {
    va1 = _mm256_loadu_si256( (__m256i *) (avg + i   ) );
    vb1 = _mm256_loadu_si256( (__m256i *) (avg + i +1) );
    vk2 = _mm256_cmpeq_epi32(va1, va1) ;                // -1
    vk2 = _mm256_add_epi32(vk2, vk2) ;                  // -2
    va2 = _mm256_add_epi32(va1, va1) ;                  // va1 * 2
    vb2 = _mm256_add_epi32(vb1, vb1) ;                  // vb1 * 2
    va3 = _mm256_add_epi32(va2, va1) ;                  // va1 * 3
    vb3 = _mm256_add_epi32(vb2, vb1) ;                  // vb1 * 3
    v31 = _mm256_add_epi32(va3, vb1) ;                  // 3*row0[i] + 1*row1[i]
    v13 = _mm256_add_epi32(va1, vb3) ;                  // 1*row0[i] + 3*row1[i]
    q31 = _mm256_srai_epi32(v31, 31) ;                  // -1 or 0 according to sign
    q13 = _mm256_srai_epi32(v13, 31) ;                  // -1 or 0 according to sign
    q31 = _mm256_sub_epi32(q31, vk2) ;                  // q31 + 2
    q13 = _mm256_sub_epi32(q13, vk2) ;                  // q13 + 2
    q31 = _mm256_add_epi32(v31, q31) ;                  // + v31
    q13 = _mm256_add_epi32(v13, q13) ;                  // + v13
    q31 = _mm256_srai_epi32(q31, 2) ;                   // IDIV4R(v31)  (dst[1, 3, 5, ...])
    q13 = _mm256_srai_epi32(q13, 2) ;                   // IDIV4R(v13)  (dst[2, 4, 6, ...])
    vlo = _mm256_unpacklo_epi32(q31, q13) ;             // 0 1 2 3 8 9 A B
    vhi = _mm256_unpackhi_epi32(q31, q13) ;             // 4 5 6 7 C D E F
    v07 = _mm256_permute2x128_si256(vlo, vhi, 0x20) ;   // 0 1 2 3 4 5 6 7 (dst[1, 2, 3, ...])
    v8f = _mm256_permute2x128_si256(vlo, vhi, 0x31) ;   // 8 9 A B C D E F (dst(9, A, B, ...])
    _mm256_storeu_si256( (__m256i *) (dst + i2    ), v07 );
    _mm256_storeu_si256( (__m256i *) (dst + i2 + 8), v8f );
    i2 = i2 + n7 + n7 ;
  }
#else
  int t1[n2], t2[n2] ;
  for(ii=0 ; ii<n2-1 ; ii++) {
    t1[ii] = 3*avg[ii] + 1*avg[ii+1] ;
    t1[ii] = IDIV4R(t1[ii]) ;
    t2[ii] = 1*avg[ii] + 3*avg[ii+1] ;
    t2[ii] = IDIV4R(t2[ii]) ;
  }
  for(i=1, ii=0 ; i<n-2 ; i+=2, ii++){
    dst[i  ] = t1[ii] ; dst[i+1] = t2[ii] ;
  }
#endif
  dst[n2+n2-1] = 5*avg[n2-1] - avg[n2-2] ; dst[n2+n2-1] = IDIV4R(dst[n2+n2-1]) ;
  if(n & 1) dst[n2+n2] = avg[n2] ;
}

// bottom row
STATIC void expand_2x2_row_0(int32_t * restrict row0, int32_t * restrict row1, int32_t * src0, int n){
  int i ;
  int ni2 = (n+1)>>1 ;
  int32_t a51[ni2] ;
  for(i=0 ; i<ni2 ; i++) {
    a51[i] = 5*row0[i] - 1*row1[i] ;  // unaverage along j to restore src0
    a51[i] = IDIV4R(a51[i]) ;
    expand_2x2_row_along_i(src0, a51, n) ;  // restore src0 along i
  }
}

// pairs of middle rows
STATIC void expand_2x2_2_rows(int32_t * restrict row0, int32_t * restrict row1, int32_t * src0, int32_t * src1, int n){
  int i ;
  int ni2 = (n+1)>>1 ;
  int32_t a31[ni2], a13[ni2] ;

  for(i=0 ; i<ni2 ; i++) {
    a31[i] = 3*row0[i] + 1*row1[i] ;  // unaverage along j to restore src0
    a31[i] = IDIV4R(a31[i]) ;
    a13[i] = 1*row0[i] + 3*row1[i] ;  // unaverage along j to restore src1
    a13[i] = IDIV4R(a13[i]) ;
  }
  expand_2x2_row_along_i(src0, a31, n) ;  // restore src0 along i
  expand_2x2_row_along_i(src1, a13, n) ;  // restore src1 along i
}

// top row if even number of rows
STATIC void expand_2x2_row_n(int32_t * restrict row0, int32_t * restrict row1, int32_t * src0, int n){
  int i ;
  int ni2 = (n+1)>>1 ;
  int32_t a15[ni2] ;
  for(i=0 ; i<ni2 ; i++) {
    a15[i] = 5*row1[i] - 1*row0[i] ;  // unaverage along j to restore src0
    a15[i] = IDIV4R(a15[i]) ;
    expand_2x2_row_along_i(src0, a15, n) ;  // restore src0 along i
  }
}

// restore a previously averaged array using linear interpolation/extrapolation
// dst[lni       *        nj] : result
// avg[((ni+1)/2)*((nj+1)/2)] : previously averaged array
// lni  : storage dimension of rows in dst
// ni   : number of averaged points in rows (may be odd)
// nj   : number of averaged rows (may be odd)
STATIC inline void expand_2x2_2D(int32_t * restrict dst, int32_t * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj){
  int i, j ;
  int ni2 = (ni+1)/2 ;
  int lni2 = lni+lni ;
  int32_t t[ni] ;

  // first row
  expand_2x2_row_0(avg, avg+ni2, dst, ni) ;
  dst += lni ;
  for(j=1 ; j<nj/2 ; j++){        // pairs of rows
    expand_2x2_2_rows(avg, avg+ni2, dst, dst+lni, ni) ;
    avg += ni2 ;
    dst += lni2 ;
  }
  // last row
  expand_2x2_row_n(avg-ni2, avg, dst, ni) ;
  if(nj & 1) expand_2x2_row_along_i(dst+lni, avg+ni2, ni) ;
}

// a[8] , b[8] -> c[4] , row length of a = 8
STATIC inline void average_rows_8x2_8(float *a, float *c){
  int i, j ;
  float t[8] ;
  for(i = 0 ; i < 8 ; i++){           // add rows
    t[i] = a[i] + a[i+8] ;
  }
  for(i =0, j=0 ; j<4 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[8,2] , b[8,2] -> c[4] , row length of a = lrow
STATIC inline void average_rows_8x2(float *a, float *c, int lrow){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 vr ;
  __m128 va, vb, v0, v1 ;
  vr = _mm256_add_ps ( _mm256_loadu_ps(a) , _mm256_loadu_ps(a + lrow) ) ;
  va = _mm256_extractf128_ps(vr, 0) ;
  vb = _mm256_extractf128_ps(vr, 1) ;
  v0 = _mm_shuffle_ps(va, vb, 136) ;
  v1 = _mm_shuffle_ps(va, vb, 221) ;
  v0 = _mm_add_ps(v0, v1) ;
  v0 = _mm_mul_ps(v0, _mm_set1_ps(.25f)) ;
  _mm_storeu_ps(c, v0) ;
#else
  int i, j ;
  float t[8] ;
  for(i = 0 ; i < 8 ; i++){           // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<4 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
#endif
}

// c[4] -> a[8,2] , b[8,2] , row length of a = lrow
STATIC inline void expand_rows_8x2(float *a, float *c, int lrow){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128 v2, va, vb ;
  __m256 v ;
  v2 = _mm_loadu_ps(c) ;
  va = _mm_unpacklo_ps(v2, v2) ;  // 1 1 0 0
  vb = _mm_unpackhi_ps(v2, v2) ;  // 3 3 2 2
  v = _mm256_insertf128_ps(v, va, 0) ;
  v = _mm256_insertf128_ps(v, vb, 1) ;
  _mm256_storeu_ps(a     , v) ;
  _mm256_storeu_ps(a+lrow, v) ;
#else
  int i, j ;
  for(i =0, j=0 ; j<4 ; j++, i+=2){
    a[i] = a[i+1] = a[i+lrow] = a[i+lrow+1] = c[j] ;
  }
#endif
}

// a[8] , b[8] -> c[4] , row length of a = lrow
STATIC inline void average_rows_8x2_avx2(float *a, float *c, int lrow){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 vr ;
  __m128 va, vb, v0, v1 ;
  vr = _mm256_add_ps ( _mm256_loadu_ps(a) , _mm256_loadu_ps(a + lrow) ) ;
  va = _mm256_extractf128_ps(vr, 0) ;
  vb = _mm256_extractf128_ps(vr, 1) ;
  v0 = _mm_shuffle_ps(va, vb, 136) ;
  v1 = _mm_shuffle_ps(va, vb, 221) ;
  v0 = _mm_add_ps(v0, v1) ;
  v0 = _mm_mul_ps(v0, _mm_set1_ps(.25f)) ;
  _mm_storeu_ps(c, v0) ;
#else
#endif
}

// a[16,2] -> c[8] , row length of a = 16
STATIC inline void average_rows_16x2_16(float *a, float *c){
  int i, j ;
  float t[16] ;
  for(i = 0 ; i < 16 ; i++){           // add rows
    t[i] = a[i] + a[i+16] ;
  }
  for(i =0, j=0 ; j<8 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[16,2] -> c[8] , row length of a = lrow
STATIC inline void average_rows_16x2(float *a, float *c, int lrow){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 vr0, vr1 ;
  __m256 vlo, vhi, v0, v1 ;
  vr0 = _mm256_add_ps ( _mm256_loadu_ps(a  ) , _mm256_loadu_ps(a + lrow    ) ) ;
  vr1 = _mm256_add_ps ( _mm256_loadu_ps(a+8) , _mm256_loadu_ps(a + lrow + 8) ) ;
  vlo = _mm256_permute2f128_ps(vr0, vr1, 32) ;
  vhi = _mm256_permute2f128_ps(vr0, vr1, 49) ;
  v0 = _mm256_shuffle_ps(vlo, vhi, 136) ;
  v1 = _mm256_shuffle_ps(vlo, vhi, 221) ;
  v0 = _mm256_add_ps (v0, v1) ;
  v0 = _mm256_mul_ps (v0, _mm256_set1_ps(.25f) ) ;
  _mm256_storeu_ps(c, v0) ;
#else
  int i, j ;
  float t[16] ;
  for(i = 0 ; i < 16 ; i++){           // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<8 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
#endif
}

// c[8] -> a[16,2] , b[16,2] , row length of a = lrow
STATIC inline void expand_rows_16x2(float *a, float *c, int lrow){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 v, va, vb, vc, vd ;
  v  = _mm256_loadu_ps(c) ;                  // 7 6 5 4 3 2 1 0
  va = _mm256_unpacklo_ps(v, v) ;            // 5 5 4 4 1 1 0 0
  vb = _mm256_unpackhi_ps(v, v) ;            // 7 7 6 6 3 3 2 2
  vc = _mm256_permute2f128_ps(va, vb, 32) ;  // 3 3 2 2 1 1 0 0
  vd = _mm256_permute2f128_ps(va, vb, 49) ;  // 7 7 6 6 5 5 4 4
  _mm256_storeu_ps(a       , vc) ;
  _mm256_storeu_ps(a+8     , vd) ;
  _mm256_storeu_ps(a+lrow  , vc) ;
  _mm256_storeu_ps(a+lrow+8, vd) ;
#else
  int i, j ;
  for(i =0, j=0 ; j<8 ; j++, i+=2){
    a[i] = a[i+1] = a[i+lrow] = a[i+lrow+1] = c[j] ;
  }
#endif
}

// a[16,2] -> c[8] , row length of a = lrow
STATIC inline void average_rows_16x2_avx2(float *a, float *c, int lrow){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 vr0, vr1 ;
  __m256 vlo, vhi, v0, v1 ;
  vr0 = _mm256_add_ps ( _mm256_loadu_ps(a  ) , _mm256_loadu_ps(a + lrow    ) ) ;
  vr1 = _mm256_add_ps ( _mm256_loadu_ps(a+8) , _mm256_loadu_ps(a + lrow + 8) ) ;
  vlo = _mm256_permute2f128_ps(vr0, vr1, 32) ;
  vhi = _mm256_permute2f128_ps(vr0, vr1, 49) ;
  v0 = _mm256_shuffle_ps(vlo, vhi, 136) ;
  v1 = _mm256_shuffle_ps(vlo, vhi, 221) ;
  v0 = _mm256_add_ps (v0, v1) ;
  v0 = _mm256_mul_ps (v0, _mm256_set1_ps(.25f) ) ;
  _mm256_storeu_ps(c, v0) ;
#else
#endif
}

// a[32,2] -> c[16] , row length of a = 32
STATIC inline void average_rows_32x2_32(float *a, float *c){
  int i, j ;
  float t[32] ;
  for(i = 0 ; i < 32 ; i++){           // add rows
    t[i] = a[i] + a[i+32] ;
  }
  for(i =0, j=0 ; j<16 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[32,2] -> c[16] , row length of a = lrowlength of a = lrow
STATIC inline void average_rows_32x2(float *a, float *c, int lrow){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 vr0, vr1 ;
  __m256 vlo, vhi, v0, v1 ;
  vr0 = _mm256_add_ps ( _mm256_loadu_ps(a  ) , _mm256_loadu_ps(a + lrow    ) ) ;
  vr1 = _mm256_add_ps ( _mm256_loadu_ps(a+8) , _mm256_loadu_ps(a + lrow + 8) ) ;
  vlo = _mm256_permute2f128_ps(vr0, vr1, 32) ;
  vhi = _mm256_permute2f128_ps(vr0, vr1, 49) ;
  v0 = _mm256_shuffle_ps(vlo, vhi, 136) ;
  v1 = _mm256_shuffle_ps(vlo, vhi, 221) ;
  v0 = _mm256_add_ps (v0, v1) ;
  v0 = _mm256_mul_ps (v0, _mm256_set1_ps(.25f) ) ;
  _mm256_storeu_ps(c, v0) ;

  a += 16 ; c += 8 ;
  vr0 = _mm256_add_ps ( _mm256_loadu_ps(a  ) , _mm256_loadu_ps(a + lrow    ) ) ;
  vr1 = _mm256_add_ps ( _mm256_loadu_ps(a+8) , _mm256_loadu_ps(a + lrow + 8) ) ;
  vlo = _mm256_permute2f128_ps(vr0, vr1, 32) ;
  vhi = _mm256_permute2f128_ps(vr0, vr1, 49) ;
  v0 = _mm256_shuffle_ps(vlo, vhi, 136) ;
  v1 = _mm256_shuffle_ps(vlo, vhi, 221) ;
  v0 = _mm256_add_ps (v0, v1) ;
  v0 = _mm256_mul_ps (v0, _mm256_set1_ps(.25f) ) ;
  _mm256_storeu_ps(c, v0) ;
#else
  int i, j ;
  float t[32] ;
  for(i = 0 ; i < 32 ; i++){           // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<16 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
#endif
}

// c[16] -> a[32,2] , b[32,2] , row length of a = lrow
STATIC inline void expand_rows_32x2(float *a, float *c, int lrow){
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m256 v, va, vb, vc, vd ;
  v  = _mm256_loadu_ps(c) ;                  // 7 6 5 4 3 2 1 0
  va = _mm256_unpacklo_ps(v, v) ;            // 5 5 4 4 1 1 0 0
  vb = _mm256_unpackhi_ps(v, v) ;            // 7 7 6 6 3 3 2 2
  vc = _mm256_permute2f128_ps(va, vb, 32) ;  // 3 3 2 2 1 1 0 0
  vd = _mm256_permute2f128_ps(va, vb, 49) ;  // 7 7 6 6 5 5 4 4
  _mm256_storeu_ps(a       , vc) ;
  _mm256_storeu_ps(a+8     , vd) ;
  _mm256_storeu_ps(a+lrow  , vc) ;
  _mm256_storeu_ps(a+lrow+8, vd) ;
  a += 16 ;                                  // next batch of 8
  c += 8 ;
  v  = _mm256_loadu_ps(c) ;                  // 7 6 5 4 3 2 1 0
  va = _mm256_unpacklo_ps(v, v) ;            // 5 5 4 4 1 1 0 0
  vb = _mm256_unpackhi_ps(v, v) ;            // 7 7 6 6 3 3 2 2
  vc = _mm256_permute2f128_ps(va, vb, 32) ;  // 3 3 2 2 1 1 0 0
  vd = _mm256_permute2f128_ps(va, vb, 49) ;  // 7 7 6 6 5 5 4 4
  _mm256_storeu_ps(a       , vc) ;
  _mm256_storeu_ps(a+8     , vd) ;
  _mm256_storeu_ps(a+lrow  , vc) ;
  _mm256_storeu_ps(a+lrow+8, vd) ;
#else
  int i, j ;
  for(i =0, j=0 ; j<16 ; j++, i+=2){
    a[i] = a[i+1] = a[i+lrow] = a[i+lrow+1] = c[j] ;
  }
#endif
}

// a[64,2] -> c[32] , row length of a = 64
STATIC inline void average_rows_64x2_64(float *a, float *c){
  int i, j ;
  float t[64] ;
  for(i = 0 ; i < 64 ; i++){           // add rows
    t[i] = a[i] + a[i+64] ;
  }
  for(i =0, j=0 ; j<32 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[64,2] -> c[32] , row length of a = lrow
STATIC inline void average_rows_64x2(float *a, float *c, int lrow){
  average_rows_32x2(a     , c     , lrow) ;
  average_rows_32x2(a + 32, c + 16, lrow) ;
}

// c[32] -> a[64,2] , b[64,2] , row length of a = lrow
STATIC inline void expand_rows_64x2(float *a, float *c, int lrow){
  expand_rows_32x2(a   , c   , lrow) ;
  expand_rows_32x2(a+32, c+16, lrow) ;
}

// a[n,2] -> c[n/2] , row length of a = lrow
STATIC inline void average_rows_nx2(float *a, float *c, int n, int lrow){
  int i, j ;
  float t[8] ;
  while(n >= 64){                        // slices of [64,2] -> [32]
    average_rows_64x2(a, c, lrow) ;
    n -= 64 ; a+= 64 ; c += 32 ;         // next slice
  }
  if(n >= 32){                           // slices of [32,2] -> [16]
    average_rows_32x2(a, c, lrow) ;
    n -= 32 ; a+= 32 ; c += 16 ;         // next slice
  }
  if(n >= 16){                           // slices of [16,2] -> [8]
    average_rows_16x2(a, c, lrow) ;
    n -= 16 ; a+= 16 ; c += 8 ;          // next slice
  }
  if(n >= 8){                            // slices of [8,2] -> [4]
    average_rows_8x2(a, c, lrow) ;
    n -= 8 ; a+= 8 ; c += 4 ;            // next slice
  }
  // process last (<8) points
  for(i = 0 ; i < n ; i++){              // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<n/2 ; j++, i+=2){    // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// c[n/2] -> a[n,2] , row length of a = lrow
STATIC inline void expand_rows_nx2(float *a, float *c, int n, int lrow){
  int i, j ;
  while(n >= 64){                        // slices of [64,2] <- [32]
    expand_rows_64x2(a, c, lrow) ;
    n -= 64 ; a+= 64 ; c += 32 ;         // next slice
  }
  if(n >= 32){                           // slices of [32,2] <- [16]
    expand_rows_32x2(a, c, lrow) ;
    n -= 32 ; a+= 32 ; c += 16 ;         // next slice
  }
  if(n >= 16){                           // slices of [16,2] <- [8]
    expand_rows_16x2(a, c, lrow) ;
    n -= 16 ; a+= 16 ; c += 8 ;          // next slice
  }
  if(n >= 8){                            // slices of [8,2] <- [4]
    expand_rows_8x2(a, c, lrow) ;
    n -= 8 ; a+= 8 ; c += 4 ;            // next slice
  }
  for(i =0, j=0 ; j<n/2 ; j++, i+=2){
    a[i] = a[i+1] = a[i+lrow] = a[i+lrow+1] = c[j] ;
  }
}

// c[32,4] -> a[64,8]
STATIC inline void expand_rows_64x8(float *a, float *c, int lrow){
  expand_rows_64x2(a    , c   , lrow) ;
  expand_rows_64x2(a+128, c+32, lrow) ;
  expand_rows_64x2(a+256, c+64, lrow) ;
  expand_rows_64x2(a+384, c+96, lrow) ;
}

// c[32,32] -> a[64,64]
STATIC inline void expand_rows_64x64(float *a, float *c, int lrow){
  int i ;
  for(i = 0 ; i < 8 ; i++){
    expand_rows_64x8(a + i*64*8, c + i*32*4, lrow) ;
  }
}

// a[64,8] -> a2[32,4] -> a4[16,2] -> a8[8] , row length of a = ltow
STATIC inline void averages_64x8(float *a, float *a2, float *a4, float *a8, int lrow){
  // 8 lines of a -> 4 lines of a2
  average_rows_64x2(a       , a2     , lrow) ;  // [64,2] -> [32]
  average_rows_64x2(a  + 128, a2 + 32, lrow) ;  // [64,2] -> [32]
  average_rows_64x2(a  + 256, a2 + 64, lrow) ;  // [64,2] -> [32]
  average_rows_64x2(a  + 384, a2 + 96, lrow) ;  // [64,2] -> [32]
  // 4 lines of a2 -> 2 lines of a4
  average_rows_32x2(a2      , a4     , 32) ;  // [32,2] -> [16]
  average_rows_32x2(a2 +  64, a4 + 16, 32) ;  // [32,2] -> [16]
  // 2 lines of a4 -> 1 line of a8
  average_rows_16x2(a4      , a8     , 16) ;  // [16,2] -> [8]
}

// a[64,64] a2[32,32]  a4[16,16]  a8[8,8] , row length of a = lrow
STATIC inline void averages_64x64(float *a, float *a2, float *a4, float *a8, int lrow){
  int i ;
  for(i = 0 ; i < 8 ; i++){   // 8 blocks of a[64,8]  a2[32,4]  a4[16,2]  a8[8]
    averages_64x8(a + i*64*8, a2 + i*32*4, a4 + i*16*2, a8 + i*8, lrow) ;
  }
}
