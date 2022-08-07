
#include <stdint.h>
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
#include <immintrin.h>
#endif

// static int b_conts = 0b01010101 ;

#define STATIC extern

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
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
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
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
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
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
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
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
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
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
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
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
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
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
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
#if defined(__x86_64__) && defined(__AVX2__) && defined(SIMD)
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
