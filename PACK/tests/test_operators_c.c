#include <stdint.h>
#include <stdio.h>

#include <misc_operators.h>
#define NP 33

static int32_t xt2[] = {-8, -7, -7, -6, -6, -5, -5, -4, -4, -3, -3, -2, -2, -1, -1,  0,  0,  0,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8} ;
static int32_t xt4[] = {-4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4} ;
static int32_t xt8[] = {-2, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  2} ;
static int32_t xr2[] = {-8, -8, -7, -7, -6, -6, -5, -5, -4, -4, -3, -3, -2, -2, -1, -1,  0,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8} ;
static int32_t xr4[] = {-4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4} ;
static int32_t xr8[] = {-2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2} ;
static int32_t xi[]  = { 999, 0, 1, 2, 3, 4, 5, 6, 7, 999 } ;

static int compare(int32_t *a, int32_t *b, int n){
  int i, e = 0 ;
  for(i=0 ; i<n ; i++) e += ( (a[i] == b[i]) ? 0 : 1 ) ;
  return e ;
}

static void e_exit(int n){
  printf("ERROR in test, aborting\n") ;
  exit(n) ;
}

int main(int argc, char **argv){
  int32_t src[NP], dst[NP+8] ;
  int32_t i, j, x, err ;
  float f ;
  float xf[NP] ;
#if defined(__x86_64__) && defined(__AVX2__)
  __m128i v128 ;
  __m256i v256 ;
  __m128i v128lo, v128hi ;
#endif
// IEEE min and max exponent test
  printf("IEEE exponents : ") ;
  for(i=0 ; i<NP ; i++) xf[i] = (i - NP/2) * 1.01f ;
  uint32_t emax = ieee_max_exponent(xf, NP) ;
  uint32_t emin = ieee_min_exponent(xf, NP) ;
  uint32_t eminmax = ieee_minmax_exponent(xf, NP) ;
  printf("x = %f -> %f -> %f -> %f -> %f, emin = %d, emax = %d, eminmax = %4.4x : ", 
         xf[0], xf[NP/2-1], xf[NP/2], xf[NP/2+1], xf[NP-1], emin, emax, eminmax) ;
  if(emin != 127 || emax != 131 || emin != (eminmax >> 8) || emax != (eminmax & 0xFF) ) {
    printf("expecting emax = 131, got %d, emin = 127, got %d, eminmax = 7f83, got %8.8x\n",
           emax, emin, eminmax) ;
    e_exit(1) ;
  }
  printf("Success\n") ;
// test float to int conversion macro and function(rounded)
  printf("float to int rounding test : ") ;
  err = 0 ;
  for(f=-3.5f ; f <= 3.5f ; f+=.25f){
    if( (f < 0.0f) && ( FLOAT_TO_INT(f) != (int)(f - .5f) ) ) err++ ;
    if( (f > 0.0f) && ( FLOAT_TO_INT(f) != (int)(f + .5f) ) ) err++ ;
    if( (f == 0.0f) && ( FLOAT_TO_INT(f) != 0 ) ) err++ ;
    if(err){
      printf("float = %10f, int = %3d, %3d\n", f, FLOAT_TO_INT(f), float_to_int(f)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
// extract 128 bits from 256 bit register
#if defined(__x86_64__) && defined(__AVX2__)
  printf("extract upper/lower 128 bits : ");
  v256 = _mm256_loadu_si256((__m256i *) (xi + 1)) ;
  v128lo = _mm_lower128(v256) ;
  v128hi = _mm_upper128(v256) ;
  _mm_storeu_si128((__m128i *)(dst  ), v128lo) ;
  _mm_storeu_si128((__m128i *)(dst+4), v128hi) ;
  for(i=0 ; i<8 ; i++) {
    if(i != dst[i]){
      printf("expecting %d, got %d\n", i, dst[i]) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
#endif
// test divide and truncate or round macros
  for(i = 0 ; i < NP ; i++) src[i] = i - NP/2;
  printf("ORIG  ") ; for(i = 0 ; i < NP ; i++) printf("%4d", src[i]) ; printf("\n\n") ;
  // DIV2T/DIV4T/DIV8T divide by 2/4/8 and truncate (round toward 0)
  // DIV2R/DIV4R/DIV8R divide by 2/4/8 and round toward +infinity or -infinity as appropriate
  // scalar and X86 SIMD versions (128 and 256 bit, 4 and 8 elements)
#if defined(__x86_64__) && defined(__AVX2__)
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV2T(src[i]) ;
  printf("IDIV2T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt2, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV2T_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-2TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt2, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV2T_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-2TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt2, NP)) ; if(err) e_exit(1) ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV2R(src[i]) ;
  printf("IDIV2R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr2, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV2R_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-2RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr2, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV2R_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-2RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr2, NP)) ; if(err) e_exit(1) ;
  printf("\n") ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV4T(src[i]) ;
  printf("IDIV4T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", IDIV4T(src[i])) ; printf(" err = %d\n", err = compare(dst, xt4, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV4T_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-4TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt4, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV4T_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-4TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt4, NP)) ; if(err) e_exit(1) ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV4R(src[i]) ;
  printf("IDIV4R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", IDIV4R(src[i])) ; printf(" err = %d\n", err = compare(dst, xr4, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV4R_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-4RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr4, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), _mm256_idiv4r_epi32(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-4RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr4, NP)) ; if(err) e_exit(1) ;
  printf("\n") ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV8T(src[i]) ;
  printf("IDIV8T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt8, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV8T_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-8TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt8, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV8T_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-8TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt8, NP)) ; if(err) e_exit(1) ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV8R(src[i]) ;
  printf("IDIV8R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr8, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV8R_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-8RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr8, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV8R_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-8RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr8, NP)) ; if(err) e_exit(1) ;
  printf("\n") ;

#else
  printf("/2     ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", src[i]/2) ; printf("\n") ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV2T(src[i]) ;
  printf("IDIV2T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt2, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV2R(src[i]) ;
  printf("IDIV2R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr2, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV4T(src[i]) ;
  printf("IDIV4T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", IDIV4T(src[i])) ; printf(" err = %d\n", err = compare(dst, xt4, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV4R(src[i]) ;
  printf("IDIV4R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", IDIV4R(src[i])) ; printf(" err = %d\n", err = compare(dst, xr4, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV8T(src[i]) ;
  printf("IDIV8T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xt8, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV8R(src[i]) ;
  printf("IDIV8R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare(dst, xr8, NP)) ; if(err) e_exit(1) ;
  printf("\n") ;
#endif
// test mask macros
  int32_t mask32 = -1 ;
  int64_t mask64 = -1l ;
  printf("LMASK32 test : ") ;
  for(i=32, mask32=-1 ; i >= 0 ; i-- , mask32 <<= 1){
    if(mask32 != (LMASK32(i))) {
      printf("LMASK32(%2d) expected %8.8x, got %8.8x\n", i, mask32, LMASK32(i)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("LMASK32Z test : ") ;
  for(i=32, mask32=-1 ; i > 0 ; i-- , mask32 <<= 1){
    if(mask32 != (LMASK32Z(i))) {
      printf("LMASK32Z(%2d) expected %8.8x, got %8.8x\n", i, mask32, LMASK32Z(i)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("RMASK32 test : ") ;
  for(i=0, mask32=0 ; i<=32 ; i++, mask32 = (mask32 << 1) | 1){
    if(mask32 != RMASK32(i)){
      printf("RMASK32(%2d) expected %8.8x, got %8.8x\n", i, mask32, RMASK32(i)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("RMASK31 test : ") ;
  for(i=0, mask32=0 ; i<32 ; i++, mask32 = (mask32 << 1) | 1){
    if(mask32 != RMASK31(i)){
      printf("RMASK31(%2d) expected %8.8x, got %8.8x\n", i, mask32, RMASK31(i)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("LMASK64 test : ") ;
  for(i=64, mask64=-1l ; i >= 0 ; i-- , mask64 <<= 1){
    if(mask64 != (LMASK64(i))) {
      printf("LMASK64(%2d) expected %16.16lx, got %16.16lx\n", i, mask64, LMASK64(i)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("LMASK64Z test : ") ;
  for(i=64, mask64=-1l ; i > 0 ; i-- , mask64 <<= 1){
    if(mask64 != (LMASK64Z(i))) {
      printf("LMASK64Z(%2d) expected %16.16lx, got %16.16lx\n", i, mask64, LMASK64Z(i)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("RMASK64 test : ") ;
  for(i=0, mask64=0l ; i<=64 ; i++, mask64 = (mask64 << 1) | 1l){
    if(mask64 != RMASK64(i)){
      printf("RMASK64(%2d) expected %16.16lx, got %16.16lx\n", i, mask64, RMASK64(i)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("RMASK63 test : ") ;
  for(i=0, mask64=0l ; i<64 ; i++, mask64 = (mask64 << 1) | 1l){
    if(mask64 != RMASK63(i)){
      printf("RMASK63(%2d) expected %16.16lx, got %16.16lx\n", i, mask64, RMASK63(i)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
// test lzcnt/lnzcnt functions
  uint32_t mask32u ;
  uint64_t mask64u ;
  printf("lzcnt_32 test : ") ;
  for(i=0, mask32u=(~0) ; i<33 ; i++, mask32u >>= 1){
    if(i != lzcnt_32(mask32u)){
      printf("lzcnt_32(%8.8x) expected %d, got %d\n", mask32, i, lzcnt_32(mask32u));
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("lnzcnt_32 test : ") ;
  for(i=32, mask32u=(~0) ; i>=0 ; i--, mask32u <<= 1){
    if(i != lnzcnt_32(mask32u)){
      printf("lnzcnt_32(%8.8x) expected %d, got %d\n", mask32, i, lnzcnt_32(mask32u));
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("lzcnt_64 test : ") ;
  for(i=0, mask64u=(~0l) ; i<65 ; i++, mask64u >>= 1){
    if(i != lzcnt_64(mask64u)){
      printf("lzcnt_64(%8.8x) expected %d, got %d\n", mask64, i, lzcnt_32(mask64u));
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("lnzcnt_64 test : ") ;
  for(i=64, mask64u=(~0l) ; i>=0 ; i--, mask64u <<= 1){
    if(i != lnzcnt_64(mask64u)){
      printf("lnzcnt_64(%8.8x) expected %d, got %d\n", mask32, i, lnzcnt_64(mask64u));
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
// tests zigzag
  printf("zigzag test : ") ;
  int32_t lasti ;
  int32_t whatp , whatm ;
  uint32_t zigp, zigm ;
  for(i=-1 ; i<0x7FFFFFFF ; i++){
    whatp = i+1 ;
    whatm = -whatp ;
    zigp = TO_ZIGZAG(whatp) ;
    zigm = TO_ZIGZAG(whatm) ;
    if(FROM_ZIGZAG(zigp) != whatp || FROM_ZIGZAG(zigm) != whatm){
      printf("ZIGZAG expected (%8.8x %8.8x), got (%8.8x %8.8x)\n", whatp, whatm, FROM_ZIGZAG(zigp), FROM_ZIGZAG(zigm)) ;
      e_exit(1) ;
    }
    lasti = whatp ;
  }
  printf("last i = %8.8x, Success\n", lasti) ;
// test negabinaty
  printf("negabinary test : ") ;
  for(i=-1 ; i<0x7FFFFFFF ; i++){
    whatp = i+1 ;
    whatm = -whatp ;
    zigp = TO_NEGABINARY(whatp) ;
    zigm = TO_NEGABINARY(whatm) ;
    if(FROM_NEGABINARY(zigp) != whatp || FROM_NEGABINARY(zigm) != whatm){
      printf("NEGABINARY expected (%8.8x %8.8x), got (%8.8x %8.8x)\n", whatp, whatm, FROM_NEGABINARY(zigp), FROM_NEGABINARY(zigm)) ;
      e_exit(1) ;
    }
    lasti = whatp ;
  }
  printf("last i = %8.8x, Success\n", lasti) ;
// tests bits needed (32 and 64 bits)
  printf("BitsNeeded_u32, BitsNeeded_32 test : ") ;
  for(i=1, whatp=1 ; i<33 ; i++, whatp *= 2){
    whatm = -whatp ;
    if(BitsNeeded_u32(whatp-1) != (i-1) || BitsNeeded_32(whatp-1) != i || BitsNeeded_32(whatm) != i ){
      printf("(%10d) BitsNeeded_u32(%8.8x) = %2d, BitsNeeded_32(%10d|%8.8x) = %2d, BitsNeeded_32(%10d|%8.8x) = %2d\n",
             whatp-1, whatp-1, BitsNeeded_u32(whatp-1), whatp-1, whatp-1, BitsNeeded_32(whatp-1), whatm, whatm, BitsNeeded_32(whatm)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("BitsNeeded_u64, BitsNeeded_64 test : ") ;
  int64_t whatp64, whatm64 ;
  for(i=1, whatp64=1 ; i<65 ; i++, whatp64 += whatp64){
    whatm64 = -whatp64 ;
    if(BitsNeeded_u64(whatp64-1) != (i-1) || BitsNeeded_64(whatp64-1) != i || BitsNeeded_64(whatm64) != i ){
      printf("(%19ld) BitsNeeded_u64(%16.16lx) = %2d, BitsNeeded_64(%19ld|%16.16lx) = %2d, BitsNeeded_64(%20ld|%16.16lx) = %2d\n",
             whatp64-1, whatp64-1, BitsNeeded_u64(whatp64-1), whatp64-1, whatp64-1, BitsNeeded_64(whatp64-1), whatm64, whatm64, BitsNeeded_64(whatm64)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
// test pop counts
  printf("pop count 32 test : ") ;
  for(i=0, whatp=0 ; i<33 ; i++, whatp = (whatp << 1) | 1){
    if(popcnt_32(whatp) != i){
      printf("popcnt_32(%8.8x) : expected %2d, got %2d\n", whatp, i, popcnt_32(whatp)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("pop count 64 test : ") ;
  for(i=0, whatp64=0 ; i<65 ; i++, whatp64 = (whatp64 << 1) | 1){
    if(popcnt_64(whatp64) != i){
      printf("popcnt_64(%8.8x) : expected %2d, got %2d\n", whatp, i, popcnt_64(whatp)) ;
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
}
