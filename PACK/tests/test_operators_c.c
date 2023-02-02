#include <stdint.h>
#include <stdio.h>

#include <rmn/misc_operators.h>
#include <rmn/misc_types.h>
#include <rmn/misc_timers.h>

#define NI 167
#define NJ 169
#define NP   33
#define NPF  15
#define NPFI 5
#define NPFJ 3

static int32_t xt2[] = {-8, -7, -7, -6, -6, -5, -5, -4, -4, -3, -3, -2, -2, -1, -1,  0,  0,  0,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8} ;
static int32_t xt4[] = {-4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4} ;
static int32_t xt8[] = {-2, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  2} ;
static int32_t xr2[] = {-8, -8, -7, -7, -6, -6, -5, -5, -4, -4, -3, -3, -2, -2, -1, -1,  0,  1,  1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8} ;
static int32_t xr4[] = {-4, -4, -4, -3, -3, -3, -3, -2, -2, -2, -2, -1, -1, -1, -1,  0,  0,  0,  1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4} ;
static int32_t xr8[] = {-2, -2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2} ;
static int32_t xi[]  = { 999, 0, 1, 2, 3, 4, 5, 6, 7, 999 } ;

static float xs1[NPF+64] ;
static float xs2[NPF+64] ;

static void init_floats(){
  int i ;
  FloatInt f ;
  float start = 1.0f;

  for(i=0 ; i<NPF/2 ; i++){
    f.i = (127 << 23) | (0x7FFFFF >> 3*i) | 0x600000 ;
    xs1[i] = f.f ;
    xs1[NPF-1-i] = -f.f ;
    if(i == NPF/2 - 1){
      f.f = start ;
      f.i |= 0x7FFFFF ;
      start = f.f ;
    }
    xs2[i] = start ;
    xs2[NPF-1-i] = -start ;
    start *= 1.321f ;
  }
  xs1[NPF/2] = 0.0f ;
  xs2[NPF/2] = 0.0f ;
//   for(i=0 ; i<NPF/2 ; i++){ printf("%11.8f", xs1[i]) ; } ;  printf("\n") ;
  for(i=0 ; i<NPF ; i++){ printf("%11.7f", xs2[i]) ; } ;  printf("\n") ;
}

static int compare_int(int32_t *a, int32_t *b, int n){
  int i, e = 0 ;
  for(i=0 ; i<n ; i++) e += ( (a[i] == b[i]) ? 0 : 1 ) ;
  return e ;
}

static void e_exit(int n){
  printf("ERROR in test, aborting\n") ;
  exit(n) ;
}

static void fill_64(float *src, int n, float f64[64]){
  int i, j ;
  for(i=0, j=0 ; i<64 ; i++, j=(j+1 >= n) ? 0 : j+1){
    f64[i] = src[j] ;
  }
}

int main(int argc, char **argv){
  int32_t src[NP], dst[NP+8] ;
  int32_t i, j, i0, j0, lni, lnj, x, err, indx ;
  float f ;
  float xf[NP], xf1[NPF+64] ;
  FloatInt fi ;
  int32_t memmask[8] ;
  int32_t blk_test[NJ][NI] ;
  int32_t blk_new[NJ][NI] ;
  int32_t blk[64] ;
  uint64_t t1, t2 ;
#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  __m128i v128 ;
  __m256i v256 ;
  __m128i v128lo, v128hi ;
#endif
  ieee_prop prop ;
  uint16_t stream16[65] ;

#if defined(__x86_64__) && defined(__AVX2__) && defined(WITH_SIMD)
  printf("vector masks (128 and 256) : ") ;
  printf("\n") ;
  for(i=0 ; i<5 ; i++) {
    v128 = _mm_memmask_si128(i) ;
    _mm_storeu_si128((__m128i *)memmask, v128 ) ;
//     printf("mask128(%d) = ", i) ;
//     for(j=0 ; j<4 ; j++) printf("%8.8x ", memmask[j]) ;  printf("\n") ;
    err = 0 ;
    for(j=0 ; j<i ; j++) if(memmask[j] != -1) err++ ;
    for(    ; j<4 ; j++) if(memmask[j] !=  0) err++ ;
    if(err > 0) e_exit(1) ;
  }
  printf("Success (128)\n") ;
  for(i=0 ; i<9 ; i++) {
    v256 = _mm256_memmask_si256(i) ;
    _mm256_storeu_si256((__m256i *)memmask, v256 ) ;
//     printf("mask256(%d) = ", i) ;
//     for(j=0 ; j<8 ; j++) printf("%8.8x ", memmask[j]) ;  printf("\n") ;
    for(j=0 ; j<i ; j++) if(memmask[j] != -1) err++ ;
    for(    ; j<8 ; j++) if(memmask[j] != 0) err++ ;
    if(err > 0) e_exit(1) ;
  }
  if(err > 0) e_exit(1) ;
  printf("Success (256)\n") ;

#endif
  printf("extract blocks : ") ;
  for(j=0 ; j<NJ ; j++){
    for(i=0 ; i<NI ; i++){
      blk_test[j][i] = (i << 8) + (j) ;
    }
  }
  for(j0=0 ; j0<NJ ; j0+=8){
    lnj = (NJ-j0 < 8) ? NJ-j0 : 8 ;
    for(i0=0 ; i0<NI ; i0+=8){
      lni = (NI-i0 < 8) ? NI-i0 : 8 ;
      indx = i0 + j0 * NI ;
      get_w32_block((void *)(&blk_test[j0][i0]), blk, lni, NI, lnj) ;
      put_w32_block((void *)(&blk_new[j0][i0]) , blk, lni, NI, lnj) ;
      err = 0 ;
      for(j=lnj-1 ; j>=0 ; j--){
        for(i=0 ; i<lni ; i++){
          if(blk[j*lni+i] != blk_test[j0+j][i0+i]) err++ ;
          if(blk_new[j0+j][i0+i] != blk_test[j0+j][i0+i]) err++ ;
        }
      }
      if(err > 0) {
        printf("\n") ;
        printf("i0 = %d, j0 = %d, indx = %d, lni = %d, lnj = %d\n", i0, j0, indx, lni, lnj) ;
        for(j=lnj-1 ; j>=0 ; j--){
          for(i=0 ; i<lni ; i++){
          printf("%8.8x|%8.8x ", blk_test[j0+j][i0+i], blk[j*lni+i] ) ;
          }
          printf("\n") ;
        }
        e_exit(1) ;
      }
    }
  }
  printf("Success\n") ;
  t1 = elapsed_cycles() ;
  for(j0=0 ; j0<NJ ; j0+=8){
    lnj = (NJ-j0 < 8) ? NJ-j0 : 8 ;
    for(i0=0 ; i0<NI ; i0+=8){
      lni = (NI-i0 < 8) ? NI-i0 : 8 ;
      indx = i0 + j0 * NI ;
      get_w32_block((void *)(&blk_test[j0][i0]), blk, lni, NI, lnj) ;
    }
  }
  t2 = elapsed_cycles() - t1 ;
  double t = t2;
  t /= (NI*NJ) ;
  printf("extract %d words in %6.2f cycles/word\n", NI*NJ, t);
  t1 = elapsed_cycles() ;
  for(j0=0 ; j0<NJ ; j0+=8){
    lnj = (NJ-j0 < 8) ? NJ-j0 : 8 ;
    for(i0=0 ; i0<NI ; i0+=8){
      lni = (NI-i0 < 8) ? NI-i0 : 8 ;
      indx = i0 + j0 * NI ;
      put_w32_block((void *)(&blk_new[j0][i0]), blk, lni, NI, lnj) ;
    }
  }
  t2 = elapsed_cycles() - t1 ;
  t = t2;
  t /= (NI*NJ) ;
  printf("insert %d words in %6.2f cycles/word\n", NI*NJ, t);

return 0;
  init_floats() ;
  printf("\n") ;
  prop = ieee_properties(xs1, NPF/2) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.7f", xs1[i]) ; printf("\n") ;
  printf("emax = %d, emin = %2.2x, mima = %d %s%s%s, n = %d\n", 
         prop.emax, prop.emin, prop.mima, prop.allp ? ", all >= 0" : "",  prop.allm ? ", all < 0" : "", prop.zero ? ", zero" : "", NPF/2 );
  prop = ieee_encode_block_16(xs1, NPF/2, 1, stream16) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.4x", stream16[i+1]) ; printf("\n") ;
  prop = ieee_decode_block_16(xf1, NPF/2, 1, stream16) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("=======================================================================================\n") ;
  prop = ieee_properties(xs1, 65) ;
  printf("expected error : prop.errf = %x\n", prop.errf) ;
  printf("=======================================================================================\n") ;
  fill_64(xs1, NPF/2, xf1) ;
  prop = ieee_properties_64(xf1) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("emax = %d, emin = %2.2x, mima = %d %s%s%s, n = %d\n", 
         prop.emax, prop.emin, prop.mima, prop.allp ? ", all >= 0" : "",  prop.allm ? ", all < 0" : "", prop.zero ? ", zero" : "", 64 );
  prop = ieee_encode_block_16(xf1, 8, 8, stream16) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.4x", stream16[i+1]) ; printf("\n") ;
  prop = ieee_decode_block_16(xf1, 8, 8, stream16) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("=======================================================================================\n") ;

  prop = ieee_properties(xs1, NPF/2+1) ;
  for(i=0 ; i<NPF/2+1 ; i++) printf("%11.7f", xs1[i]) ; printf("\n") ;
  printf("emax = %d, emin = %2.2x, mima = %d %s%s%s, n = %d\n", 
         prop.emax, prop.emin, prop.mima, prop.allp ? ", all >= 0" : "",  prop.allm ? ", all < 0" : "", prop.zero ? ", zero" : "", NPF/2+1 );
  prop = ieee_encode_block_16(xs1, NPF/2+1, 1, stream16) ;
  for(i=0 ; i<NPF/2+1 ; i++) printf("%11.4x", stream16[i+1]) ; printf("\n") ;
  prop = ieee_decode_block_16(xf1, NPF/2+1, 1, stream16) ;
  for(i=0 ; i<NPF/2+1 ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("=======================================================================================\n") ;

  fill_64(xs1, NPF/2+1, xf1) ;
  prop = ieee_properties_64(xf1) ;
  for(i=0 ; i<NPF/2+1 ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("emax = %d, emin = %2.2x, mima = %d %s%s%s, n = %d\n", 
         prop.emax, prop.emin, prop.mima, prop.allp ? ", all >= 0" : "",  prop.allm ? ", all < 0" : "", prop.zero ? ", zero" : "", 64 );
  prop = ieee_encode_block_16(xf1, 8, 8, stream16) ;
  for(i=0 ; i<NPF/2+1 ; i++) printf("%11.4x", stream16[i+1]) ; printf("\n") ;
  prop = ieee_decode_block_16(xf1, 8, 8, stream16) ;
  for(i=0 ; i<NPF/2+1 ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("=======================================================================================\n") ;

  prop = ieee_properties(xs1, NPF) ;
  for(i=0 ; i<NPF ; i++) printf("%11.7f", xs1[i]) ; printf("\n") ;
  printf("emax = %d, emin = %2.2x, mima = %d %s%s%s, n = %d\n", 
         prop.emax, prop.emin, prop.mima, prop.allp ? ", all >= 0" : "",  prop.allm ? ", all < 0" : "", prop.zero ? ", zero" : "", NPF );
  prop = ieee_encode_block_16(xs1, NPFI, NPFJ, stream16) ;
  for(i=0 ; i<NPF ; i++) printf("%11.4x", stream16[i+1]) ; printf("\n") ;
  prop = ieee_decode_block_16(xf1, NPFI, NPFJ, stream16) ;
  for(i=0 ; i<NPF ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("=======================================================================================\n") ;

  prop = ieee_properties(xs1+NPF/2+1, NPF/2) ;
  for(i=NPF/2+1 ; i<NPF ; i++) printf("%11.7f", xs1[i]) ; printf("\n") ;
  printf("emax = %d, emin = %2.2x, mima = %d %s%s%s, n = %d\n", 
         prop.emax, prop.emin, prop.mima, prop.allp ? ", all >= 0" : "",  prop.allm ? ", all < 0" : "", prop.zero ? ", zero" : "", NPF/2  );
  prop = ieee_encode_block_16(xs1+NPF/2+1, NPF/2, 1, stream16) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.4x", stream16[i+1]) ; printf("\n") ;
  prop = ieee_decode_block_16(xf1, NPF/2, 1, stream16) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("=======================================================================================\n") ;

  fill_64(xs1+NPF/2+1, NPF/2, xf1) ;
  prop = ieee_properties_64(xf1) ;
  for(i=NPF/2+1 ; i<NPF ; i++) printf("%11.7f", xs1[i]) ; printf("\n") ;
  printf("emax = %d, emin = %2.2x, mima = %d %s%s%s, n = %d\n", 
         prop.emax, prop.emin, prop.mima, prop.allp ? ", all >= 0" : "",  prop.allm ? ", all < 0" : "", prop.zero ? ", zero" : "", 64 );
  prop = ieee_encode_block_16(xf1, 8, 8, stream16) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.4x", stream16[i+1]) ; printf("\n") ;
  prop = ieee_decode_block_16(xf1, 8, 8, stream16) ;
  for(i=0 ; i<NPF/2 ; i++) printf("%11.7f", xf1[i]) ; printf("\n") ;
  printf("=======================================================================================\n") ;
  prop = ieee_encode_block_16(xf1, 7, 9, stream16) ;
  printf("expected error : prop.errf = %x\n", prop.errf) ;
  printf("=======================================================================================\n") ;
  prop = ieee_decode_block_16(xf1, 16, 4, stream16) ;
  printf("expected error : prop.errf = %x\n", prop.errf) ;
  printf("=======================================================================================\n") ;
return 0 ;
// IEEE min and max exponent test
  printf("IEEE properties : ") ;
  for(i=0 ; i<NP ; i++) xf[i] = (i - NP/2) * 1.01f ;
  uint32_t emax = ieee_max_exponent(xf, NP) ;
  uint32_t emin = ieee_min_exponent(xf, NP) ;
  printf("x = %f -> %f -> %f -> %f -> %f, emin = %d, emax = %d : ", 
         xf[0], xf[NP/2-1], xf[NP/2], xf[NP/2+1], xf[NP-1], emin, emax) ;
  if(emin != 127 || emax != 131) {
    printf("expecting emax = 131, got %d, emin = 127, got %d\n",
           emax, emin) ;
    e_exit(1) ;
  }
  prop = ieee_properties(xf, NP);   // test full array
  printf("Success\n") ;
  printf("properties full array  (size = %ld bits), emin = %d, emax = %d, allp = %d, allm = %d : ", 8*sizeof(prop), prop.emin, prop.emax, prop.allp, prop.allm) ;
  if(prop.emin != 127 || prop.emax != 131 || prop.allp != 0 || prop.allm != 0) {
    printf("\nexpected(got) emin = 127(%d), emax = 131(%d), allp = 1(%d), allm = 0(%d) : ", prop.emin, prop.emax, prop.allp, prop.allm) ;
    e_exit(1) ;
  }
  printf("Success\n") ;
  prop = ieee_properties(xf, NP/2);           // test bottom of array (all negative)
  printf("properties bottom half (size = %ld bits), emin = %d, emax = %d, allp = %d, allm = %d : ", 8*sizeof(prop), prop.emin, prop.emax, prop.allp, prop.allm) ;
  if(prop.emin != 127 || prop.emax != 131 || prop.allp != 0 || prop.allm != 1) {
    printf("\nexpected(got) emin = 127(%d), emax = 131(%d), allp = 1(%d), allm = 0(%d) : ", prop.emin, prop.emax, prop.allp, prop.allm) ;
    e_exit(1) ;
  }
  printf("Success\n") ;
  prop = ieee_properties(xf+NP/2, NP/2+1);    // test top of array (all non negative)
  printf("properties top half    (size = %ld bits), emin = %d, emax = %d, allp = %d, allm = %d : ", 8*sizeof(prop), prop.emin, prop.emax, prop.allp, prop.allm) ;
  if(prop.emin != 127 || prop.emax != 131 || prop.allp != 1 || prop.allm != 0) {
    printf("\nexpected(got) emin = 127(%d), emax = 131(%d), allp = 1(%d), allm = 0(%d) : ", prop.emin, prop.emax, prop.allp, prop.allm) ;
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

// extract 128 bits from 256 bit register (X86_64)
#if defined(__x86_64__) && defined(__AVX2__)
  printf("AVX2 SIMD extract upper/lower 128 bits : ");
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
  printf("IDIV2T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt2, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV2T_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-2TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt2, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV2T_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-2TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt2, NP)) ; if(err) e_exit(1) ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV2R(src[i]) ;
  printf("IDIV2R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr2, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV2R_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-2RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr2, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV2R_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-2RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr2, NP)) ; if(err) e_exit(1) ;
  printf("\n") ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV4T(src[i]) ;
  printf("IDIV4T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", IDIV4T(src[i])) ; printf(" err = %d\n", err = compare_int(dst, xt4, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV4T_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-4TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt4, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV4T_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-4TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt4, NP)) ; if(err) e_exit(1) ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV4R(src[i]) ;
  printf("IDIV4R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", IDIV4R(src[i])) ; printf(" err = %d\n", err = compare_int(dst, xr4, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV4R_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-4RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr4, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), _mm256_idiv4r_epi32(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-4RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr4, NP)) ; if(err) e_exit(1) ;
  printf("\n") ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV8T(src[i]) ;
  printf("IDIV8T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt8, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV8T_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-8TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt8, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV8T_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-8TV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt8, NP)) ; if(err) e_exit(1) ;

  for(i = 0 ; i < NP ; i++) dst[i] = IDIV8R(src[i]) ;
  printf("IDIV8R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr8, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=4) _mm_storeu_si128((__m128i *)(dst +  i), IDIV8R_128(_mm_loadu_si128((__m128i *)(src + i)))) ;
  printf("128-8RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr8, NP)) ; if(err) e_exit(1) ;
  for(i=0 ; i<NP ; i+=8) _mm256_storeu_si256((__m256i *)(dst + i), IDIV8R_256(_mm256_loadu_si256((__m256i *)(src + i)))) ;
  printf("256-8RV") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr8, NP)) ; if(err) e_exit(1) ;
  printf("\n") ;

#else
  printf("/2     ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", src[i]/2) ; printf("\n") ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV2T(src[i]) ;
  printf("IDIV2T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt2, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV2R(src[i]) ;
  printf("IDIV2R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr2, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV4T(src[i]) ;
  printf("IDIV4T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", IDIV4T(src[i])) ; printf(" err = %d\n", err = compare_int(dst, xt4, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV4R(src[i]) ;
  printf("IDIV4R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", IDIV4R(src[i])) ; printf(" err = %d\n", err = compare_int(dst, xr4, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV8T(src[i]) ;
  printf("IDIV8T ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xt8, NP)) ; if(err) e_exit(1) ;
  for(i = 0 ; i < NP ; i++) dst[i] = IDIV8R(src[i]) ;
  printf("IDIV8R ") ; for(i = 0 ; i < NP ; i++) printf("%3d,", dst[i]) ; printf(" err = %d\n", err = compare_int(dst, xr8, NP)) ; if(err) e_exit(1) ;
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
      printf("lzcnt_64(%16.16lx) expected %d, got %d\n", mask64, i, lzcnt_32(mask64u));
      e_exit(1) ;
    }
  }
  printf("Success\n") ;
  printf("lnzcnt_64 test : ") ;
  for(i=64, mask64u=(~0l) ; i>=0 ; i--, mask64u <<= 1){
    if(i != lnzcnt_64(mask64u)){
      printf("lnzcnt_64(%16.16lx) expected %d, got %d\n", mask64u, i, lnzcnt_64(mask64u));
      e_exit(1) ;
    }
  }
  printf("Success\n") ;\

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
