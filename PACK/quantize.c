#if defined(__SSE4_1__)
#include <x86intrin.h>
static unsigned char shuf[] = {  0,  1,  4,  5,  8,  9, 12, 13,128,128,128,128,128,128,128,128};
#endif

#define VL 4
// nbits is "hidden" through round and mask (it is assumed that nbits <= 16)
void quantize(unsigned short int *iz, float *z, int n, float fac, int round, int mask, int minimum) {
  int i;
#if defined(__SSE4_1__)
  __m128  x0, vfac;
  __m128i i0, vrnd, vmsk, vshf;
#else
  int j, it;
#endif

#if defined(__SSE4_1__)
  // load constants
  vfac = _mm_set1_ps(fac);                // scaling factor
  vrnd = _mm_set1_epi32(round-minimum);   // round - minimum
  vmsk = _mm_set1_epi32(mask);            // clipping value
  vshf = _mm_loadu_si128((__m128i*)shuf); // 8 bit vector to extract 16 bit tokens
#else
  round = round - minimum;
#endif
  for (i=0 ; i<n ; i+=VL){
#if defined(__SSE4_1__)
    x0  = _mm_loadu_ps(&z[i]) ;       // load 4 floats from z
    x0 = _mm_mul_ps(x0,vfac);         // scale
    i0 = _mm_cvtps_epi32(x0);         // convert to integer
    i0 = _mm_add_epi32(i0,vrnd);      // add round - minimum
    i0 = _mm_srli_epi32(i0,8);        // shift right
    i0 = _mm_min_epi32(i0,vmsk);      // min( ... , mask)
    i0 = _mm_shuffle_epi8(i0,vshf);   // move even 16 bit tokens into lower 64 bits
    _mm_storel_epi64((__m128i*)&iz[i],i0);  // store 4 shorts into iz
#else
    for(j=0;j<VL;j++) {
      it = fac * z[j+i];         // scale and convert to integer (my be negative)
      it = it + round;           // add round - minimum (should never be negative)
      it = it >> 8;              // shift right
      iz[j+i] = (it > mask) ? mask : it; // min( ... , mask) and store
    }
#endif
  }
}

// nbits MUST BE <= 16 for this routine to work
void fast_quantize(unsigned short int *iz, float *z, int n, int nbits, float maxval, float minval) {
  union {
    unsigned int i;
    float f;
  } m1, m2, m3;
  int exp1, exp2, exp3, i, it;
  float fac, range;
  int mask_trunc ;
  int mask_nbits ;
  int offset, irange;
  int round ;

  range = maxval - minval;
  round = 1 << (23 - nbits);
  mask_nbits = ~( -1 << nbits );
  mask_trunc = ~( -1 << (24 - nbits) );
  m1.f = maxval;
  m2.f = minval;
  m3.f = range;
  exp1 = 0xFF & (m1.i >> 23) ;
  exp2 = 0xFF & (m2.i >> 23) ;
  exp3 = 0xFF & (m3.i >> 23) ;
  exp1 = (exp1 > exp2) ? exp1 : exp2 ;       // largest exponent ( max, min )
  exp1 = (exp1 > exp3) ? exp1 : exp3 ;       // largest exponent ( max, min, range )
  m1.i = (127 + (23 - (exp1 - 127))) << 23;  // factor to bring largest exponent to 23
  fac = m1.f;

  offset = minval * fac ;
  offset = offset & mask_trunc ;  // drop lower (24 - nbits) bits
  round = round - offset;         // integer offset, (reflects rounded minimum value)
  irange = range * fac ;          // compute right shift count needed to make irange <= mask_nbits

  for ( i=0 ; i<n ; i++){
    it = z[i] * fac;          // "normalize" to largest exponent
    it = (it + round) >> (24-nbits);   // remove offset, add rounding term
    iz[i] = (it > mask_nbits) ? mask_nbits : it; // clamp values at FFFF
  }
}