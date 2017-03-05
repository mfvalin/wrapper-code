
typedef struct {
  int o;    // offset
  int e;    // exponent
} PackHeader;

typedef union{
    unsigned int i;
    float f;
} FltInt;

#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )

#if defined(ALL)
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
// if 8 < nbits <= 16, unsigned shorts will be stored
// if 0 < nbits <= 8, unsigned bytes will be stored
#endif

void fast_quantize(void *iz, float *z, int n, int nbits, PackHeader *p ) {
  FltInt m1, m2, m3;
  int exp1, exp2, exp3, i, j, it;
  float fac, range;
  int mask_trunc ;
  int mask_nbits ;
  int offset, irange;
  int round = 0 ;
  unsigned char *izb = (unsigned char *) iz;
  unsigned short *izs = (unsigned short *) iz;
  float mi[8], ma[8], su[8];
  float maxval, minval;

  if(nbits <= 0 ) return;

  for(i=0; i<8 ; i++) mi[i] = z[i];
  for(i=0; i<8 ; i++) ma[i] = z[i];
  for(i=0; i<8 ; i++) su[i] = z[i];

  for(j=(n&7) ; j<n ; j+=8){
    for(i=0; i<8 ; i++){
      ma[i] = MAX(ma[i] , z[j+i]);
      mi[i] = MIN(mi[i] , z[j+i]);
      su[i] = su[i] + z[i];
    }
  }
  for(i=1; i<8 ; i++) {
    ma[0] = MAX(ma[0] , ma[i]);
    mi[0] = MIN(mi[0] , mi[i]);
    su[0] += su[i];
  }
  maxval = ma[0];
  minval = mi[0];

  range = maxval - minval;
  round = 1 << (23 - nbits);
  mask_nbits = ~( -1 << nbits );
  mask_trunc = ( -1 << (24 - nbits) );
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

  if(nbits > 8) {               // produce a stream of unsigned shorts
    for ( i=0 ; i<n ; i++){
      it = z[i] * fac;          // "normalize" to largest exponent
      it = (it + round) >> (24-nbits);   // remove offset, add rounding term
      izs[i] = (it > mask_nbits) ? mask_nbits : it; // clamp values at FFFF
    }
  }else{                        // produce a stream of unsigned bytes
    for ( i=0 ; i<n ; i++){
      it = z[i] * fac;          // "normalize" to largest exponent
      it = (it + round) >> (24-nbits);   // remove offset, add rounding term
      izb[i] = (it > mask_nbits) ? mask_nbits : it; // clamp values at FFFF
    }
  }
  p->o = offset;
  p->e = exp1;
}
void fast_unquantize(void *iz, float *z, int n, int nbits, PackHeader *p) {
  FltInt m1;
  float fac;
  int i, t;
  unsigned char *izb = (unsigned char *) iz;
  unsigned short *izs = (unsigned short *) iz;
  int offset = p->o;
  int exp = p->e;

  m1.i = (exp - 23) << 23;  // inverse of factor to bring largest exponent to 23
  fac = m1.f;
  if(nbits > 8) {
    for ( i=0 ; i<n ; i++){
      t = izs[i] << (24-nbits) ;
      z[i] = (t + offset) * fac;
    }
  }else{
    for ( i=0 ; i<n ; i++){
      t = izb[i] << (24-nbits) ;
      z[i] = (t + offset) * fac;
    }
  }

}
#if defined(SELF_TEST)
#define NTEST 32800
#define NBITS 14
#define ABS(a) ((a) > 0 ? a : -(a))
#include <stdio.h>
main(){
  float zi[NTEST] ;
  float zo[NTEST] ;
  short iz[NTEST] ;
  PackHeader p;
  int i;
  float toler;
  double avg = 0.0;
  int error=0;

  for(i=0 ; i<NTEST ; i++) { zi[i] = .00001 + .012345 * (i - NTEST/2 ) ;  zo[i] = -99999.0 ; }
  toler = (zi[NTEST-1] - zi[0]);
  i = 1 << (NBITS);
  toler /= i; //  toler *= 2;
  fast_quantize(iz, zi, NTEST, NBITS, &p );
  printf("offset = %10d \n",p.o);
  fast_unquantize(iz, zo, NTEST, NBITS, &p);
  for(i=0 ; i<NTEST ; i++) {
    if(ABS(zo[i]-zi[i])>toler) error++;
//     printf("%3d %10.5f %10.5f %10.5f %10.5f (%8.5f) %s\n",i,zi[i],zo[i],(zo[i]-zi[i]) / zi[i],(zo[i]-zi[i]) , toler, (ABS(zo[i]-zi[i])>toler)?"*":" ");
    avg += (zo[i]-zi[i]);
  }
  avg /= NTEST;
   printf("from %15.8f to %15.8f\n",zi[0],zi[NTEST-1]);
  printf("bias = %15.8f, toler = %7.5f, bias/toler = %6.4f, toler exceeded=%d\n",avg,toler,avg/toler,error);
}
#endif
