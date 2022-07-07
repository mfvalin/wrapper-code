/* Hopefully useful routines for C and FORTRAN
 * Copyright (C) 2020  Recherche en Prevision Numerique
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */
#include <stdint.h>

typedef struct {
  int o;      // quantization offset
  int e;      // largest exponent (min, max, range) (with bias removed)
  int nbits ; // number of useful bits in quantized token
//   float fac ;
} PackHeader;

typedef union{
  int32_t i ;
  float f ;
} FloatInt;

typedef union{
  int64_t l ;
  double d ;
} DoubleLong;

#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )

void float_quantize_prep(float *z, int n, int nbits, PackHeader *p, float maxval, float minval) {
  FloatInt m1, m2, m3;
  DoubleLong m0 ;
  int exp1, exp2, exp3 ;
  float fac32, range;
  double fac64 ;
  int mask_trunc ;
  int offset ;

  if(nbits <= 0 ) return;
  if(nbits > 24 ) nbits = 24 ;               // no more than 24 bits will be kept
  p->nbits = nbits ;

  range = maxval - minval;
  mask_trunc = ( -1 << (24 - nbits) ) ;      // truncation mask
  m1.f = maxval ;
  m2.f = minval ;
  m3.f = range ;
  exp1 = 0xFF & (m1.i >> 23) ;               // exponents for max, min, range
  exp2 = 0xFF & (m2.i >> 23) ;
  exp3 = 0xFF & (m3.i >> 23) ;
  exp1 = (exp2 > exp1) ? exp2 : exp1 ;       // largest exponent ( max, min ) (biased)
  exp1 = (exp3 > exp1) ? exp3 : exp1 ;       // largest exponent ( max, min, range ) (biased)
  p->e = exp1 - 127 ;                        // largest exponent (with bias removed)

  
//   if(m1.i < 1 || m1.i > 254){
//     printf(" m1.i = %d ", m1.i) ;
//     m1.i = 0 ;
//   }
  if(exp1 < 255 && exp1 > 0){                            // O.K. for IEEE float 32
    m1.i = (127 + (23 - (exp1 - 127))) ;       // factor to bring largest exponent to 23
    m1.i <<= 23 ;
    fac32 = m1.f;                            // normalizing multiplier
    offset = minval * fac32 ;                // quantized minimum value
    printf(" fac32 = %10.4g(%d) ", fac32, exp1);
  }else{                                     // must use IEEE float 64
    m0.l = (1023 + (23 - (exp1 - 127)));     // factor to bring largest exponent to 23
    m0.l = m0.l << 52 ;
    fac64 = m0.d;                            // normalizing multiplier
    offset = minval * fac64 ;                // quantized minimum value
    printf(" fac64 = %10.4g(%d) ", fac64, exp1);
  }

  offset = offset & mask_trunc ;    // drop lower (24 - nbits) bits
  p->o = offset ;                   // truncated minimum value
}

void float_quantize_prep_32(float *z, int n, int nbits, PackHeader *p, float maxval, float minval) {
  FloatInt m1, m2, m3;
  int exp1, exp2, exp3 ;
  float fac32, range;
  int mask_trunc ;
  int offset ;

  if(nbits <= 0 ) return;
  if(nbits > 24 ) nbits = 24 ;               // no more than 24 bits will be kept
  p->nbits = nbits ;

  range = maxval - minval;
  mask_trunc = ( -1 << (24 - nbits) ) ;      // truncation mask
  m1.f = maxval ;
  m2.f = minval ;
  m3.f = range ;
  exp1 = 0xFF & (m1.i >> 23) ;               // exponents for max, min, range
  exp2 = 0xFF & (m2.i >> 23) ;
  exp3 = 0xFF & (m3.i >> 23) ;
  exp1 = (exp2 > exp1) ? exp2 : exp1 ;       // largest exponent ( max, min ) (biased)
  exp1 = (exp3 > exp1) ? exp3 : exp1 ;       // largest exponent ( max, min, range ) (biased)
  p->e = exp1 - 127 ;                        // largest exponent (with bias removed)

  m1.i = (127 + (23 - (exp1 - 127))) ;       // factor to bring largest exponent to 23
  if(m1.i < 1 || m1.i > 254){
    printf(" EXP1 = %d ", exp1) ;
  }
  m1.i <<= 23 ;
  fac32 = m1.f;                              // normalizing multiplier

  offset = minval * fac32 ;         // quantized minimum value
  offset = offset & mask_trunc ;    // drop lower (24 - nbits) bits
  p->o = offset ;                   // truncated minimum value
}

void float_quantize_prep_64(float *z, int n, int nbits, PackHeader *p, float maxval, float minval) {
  FloatInt m1, m2, m3;
  DoubleLong m0 ;
  int exp1, exp2, exp3 ;
  float range;
  double fac64 ;
  int mask_trunc ;
  int offset ;

  if(nbits <= 0 ) return;
  if(nbits > 24 ) nbits = 24 ;               // no more than 24 bits will be kept
  p->nbits = nbits ;

  range = maxval - minval;
  mask_trunc = ( -1 << (24 - nbits) ) ;      // truncation mask
  m1.f = maxval ;
  m2.f = minval ;
  m3.f = range ;
  exp1 = 0xFF & (m1.i >> 23) ;               // exponents for max, min, range
  exp2 = 0xFF & (m2.i >> 23) ;
  exp3 = 0xFF & (m3.i >> 23) ;
  exp1 = (exp2 > exp1) ? exp2 : exp1 ;       // largest exponent ( max, min ) (biased)
  exp1 = (exp3 > exp1) ? exp3 : exp1 ;       // largest exponent ( max, min, range ) (biased)
  p->e = exp1 - 127 ;                        // largest exponent (with bias removed)

  m0.l = (1023 + (23 - (exp1 - 127)));       // factor to bring largest exponent to 23
  if(m0.l < 1 || m0.l > 254){
    printf(" EXP1 = %d ", exp1) ;
  }
  m0.l = m0.l << 52 ;
  fac64 = m0.d;                              // normalizing multiplier

  offset = minval * fac64 ;       // quantized minimum value
  offset = offset & mask_trunc ;  // drop lower (24 - nbits) bits
  p->o = offset ;                 // truncated minimum value
}

// quantizer for 32 bit floats producing a stream of unsigned 32 bit integers
// nbits is assumed to be <= 24
void float_quantize(void *iz, float *z, int n, PackHeader *p ) {
  int i, it;
  int offset = p->o ;   // offset reflecting minimum value
  int exp1   = p->e ;   // largest exponent in original float values (bias removed)
  int round = 0 ;
  int *izw = (int   *) iz;
  int nbits = p->nbits ;

  if(nbits <= 0 ) return ;
  if(nbits > 24 ) nbits = 24 ;         // no more than 24 bits will be kept
  if(nbits < 24)                       // no rounding if nbits = 24
    round = 1 << (23 - nbits);         // rounding for quantized value
  round = round - offset;              // combine round and offset

  if(exp1 > 1 && exp1 < 127){
    FloatInt m1;
    float fac32 ;
    m1.i = (127 + (23 - exp1)) << 23 ;   // factor to bring largest exponent to 23
    fac32 = m1.f ;
    for ( i=0 ; i<n ; i++){
      it = z[i] * fac32;                 // "normalize" to largest exponent
      it = (it + round) >> (24-nbits);   // remove offset, add rounding term
      izw[i] = (it < 0) ? 0 : it ;       // quantized result MUST be >= 0
    }
  }else{
    DoubleLong m1;
    double fac64 ;
    m1.l = (1023 + (23 - exp1)) ;        // factor to bring largest exponent to 23
    m1.l = m1.l << 52 ;
    fac64 = m1.d ;
    for ( i=0 ; i<n ; i++){
      it = z[i] * fac64;                 // "normalize" to largest exponent
      it = (it + round) >> (24-nbits);   // remove offset, add rounding term
      izw[i] = (it < 0) ? 0 : it ;       // quantized result MUST be >= 0
    }
  }
}

// quantizer for 32 bit floats producing a stream of unsigned 32 bit integers
// nbits is assumed to be <= 24
void float_quantize_32(void *iz, float *z, int n, int nbits, PackHeader *p ) {
  int i, it;
  FloatInt m1;
  float fac32 ;
  int offset = p->o ;   // offset reflecting minimum value
  int exp1   = p->e ;   // largest exponent in original float values (bias removed)
  int round = 0 ;
  int *izw = (int   *) iz;

  if(nbits <= 0 ) return ;
  if(nbits > 24 ) nbits = 24 ;         // no more than 24 bits will be kept

  m1.i = (127 + (23 - exp1)) << 23;  // factor to bring largest exponent to 23
  // if (127 + (23 - (exp1 - 127))) > 254 there is a problem with the exponent of fac32 (> 127)
  fac32 = ((127 + (23 - exp1)) > 254) ? 0.0f : m1.f;

  if(nbits < 24)                       // no rounding if nbits = 24
    round = 1 << (23 - nbits);         // rounding for quantized value
  round = round - offset;              // combine round and offset

  for ( i=0 ; i<n ; i++){
    it = z[i] * fac32;                   // "normalize" to largest exponent
    it = (it + round) >> (24-nbits);   // remove offset, add rounding term
    izw[i] = (it < 0) ? 0 : it ;       // quantized result MUST be >= 0
  }
}

// quantizer for 32 bit floats producing a stream of unsigned 32 bit integers
// nbits is assumed to be <= 24
void float_quantize_64(void *iz, float *z, int n, int nbits, PackHeader *p ) {
  int i, it;
  DoubleLong m1;
  double fac64 ;
  int offset = p->o ;   // offset reflecting minimum value
  int exp1   = p->e ;   // largest exponent in original float values (bias removed)
  int round = 0 ;
  int *izw = (int   *) iz;

  if(nbits <= 0 ) return ;
  if(nbits > 24 ) nbits = 24 ;         // no more than 24 bits will be kept

  m1.l = (1023 + (23 - exp1)) ;        // factor to bring largest exponent to 23
  m1.l = m1.l << 52 ;
  fac64 = m1.d ;

  if(nbits < 24)                       // no rounding if nbits = 24
    round = 1 << (23 - nbits);         // rounding for quantized value
  round = round - offset;              // combine round and offset

  for ( i=0 ; i<n ; i++){
    it = z[i] * fac64;                 // "normalize" to largest exponent
    it = (it + round) >> (24-nbits);   // remove offset, add rounding term
    izw[i] = (it < 0) ? 0 : it ;       // quantized result MUST be >= 0
  }
}

void float_unquantize(void *iz, float *z, int n, PackHeader *p) {
  int i, t;
  unsigned int   *izw = (unsigned int   *) iz;
  int offset = p->o;
  int exp = p->e;
  int nbits = p->nbits ;

  if(exp > 1 && exp < 127){
    FloatInt m1;
    float fac32;
    m1.i = (127 + exp - 23) ;  // inverse of factor to bring largest exponent to 23
    m1.i <<= 23 ;
    fac32 = m1.f;
    for ( i=0 ; i<n ; i++){
      t = izw[i] << (24-nbits) ;
      z[i] = (t + offset) * fac32;
  }
  }else{
    DoubleLong m1;
    double fac64;
    exp = (exp > 127) ? 127 : exp ;
    m1.l = (exp + 1023 - 23) ;  // inverse of factor to bring largest exponent to 23
    m1.l = m1.l << 52 ;
    fac64 = m1.d;
    for ( i=0 ; i<n ; i++){
      t = izw[i] << (24-nbits) ;
      z[i] = (t + offset) * fac64;
    }
  }
}

void float_unquantize_32(void *iz, float *z, int n, int nbits, PackHeader *p) {
  FloatInt m1;
  float fac32;
  int i, t;
  unsigned int   *izw = (unsigned int   *) iz;
  int offset = p->o;
  int exp = p->e;

  m1.i = (127 + exp - 23) ;  // inverse of factor to bring largest exponent to 23
//   if(m1.i > 254 || m1.i < 1) m1.i = 0 ;
  m1.i <<= 23 ;
  fac32 = m1.f;
  for ( i=0 ; i<n ; i++){
    t = izw[i] << (24-nbits) ;
    z[i] = (t + offset) * fac32;
  }
}

void float_unquantize_64(void *iz, float *z, int n, int nbits, PackHeader *p) {
  DoubleLong m1;
  double fac64;
  int i, t;
  unsigned int   *izw = (unsigned int   *) iz;
  int offset = p->o;
  int exp = p->e;

  m1.l = (exp + 1023 - 23) ;  // inverse of factor to bring largest exponent to 23
  m1.l = m1.l << 52 ;
  fac64 = m1.d;
  for ( i=0 ; i<n ; i++){
    t = izw[i] << (24-nbits) ;
    z[i] = (t + offset) * fac64;
  }
}
#if defined(SELF_TEST)
#define NPTS 32800
#define ABS(a) ((a) > 0 ? a : -(a))
#include <stdio.h>
int main(){
  float zi[NPTS] ;
  float zo[NPTS] ;
  short iz[NPTS] ;
  int iw[NPTS] ;
  PackHeader p;
  int i, j;
  double toler;
  double avg, avgi, avgo ;
  int error=0;
  float minval, maxval ;
  double delta, maxdelta ;
  double avgdelta = 0.0 ;
  int NBITS = 22 ;

  for(i=0 ; i<NPTS ; i++) { zi[i] = .00001 + .012345 * (i - NPTS/2 ) ;  zo[i] = -99999.0 ; }
  avgi = 0.0 ;
  for(i=0 ; i<NPTS ; i++) { zi[i] *= 1.6E+36 ; avgi += zi[i] ;}
  avgi /= NPTS ;
  printf("avgi = %15.8g\n\n", avgi) ;
//   for(i=0 ; i<NPTS ; i++) { zi[i] *= 8.4040326E+35f ; }
  minval = maxval = zi[0] ;
  for(i=1 ; i<NPTS ; i++) { 
    maxval = (zi[i] > maxval) ? zi[i] : maxval ; minval = (zi[i] < minval) ? zi[i] : minval ;
  }

  for(j = 1 ; j < 25 ; j++) {
    NBITS = j ;
    i = 1 << (NBITS);
//     toler = (zi[NPTS-1] - zi[0]);
    toler = zi[NPTS-1]; toler -= zi[0] ;
    toler /= i; //  toler *= 2;

    float_quantize_prep(zi, NPTS, NBITS, &p, maxval, minval );
//     float_quantize_prep_64(zi, NPTS, NBITS, &p, maxval, minval );
//     printf("offset = %10d, exp = %d, nbits = %d, ", p.o, p.e, NBITS);
    if(p.e > 0 && p.e < 127){
      printf(" [32] ");
//       float_quantize_32(iw, zi, NPTS, NBITS, &p ) ;
      float_quantize(iw, zi, NPTS, &p ) ;
//       float_unquantize_32(iw, zo, NPTS, NBITS, &p);
      float_unquantize(iw, zo, NPTS, &p);
    }else{
      printf(" [64] ");
//       float_quantize_64(iw, zi, NPTS, NBITS, &p ) ;
      float_quantize(iw, zi, NPTS, &p ) ;
//       float_unquantize_64(iw, zo, NPTS, NBITS, &p);
      float_unquantize(iw, zo, NPTS, &p);
    }
    printf("nbits = %2d, max, min, rng = %9.5g %9.5g %9.5g, ", NBITS, maxval, minval, maxval-minval);
    avgo = 0.0 ;
    for(i=0 ; i<NPTS ; i++) avgo += zo[i] ;
    avgo /= NPTS ;
    printf("avgo = %15.8g, ", avgo) ;
    error = 0 ; maxdelta = 0.0 ; avgdelta = 0.0 ; avg = 0.0 ;
    for(i=0 ; i<NPTS ; i++) {
//       delta = ABS(zo[i]-zi[i]) ;
      delta = zo[i] ; delta -= zi[i] ;
      avg += delta ;
      delta = ABS(delta) ;
      avgdelta += delta ;
      maxdelta = (delta > maxdelta) ? delta : maxdelta ;
      if(delta > toler) error++;
    }
    avg /= NPTS;
//     printf("%d points from %15.5f to %15.5f, maxdelta = %15.8f, avgdelta = %15.8f\n",
//            NPTS, zi[0], zi[NPTS-1], maxdelta, avgdelta/NPTS);
    printf("maxdelta = %10.4g, avgdelta = %10.4g, ",
           maxdelta, avgdelta/NPTS);
    printf("bias = %15.8g, toler = %9.5g, bias/toler = %8.4f, maxdelta/toler = %6.4f, exceed=%d\n\n",
          avg, toler, avg/toler, maxdelta/toler, error);
  }
return 0 ;
#if 0
  fast_quantize(iz, zi, NPTS, NBITS, &p );
  printf("offset = %10d, exp = %d \n",p.o,p.e);
  fast_unquantize(iz, zo, NPTS, NBITS, &p);
  error = 0 ;
  for(i=0 ; i<NPTS ; i++) {
    if(ABS(zo[i]-zi[i])>toler) error++;
//     printf("%3d %10.5f %10.5f %10.5f %10.5f (%8.5f) %s\n",i,zi[i],zo[i],(zo[i]-zi[i]) / zi[i],(zo[i]-zi[i]) , toler, (ABS(zo[i]-zi[i])>toler)?"*":" ");
    avg += (zo[i]-zi[i]);
  }
  avg /= NPTS;
   printf("from %15.8f to %15.8f\n",zi[0],zi[NPTS-1]);
  printf("bias = %15.8f, toler = %7.5f, bias/toler = %6.4f, toler exceeded=%d\n",avg,toler,avg/toler,error);
#endif
}
#endif
