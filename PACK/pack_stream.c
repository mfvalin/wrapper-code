/*
 * Copyright (C) 2021  Environnement Canada
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * Author:
 *     M. Valin,   Recherche en Prevision Numerique, 2021
 */
#include <stdio.h>
#include <stdint.h>

#include <pack_stream.h>

// extract n unsigned tokens from bit stream src, store into 64 bit unsigned integers
int pstream_unpack_u64(pstream *ps, uint32_t *src, uint64_t *dest, uint32_t n) {
  uint64_t min64 ;
  int nbt, nbm, i, npt, sgn ;
  uint64_t t64 ;

  pstream_init(ps, src) ;
  sgn = pstream_get_32(ps, 16) ;
  if(sgn != 0xABBA)  return 0  ;       // invalid header
  nbt = pstream_get_32(ps,  8) ;
  nbm = pstream_get_32(ps,  8) ;
  npt = pstream_get_32(ps, 32) ;
  if(npt < n) return 0 ;               // wrong data count
  min64 = pstream_get_64(ps, nbm) ;

  for(i=0 ; i<n ; i++){
    t64 = pstream_get_64(ps, nbt) ;
    dest[i] = t64 + min64 ;
  }
  return n ;    // number of tokens extracted, should be equal to n
}

// extract n signed tokens from bit stream src, store into 64 bit signed integers
int pstream_unpack_i64(pstream *ps, uint32_t *src, int64_t *dest, uint32_t n) {
  int64_t min64 ;
  int nbt, nbm, i, npt, sgn ;
  uint64_t t64 ;

  pstream_init(ps, src) ;
  sgn = pstream_get_32(ps, 16) ;
  if(sgn != 0xABBA)  return 0  ;       // invalid header
  nbt = pstream_get_32(ps,  8) ;
  nbm = pstream_get_32(ps,  8) ;
  npt = pstream_get_32(ps, 32) ;
  if(npt < n) return 0 ;               // wrong data count
  t64 = pstream_get_64(ps, nbm) ;
  min64 = t64 ;
  if(nbm == 32){                       // propagate sign
    min64 <<= 32 ;
    min64 >>= 32 ;
  }

  for(i=0 ; i<n ; i++){
    t64 = pstream_get_64(ps, nbt) ;
    dest[i] = t64 + min64 ;
  }
  return n ;    // number of tokens extracted, should be equal to n
}

// extract n unsigned tokens from bit stream src, store into 32 bit unsigned integers
int pstream_unpack_u32(pstream *ps, uint32_t *src, uint32_t *dest, uint32_t n) {
  uint32_t min32 ;
  int nbt, nbm, i, npt, sgn ;
  uint32_t t32 ;

  pstream_init(ps, src) ;
  sgn = pstream_get_32(ps, 16) ;
  if(sgn != 0xABBA)  return 0  ;       // invalid header
  nbt = pstream_get_32(ps,  8) ;
  if(nbt > 32) return 0 ;              // this is not a 32 bit stream
  nbm = pstream_get_32(ps,  8) ;
  if(nbm > 32) return 0 ;              // this is not a 32 bit stream
  npt = pstream_get_32(ps, 32) ;
  if(npt < n) return 0 ;              // wrong data count
  min32 = pstream_get_32(ps, nbm) ;

  for(i=0 ; i<n ; i++){
    t32 = pstream_get_32(ps, nbt) ;
    dest[i] = t32 + min32 ;
  }
  return n ;
}

// extract n signed tokens from bit stream src, store into 32 bit signed integers
int pstream_unpack_i32(pstream *ps, uint32_t *src, int32_t *dest, uint32_t n) {
  int32_t min32 ;
  int nbt, nbm, i, npt, sgn ;
  uint32_t t32 ;

  pstream_init(ps, src) ;
  sgn = pstream_get_32(ps, 16) ;
  if(sgn != 0xABBA)  return 0  ;       // invalid header
  nbt = pstream_get_32(ps,  8) ;
  if(nbt > 32) return 0 ;              // this is not a 32 bit stream
  nbm = pstream_get_32(ps,  8) ;
  if(nbm > 32) return 0 ;              // this is not a 32 bit stream
  npt = pstream_get_32(ps, 32) ;
  if(npt < n) return 0 ;              // wrong data count
  t32 = pstream_get_32(ps, nbm) ;
  min32 = (int32_t) t32 ;

  for(i=0 ; i<n ; i++){
    t32 = pstream_get_32(ps, nbt) ;
    dest[i] = t32 + min32 ;
  }
  return n ;
}

int pstream_pack_u64(pstream *ps, uint64_t *src, uint32_t *dest, uint32_t n, int bavail) {
  uint64_t min64 = *src ;
  uint64_t max64 = *src ;
  uint64_t range = 0 ;
  int32_t  i, nbt ;
  uint32_t token, mask ;
  uint64_t tok64, mask64 ;
  uint64_t bneeded ;

  for(i = 1 ; i < n ; i++) {   // find minimum and maximum values
    min64 = (src[i] < min64) ? src[i] : min64 ;
    max64 = (src[i] > max64) ? src[i] : max64 ;
  }     // find extrema
  range = max64 - min64 ;
  nbt = 0 ;
  while(range > 0) {           // compute number of bits necessary to encode the data range
    range >>= 1 ;
    nbt ++ ;
  }
  range = max64 - min64 ;
  mask = 0 ; mask = ~mask ; mask = mask >> (32 - nbt) ;  // nbt bits set to 1 on the right

  bneeded = nbt ;
  bneeded = bneeded * n ;               // bits needed = nbt * number of points
  bneeded += (16 + 8 + 8 + 32 + 64) ;   // add header size
  bneeded = (bneeded + 31) / 32 ;       // in uint32_t units
  if(bneeded > bavail) return -1 ;      // ERROR: not enough space for bit stream

  pstream_init(ps, dest) ;              // initialize stream packing structure
  pstream_put_32(ps, 0xABBA, 16) ;      // signature
  pstream_put_32(ps, nbt, 8) ;          // number of bits per token
  pstream_put_32(ps, 64, 8) ;           // number of bits for minimum value is 64
  pstream_put_32(ps, n, 32) ;           // number of elements
  pstream_put_64(ps, min64, 64) ;       // minimum value

  if(nbt <= 32) {   // elements 32 bits or less
    for(i = 0 ; i < n ; i++){
      token = src[i] - min64 ; token &= mask ;
      pstream_put_32(ps, token, nbt) ;
    }
  }else{            // elements wider thatn 32 bits
    mask64 = 0 ; mask64 = ~mask64 ; mask64 = mask64 >> (64 - nbt) ;
    for(i = 0 ; i < n ; i++){
      tok64 = src[i] - min64 ; tok64 &= mask64 ;
      pstream_put_64(ps, tok64, nbt) ;
    }
  }

  pstream_flush(ps) ;
  return (ps->p - ps->s) ;
}

int pstream_pack_i64(pstream *ps, int64_t *src, uint32_t *dest, uint32_t n, int bavail) {
  int64_t min64 = *src ;
  int64_t max64 = *src ;
  uint64_t range = 0 ;
  int32_t  i, nbt ;
  uint32_t token, mask ;
  uint64_t tok64, mask64 ;
  uint64_t bneeded ;

  for(i = 1 ; i < n ; i++) {   // find minimum and maximum values
    min64 = (src[i] < min64) ? src[i] : min64 ;
    max64 = (src[i] > max64) ? src[i] : max64 ;
  }     // find extrema
  range = max64 - min64 ;
  nbt = 0 ;
  while(range > 0) {           // compute number of bits necessary to encode the data range
    range >>= 1 ;
    nbt ++ ;
  }
  range = max64 - min64 ;
  mask = 0 ; mask = ~mask ; mask = mask >> (32 - nbt) ;  // nbt bits set to 1 on the right

  bneeded = nbt ;
  bneeded = bneeded * n ;               // bits needed = nbt * number of points
  bneeded += (16 + 8 + 8 + 32 + 64) ;   // add header size
  bneeded = (bneeded + 31) / 32 ;       // in uint32_t units
  if(bneeded > bavail) return -1 ;      // ERROR: not enough space for bit stream

  pstream_init(ps, dest) ;              // initialize stream packing structure
  pstream_put_32(ps, 0xABBA, 16) ;      // signature
  pstream_put_32(ps, nbt, 8) ;          // number of bits per token
  pstream_put_32(ps, 64, 8) ;           // number of bits for minimum value is 64
  pstream_put_32(ps, n, 32) ;           // number of elements
  pstream_put_64(ps, min64, 64) ;       // minimum value

  if(nbt <= 32) {   // 32 bits or less
    for(i = 0 ; i < n ; i++){
      token = src[i] - min64 ; token &= mask ;
      pstream_put_32(ps, token, nbt) ;
    }
  }else{
    mask64 = 0 ; mask64 = ~mask64 ; mask64 = mask64 >> (64 - nbt) ;
    for(i = 0 ; i < n ; i++){
      tok64 = src[i] - min64 ; tok64 &= mask64 ;
      pstream_put_64(ps, tok64, nbt) ;
    }
  }

  pstream_flush(ps) ;
  return (ps->p - ps->s) ;
}

int pstream_pack_u32(pstream *ps, uint32_t *src, uint32_t *dest, uint32_t n, int bavail) {
  uint32_t min32 = *src ;
  uint32_t max32 = *src ;
  uint32_t range = 0 ;
  int32_t  i, nbt ;
  uint32_t token, mask ;
  uint64_t bneeded ;

  for(i = 1 ; i < n ; i++) {   // find minimum and maximum values
    min32 = (src[i] < min32) ? src[i] : min32 ;
    max32 = (src[i] > max32) ? src[i] : max32 ;
  }     // find extrema
  range = max32 - min32 ;
  nbt = 0 ;
  while(range > 0) {           // compute number of bits necessary to encode the data range
    range >>= 1 ;
    nbt ++ ;
  }
  range = max32 - min32 ;
  mask = 0 ; mask = ~mask ; mask = mask >> (32 - nbt) ;  // nbt bits set to 1 on the right

  bneeded = nbt ;
  bneeded = bneeded * n ;               // bits needed = nbt * number of points
  bneeded += (16 + 8 + 8 + 32 + 64) ;   // add header size
  bneeded = (bneeded + 31) / 32 ;       // in uint32_t units
  if(bneeded > bavail) return -1 ;      // ERROR: not enough space for bit stream

  pstream_init(ps, dest) ;              // initialize stream packing structure
  pstream_put_32(ps, 0xABBA, 16) ;      // signature
  pstream_put_32(ps, nbt, 8) ;          // number of bits per token
  pstream_put_32(ps, 32, 8) ;           // number of bits for minimum value is 32
  pstream_put_32(ps, n, 32) ;           // number of elements
  pstream_put_32(ps, min32, 32) ;       // minimum value

  for(i = 0 ; i < n ; i++){
    token = src[i] - min32 ; token &= mask ;
    pstream_put_32(ps, token, nbt) ;
  }

  pstream_flush(ps) ;
  return (ps->p - ps->s) ;
}

int pstream_pack_i32(pstream *ps, int32_t *src, uint32_t *dest, uint32_t n, int bavail) {
  int32_t min32 = *src ;
  int32_t max32 = *src ;
  uint32_t range = 0 ;
  int32_t  i, nbt ;
  uint32_t token, mask ;
  uint64_t bneeded ;

  for(i = 1 ; i < n ; i++) {   // find minimum and maximum values
    min32 = (src[i] < min32) ? src[i] : min32 ;
    max32 = (src[i] > max32) ? src[i] : max32 ;
  }     // find extrema
  range = max32 - min32 ;
  nbt = 0 ;
  while(range > 0) {           // compute number of bits necessary to encode the data range
    range >>= 1 ;
    nbt ++ ;
  }
  range = max32 - min32 ;
  mask = 0 ; mask = ~mask ; mask = mask >> (32 - nbt) ;  // nbt bits set to 1 on the right

  bneeded = nbt ;
  bneeded = bneeded * n ;               // bits needed = nbt * number of points
  bneeded += (16 + 8 + 8 + 32 + 64) ;   // add header size
  bneeded = (bneeded + 31) / 32 ;       // in uint32_t units
  if(bneeded > bavail) return -1 ;      // ERROR: not enough space for bit stream

  pstream_init(ps, dest) ;
  pstream_put_32(ps, 0xABBA, 16) ;      // signature
  pstream_put_32(ps, nbt, 8) ;          // number of bits per token
  pstream_put_32(ps, 32, 8) ;           // number of bits for minimum value is 32
  pstream_put_32(ps, n, 32) ;           // number of elements
  pstream_put_32(ps, min32, 32) ;       // minimum value

  for(i = 0 ; i < n ; i++){
    token = src[i] - min32 ; token &= mask ;
    pstream_put_32(ps, token, nbt) ;
  }

  pstream_flush(ps) ;
  return (ps->p - ps->s) ;
}

#if defined(SELF_TEST)
#define NPTS 17
int main(){
  uint64_t pak64[NPTS] ;
  int64_t pak64m[NPTS] ;
  uint64_t pak64a[NPTS] ;
  int64_t pak64b[NPTS] ;
  uint32_t pak32[NPTS] ;
  int32_t pak32m[NPTS] ;
  uint32_t pak32a[NPTS] ;
  int32_t pak32b[NPTS] ;
//   uint32_t packed[NPTS * 4] ;
  uint32_t packed1[NPTS * 4] ;
  uint32_t packed2[NPTS * 4] ;
  uint32_t packed3[NPTS * 4] ;
  uint32_t packed4[NPTS * 4] ;
  int i, n2 ;
  pstream ps ;

  for(i=0 ; i< NPTS ; i++) pak64[i] = i + 0xFFFF00000000L ;
  for(i=0 ; i< NPTS ; i++) pak64m[i] = i - 0x0FFF00000000L ;
  pak64[NPTS-1] += 0x80 ;
  pak64m[NPTS-1] += 0x80 ;

  for(i=0 ; i< NPTS ; i++) pak32[i] = i + 0xFFFF0000L ;
  for(i=0 ; i< NPTS ; i++) pak32m[i] = i - 0x0FFF0000L ;
  pak32[NPTS-1] += 0x80 ;
  pak32m[NPTS-1] += 0x80 ;
//   n2 = stream_pack_u64(pak64, packed, NPTS) ;
//   printf("n = %d\n",n2) ;
//   for(i=0 ; i<n2 ; i++) printf("%8.8x",packed[i]) ;
//   printf("\n") ;
  n2 = pstream_pack_u64(&ps, pak64, packed1, NPTS, NPTS * 4) ;
  printf("n = %d\n",n2) ;
  for(i=0 ; i<n2 ; i++) printf("%8.8x",packed1[i]) ;
  printf("\n") ;

  n2 = pstream_pack_i64(&ps, pak64m, packed2, NPTS, NPTS * 4) ;
  printf("n = %d\n",n2) ;
  for(i=0 ; i<n2 ; i++) printf("%8.8x",packed2[i]) ;
  printf("\n") ;

  n2 = pstream_unpack_u64(&ps, packed1, pak64a, NPTS) ;  // extract 64 bit stream into 64 bit unsigned tokens
  n2 = pstream_unpack_i64(&ps, packed2, pak64b, NPTS) ;  // extract 64 bit stream into 32 bit signed tokens
  for(i=0 ; i<NPTS ; i++){
    printf("%16.16lx %16.16lx %16.16lx %16.16lx \n",pak64[i],pak64a[i],pak64m[i],pak64b[i]);
  }

  n2 = pstream_pack_u32(&ps, pak32, packed3, NPTS, NPTS * 4) ;
  printf("n = %d\n",n2) ;
  for(i=0 ; i<n2 ; i++) printf("%8.8x",packed3[i]) ;
  printf("\n") ;

  n2 = pstream_pack_i32(&ps, pak32m, packed4, NPTS, NPTS * 4) ;
  printf("n = %d\n",n2) ;
  for(i=0 ; i<n2 ; i++) printf("%8.8x",packed4[i]) ;
  printf("\n") ;

//   pstream_rewind(&ps) ;
  n2 = pstream_unpack_u32(&ps, packed3, pak32a, NPTS) ;  // extract from 32 bit stream into 32 bit unsigned tokens
  n2 = pstream_unpack_u64(&ps, packed3, pak64a, NPTS) ;  // extract from 32 bit stream into 64 bit unsigned tokens
  n2 = pstream_unpack_i32(&ps, packed4, pak32b, NPTS) ;  // extract from 32 bit stream into 32 bit signed tokens
  n2 = pstream_unpack_i64(&ps, packed4, pak64b, NPTS) ;  // extract from 32 bit stream into 64 bit signed tokens
  for(i=0 ; i<NPTS ; i++){
    printf("%8.8x %8.8x %16.16lx %8.8x %8.8x %16.16lx\n",pak32[i],pak32a[i],pak64a[i],pak32m[i],pak32b[i],pak64b[i]);
  }
  n2 = pstream_unpack_u32(&ps, packed1, pak32a, NPTS) ;  // extract from 64 bit stream into 32 bit unsigned tokens
  if(n2 <= 0) printf("expected failure from 64 bit stream into 32 bit unsigned tokens\n") ;
  else printf("unexpectred n2 = %d from 64 bit stream into 32 bit unsigned tokens\n",n2);
  n2 = pstream_unpack_i32(&ps, packed2, pak32b, NPTS) ;  // extract from 64 bit stream into 32 bit signed tokens
  if(n2 <= 0) printf("expected failure from 64 bit stream into 32 bit signed tokens\n") ;
  else printf("unexpectred n2 = %d from 64 bit stream into 32 bit signed tokens\n",n2);

  return 0 ;
}
#endif
