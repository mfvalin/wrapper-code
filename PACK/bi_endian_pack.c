//
// Copyright (C) 2022  Environnement Canada
//
// This is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// Author:
//     M. Valin,   Recherche en Prevision Numerique, 2022
//
// set of functions to help manage insertion/extraction of integers into/from a bit stream
// the bit stream is a sequence of unsigned 32 bit integers
// both Big and Little endian style insertion/extraction support is provided
//
#include <stdint.h>

#define STATIC extern
#include <rmn/bi_endian_pack.h>

// little endian style insertion of values into a bit stream
// p     : stream                                [INOUT]
// w32   : array of values to insert             [IN]
// nbits : number of bits kept for each value    [IN]
// nw    : number of values from w32             [IN}
int  LeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw){
  int i = 0, n = (nw < 0) ? -nw : nw ;
  uint64_t  accum = p->accum ;
  int32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;
  uint32_t mask = RMask(nbits) ;

  if(insert < 0) return 0;      // ERROR: not in insert mode

  if(nbits <= 16) {       // process values two at a time
    uint32_t t, nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      // little endian => upper part [i+1] | lower part [i]
      t  = (w32[i  ] & mask) | ((w32[i+1] & mask) << nbits) ;
      LE64_PUT_NBITS(accum, insert, t, nb, stream) ;   // insert a pair of values
    }
  }
  for(    ; i<n ; i++){
    LE64_PUT_NBITS(accum, insert, w32[i], nbits, stream) ;
  }
  if(nw <= 0) LE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
  return n ;
}

// big endian style insertion of values into a bit stream
// p     : stream                                [INOUT]
// w32   : array of values to insert             [IN]
// nbits : number of bits kept for each value    [IN]
// nw    : number of values from w32             [IN}
int  BeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw){
  int i = 0, n = (nw < 0) ? -nw : nw ;
  uint64_t  accum = p->accum ;
  int32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;

  if(insert < 0) return 0;      // ERROR: not in insert mode

  if(nbits <= 16) {       // process values two at a time
    uint32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      // big endian => upper part [i] | lower part [i+1]
      t  = (w32[i+1] & mask) | ((w32[i  ] & mask) << nbits) ;
      BE64_PUT_NBITS(accum, insert, t, nb, stream) ;   // insert a pair of values
    }
  }
  for(    ; i<n ; i++){
    BE64_PUT_NBITS(accum, insert, w32[i], nbits, stream) ;
  }
  if(nw <= 0) BE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
  return n ;
}

// little endian style extraction of unsigned values from a bit stream
// p     : stream                                [INOUT]
// w32   : array of unsigned values extracted    [OUT]
// nbits : number of bits kept for each value    [IN]
// n     : number of values from w32             [IN}
int  LeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n){
  int i = 0 ;
  uint64_t  accum = p->accum ;
  int32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(xtract < 0) return 0;      // ERROR: not in extract mode

  if(nbits <= 16) {       // process values two at a time
    uint32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      LE64_GET_NBITS(accum, xtract, t, nb, stream) ;   // get a pair of values
      w32[i  ] = t & mask ;                    // little endian means lower part first
      w32[i+1] = t >> nbits ;                  // then upper part
    }
  }
  for(    ; i<n ; i++){
    LE64_GET_NBITS(accum, xtract, w32[i], nbits, stream) ;
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
  return n ;
}

// little endian style extraction of signed values from a bit stream
// p     : stream                                [INOUT]
// w32   : array of unsigned values extracted    [OUT]
// nbits : number of bits kept for each value    [IN]
// n     : number of values from w32             [IN}
int  LeStreamXtractSigned(bitstream *p, int32_t *w32, int nbits, int n){
  int i = 0 ;
  int64_t  accum = (int64_t)p->accum ;
  int32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(xtract < 0) return 0;      // ERROR: not in extract mode

  if(nbits <= 16) {       // process values two at a time
    int32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      LE64_GET_NBITS(accum, xtract, t, nb, stream) ;   // get a pair of values
      // use shift to propagate sign
      w32[i  ] = (t << (32-nbits)) >> (32-nbits) ;     // little endian means lower part first
      w32[i+1] = t >> nbits ;                          // then upper part
    }
  }
  for(    ; i<n ; i++){
    LE64_GET_NBITS(accum, xtract, w32[i], nbits, stream) ;
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
  return n ;
}

// big endian style extraction of unsigned values from a bit stream
// p     : stream                                [INOUT]
// w32   : array of unsigned values extracted    [OUT]
// nbits : number of bits kept for each value    [IN]
// n     : number of values from w32             [IN}
int  BeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n){
  int i = 0 ;
  uint64_t  accum = p->accum ;
  int32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(xtract < 0) return 0;      // ERROR: not in extract mode

  if(nbits <= 16) {       // process values two at a time
    uint32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      BE64_GET_NBITS(accum, xtract, t, nb, stream) ;      // get a pair of values
      w32[i  ] = t >> nbits ;                     // big endian means upper part first
      w32[i+1] = t & mask ;                       // then lower part
    }
  }
  for(    ; i<n ; i++){
    BE64_GET_NBITS(accum, xtract, w32[i], nbits, stream) ;
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
  return n ;
}

// big endian style extraction of signed values from a bit stream
// p     : stream                                [INOUT]
// w32   : array of signed values extracted      [OUT]
// nbits : number of bits kept for each value    [IN]
// n     : number of values from w32             [IN}
int  BeStreamXtractSigned(bitstream *p, int32_t *w32, int nbits, int n){
  int i = 0 ;
  int64_t  accum = (int64_t)p->accum ;
  int32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(xtract < 0) return 0;      // ERROR: not in extract mode

  if(nbits <= 16) {       // process values two at a time
    int32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      BE64_GET_NBITS(accum, xtract, t, nb, stream) ;         // get a pair of values
      // use shift to propagate sign
      w32[i  ] = t >> nbits ;                        // big endian means upper part first
      w32[i+1] = (t << (32-nbits)) >> (32-nbits) ;   // then lower part
    }
  }
  for(    ; i<n ; i++){
    BE64_GET_NBITS(accum, xtract, w32[i], nbits, stream) ;
  }
  p->accum = (uint64_t)accum ;
  p->xtract = xtract ;
  p->stream = stream ;
  return n ;
}

// little endian style insertion of values into a bit stream
// p     : stream                                [INOUT]
// w32   : array of values inserted              [IN]
// nbits : array of nuber of bits to keep        [IN]
// n     : number of values to insert            [IN]
// n[i], nbits[i] is an associated (nbits , n) pair
// n[i] <= 0 marks the end of the pair list
int  LeStreamInsertM(bitstream *p, uint32_t *w32, int *nbits, int *n){
  int i, nw = 0 ;
  uint64_t  accum = p->accum ;
  int32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;

  if(insert < 0) return 0;      // ERROR: not in insert mode

  while(n[0] > 0 ){     // loop until end of list (non positive value)
    nw += n[0] ;
    for(i=0 ; i<n[0] ; i++){
      LE64_PUT_NBITS(accum, insert, w32[i], nbits[0], stream) ;
    }
    nbits++ ;
    n++ ;       // next pair
  }
  if(n[0] == -1) LE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
  return nw ;
}

// big endian style insertion of values into a bit stream
// p     : stream                                [INOUT]
// w32   : array of values inserted              [IN]
// nbits : array of nuber of bits to keep        [IN]
// n     : number of values to insert            [IN]
// n[i], nbits[i] is an associated (nbits , n) pair
// n[i] <= 0 marks the end of the pair list
int  BeStreamInsertM(bitstream *p, uint32_t *w32, int *nbits, int *n){
  int i, nw = 0 ;
  uint64_t  accum = p->accum ;
  int32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;

  if(insert < 0) return 0;      // ERROR: not in insert mode

  while(n[0] > 0 ){     // loop until end of list (non positive value)
    nw += n[0] ;
    for(i=0 ; i<n[0] ; i++){
      BE64_PUT_NBITS(accum, insert, w32[i], nbits[0], stream) ;
    }
    nbits++ ;
    n++ ;       // next pair
  }
  if(n[0] == -1) BE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
  return nw ;
}

// little endian style extraction of unsigned values from a bit stream
// p     : stream                                [INOUT]
// w32   : array of unsigned values extracted    [OUT]
// nbits : array of nuber of bits to keep        [IN]
// n     : number of unsigned values to extract  [IN]
// n[i], nbits[i] is an associated (nbits , n) pair
// n[i] <= 0 marks the end of the pair list
int  LeStreamXtractM(bitstream *p, uint32_t *w32, int *nbits, int *n){
  int i, nw = 0 ;
  uint64_t  accum = p->accum ;
  int32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(xtract < 0) return 0;      // ERROR: not in extract mode

  while(n[0] > 0 ){     // loop until end of list (non positive value)
    nw += n[0] ;
    for(i=0 ; i<n[0] ; i++){
      LE64_GET_NBITS(accum, xtract, w32[i], nbits[0], stream) ;
    }
    nbits++ ;
    n++ ;       // next pair
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
  return nw ;
}

// big endian style extraction of unsigned values from a bit stream
// p     : stream                                [INOUT]
// w32   : array of unsigned values extracted    [OUT]
// nbits : array of nuber of bits to keep        [IN]
// n     : number of unsigned values to extract  [IN]
// n[i], nbits[i] is an associated (nbits , n) pair
// n[i] <= 0 marks the end of the pair list
int  BeStreamXtractM(bitstream *p, uint32_t *w32, int *nbits, int *n){
  int i, nw = 0 ;
  uint64_t  accum = p->accum ;
  int32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(xtract < 0) return 0;      // ERROR: not in extract mode

  while(n[0] > 0 ){     // loop until end of list (non positive value)
    nw += n[0] ;
    for(i=0 ; i<n[0] ; i++){
      BE64_GET_NBITS(accum, xtract, w32[i], nbits[0], stream) ;
    }
    nbits++ ;
    n++ ;       // next pair
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
  return nw ;
}

