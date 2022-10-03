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
//     M. Valin,   Recherche en Prevision Numerique, 2022/09
//
#include <stdint.h>

#define STATIC extern
#include <bi_endian_pack.h>

void  LeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw){
  int i = 0, n = (nw < 0) ? -nw : nw ;
  uint64_t  accum = p->accum ;
  uint32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;
  uint32_t mask = RMask(nbits) ;

  if(nbits <= 16) {
    uint32_t t, nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      // little endian => upper part [i+1] | lower part [i]
      t  = (w32[i  ] & mask) | ((w32[i+1] & mask) << nbits) ;
      LE64_PUT_NBITS(accum, insert, t, nb) ;
    }
  }
  for(    ; i<n ; i++){
    LE64_PUT_NBITS(accum, insert, w32[i], nbits) ;
  }
  if(nw <= 0) LE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
}

void  BeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw){
  int i = 0, n = (nw < 0) ? -nw : nw ;
  uint64_t  accum = p->accum ;
  uint32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;

  if(nbits <= 16) {
    uint32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      // big endian => upper part [i] | lower part [i+1]
      t  = (w32[i+1] & mask) | ((w32[i  ] & mask) << nbits) ;
      BE64_PUT_NBITS(accum, insert, t, nb) ;
    }
  }
  for(    ; i<n ; i++){
    BE64_PUT_NBITS(accum, insert, w32[i], nbits) ;
  }
  if(nw <= 0) BE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
}

void  LeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n){
  int i = 0 ;
  uint64_t  accum = p->accum ;
  uint32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(nbits <= 16) {
    uint32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      LE64_GET_NBITS(accum, xtract, t, nb) ;
      w32[i  ] = t & mask ;    // little endian means lower part first
      w32[i+1] = t >> nbits ;  // then upper part
    }
  }
  for(    ; i<n ; i++){
    LE64_GET_NBITS(accum, xtract, w32[i], nbits) ;
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
}

void  BeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n){
  int i = 0 ;
  uint64_t  accum = p->accum ;
  uint32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(nbits <= 16) {
    uint32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      BE64_GET_NBITS(accum, xtract, t, nb) ;
      w32[i  ] = t >> nbits ;    // big endian means upper part first
      w32[i+1] = t & mask ;      // then lower part
    }
  }
  for(    ; i<n ; i++){
    BE64_GET_NBITS(accum, xtract, w32[i], nbits) ;
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
}

void  BeStreamXtractSigned(bitstream *p, int32_t *w32, int nbits, int n){
  int i = 0 ;
  int64_t  accum = (int64_t)p->accum ;
  uint32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;

  if(nbits <= 16) {
    int32_t t, mask = RMask(nbits), nb = nbits + nbits ;
    for(    ; i<n-1 ; i+=2){
      BE64_GET_NBITS(accum, xtract, t, nb) ;
      w32[i  ] = t >> nbits ;    // big endian means upper part first
//       w32[i+1] = t & mask ;      // then lower part
      w32[i+1] = (t << (32-nbits)) >> (32-nbits) ;      // then lower part
    }
  }
  for(    ; i<n ; i++){
    BE64_GET_NBITS(accum, xtract, w32[i], nbits) ;
  }
  p->accum = (uint64_t)accum ;
  p->xtract = xtract ;
  p->stream = stream ;
}

void  LeStreamInsertM(bitstream *p, uint32_t *w32, int *nbits, int *n){
  int i ;
  uint64_t  accum = p->accum ;
  uint32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;
  while(n[0] > 0 ){
    for(i=0 ; i<n[0] ; i++){
      LE64_PUT_NBITS(accum, insert, w32[i], nbits[0]) ;
    }
    nbits++ ;
    n++ ;
  }
  if(n[0] == -1) LE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
}

void  BeStreamInsertM(bitstream *p, uint32_t *w32, int *nbits, int *n){
  int i ;
  uint64_t  accum = p->accum ;
  uint32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;
  while(n[0] > 0 ){
    for(i=0 ; i<n[0] ; i++){
      BE64_PUT_NBITS(accum, insert, w32[i], nbits[0]) ;
    }
    nbits++ ;
    n++ ;
  }
  if(n[0] == -1) BE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
}

void  LeStreamXtractM(bitstream *p, uint32_t *w32, int *nbits, int *n){
  int i ;
  uint64_t  accum = p->accum ;
  uint32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;
  while(n[0] > 0 ){
    for(i=0 ; i<n[0] ; i++){
      LE64_GET_NBITS(accum, xtract, w32[i], nbits[0]) ;
    }
    nbits++ ;
    n++ ;
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
}

void  BeStreamXtractM(bitstream *p, uint32_t *w32, int *nbits, int *n){
  int i ;
  uint64_t  accum = p->accum ;
  uint32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;
  while(n[0] > 0 ){
    for(i=0 ; i<n[0] ; i++){
      BE64_GET_NBITS(accum, xtract, w32[i], nbits[0]) ;
    }
    nbits++ ;
    n++ ;
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
}

