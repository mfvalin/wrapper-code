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
  int i, n = (nw < 0) ? -nw : nw ;
  uint64_t  accum = p->accum ;
  uint32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;
  for(i=0 ; i<n ; i++){
    LE64_PUT_NBITS(accum, insert, w32[i], nbits) ;
  }
  if(nw <= 0) LE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
}

void  BeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nw){
  int i, n = (nw < 0) ? -nw : nw ;
  uint64_t  accum = p->accum ;
  uint32_t  insert = p->insert ;
  uint32_t *stream = p->stream ;
  for(i=0 ; i<n ; i++){
//     printf("accum = %16.16lx, token = %8.8x", accum, w32[i]) ;
    BE64_PUT_NBITS(accum, insert, w32[i], nbits) ;
//     printf(", insert = %3d, accum = %16.16lx, stream - start = %ld \n", insert, accum, stream - p->start) ;
  }
  if(nw <= 0) BE64_INSERT_FINAL(accum, insert, stream) ;
  p->accum = accum ;
  p->insert = insert ;
  p->stream = stream ;
}

void  LeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n){
  int i ;
  uint64_t  accum = p->accum ;
  uint32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;
  for(i=0 ; i<n ; i++){
    LE64_GET_NBITS(accum, xtract, w32[i], nbits) ;
  }
  p->accum = accum ;
  p->xtract = xtract ;
  p->stream = stream ;
}

void  BeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n){
  int i ;
  uint64_t  accum = p->accum ;
  uint32_t  xtract = p->xtract ;
  uint32_t *stream = p->stream ;
  for(i=0 ; i<n ; i++){
//     printf("accum = %16.16lx, xtract = %2d, stream = %8.8x", accum, xtract, *stream) ;
    BE64_GET_NBITS(accum, xtract, w32[i], nbits) ;
//     printf(", accum = %16.16lx, w32[i] = %8.8x, xtract = %d, stream = %8.8x\n", accum, w32[i], xtract, *stream) ;
  }
  p->accum = accum ;
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

