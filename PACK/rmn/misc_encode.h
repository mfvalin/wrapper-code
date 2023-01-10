/* Hopefully useful routines for C and FORTRAN
 * Copyright (C) 2020-2023  Recherche en Prevision Numerique
 *
 * This is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * Author : M. Valin (RPN-SI)
 */

#include <rmn/bi_endian_pack.h>

uint32_t stream_get_block_8x8(uint32_t * restrict src, int lni, uint32_t * restrict block, uint32_t * restrict po, uint32_t * restrict gain );
void stream_encode_init(bitstream *bstream, void *buffer, size_t bufsize);
uint32_t stream_encode_ublock(uint32_t *src, int nx, int ny, int nbits, bitstream *bstream);
uint32_t stream_decode_ublock(uint32_t *dst, bitstream *bstream);
