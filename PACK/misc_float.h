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
//     M. Valin,   Recherche en Prevision Numerique, august 2022
//

void fp32_bf16(void *f32, uint32_t *u16, int32_t n) ;
void bf16_fp32(void *f32, uint32_t *u16, int32_t n) ;

void fp32_nf16(void *f32, uint32_t *u16, int32_t n, int32_t nexp) ;
void nf16_fp32(void *f32, uint32_t *u16, int32_t n, int nexp) ;

void fp32_bf24(void *f32, uint32_t *u24, int32_t n) ;
void bf24_fp32(void *f32, uint32_t *u24, int32_t n) ;

void u32_u24(uint32_t *u32, uint32_t *u24, int32_t n) ;
void u24u32(uint32_t *u32, uint32_t *u24, int32_t n) ;

void i32_u24(int32_t *i32, uint32_t *u24, int32_t n) ;
void u24i32(int32_t *i32, uint32_t *u24, int32_t n) ;
