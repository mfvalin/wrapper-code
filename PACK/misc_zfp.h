//  Copyright (C) 2022  Environnement Canada
//
//  This is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation,
//  version 2.1 of the License.
//
//  This software is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//

#include <zfp.h>

int ZfpCompress_set_debug(int flag);
int ZfpCompress_set_diag(int flag);
void *ZfpCompress(void* array, int nx, int ny, int nz, float maxerror, int precision , int *Ssize);
int32_t ZfpArrayDims(void *Stream, int *d, int *s, int Ssize);
int32_t ZfpExpand(void* array, int nx, int ny, int nz, void *Stream, int Ssize);
