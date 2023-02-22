/*
Hopefully useful code for C and Fortran
Copyright (C) 2022  Recherche en Prevision Numerique

This code is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation,
version 2.1 of the License.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.
*/

#if ! defined(MISC_PROPERTIES)
#define MISC_PROPERTIES

#include <rmn/misc_types.h>

ieee_prop_64 get_ieee_properties(void *restrict blk, int ni, int nj);

#endif
