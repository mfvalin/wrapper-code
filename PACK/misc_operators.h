// Hopefully useful code for C
// Copyright (C) 2022  Recherche en Prevision Numerique
//
// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//

#if ! defined(MISC_OPERATORS)
#define MISC_OPERATORS

#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define ABS(val)  ((val) < 0) ? (-(val)) : (val)
#define MINMAX(min,max,src,n)  { int i=1 ; min = max = src[0] ; while(i++ < n) { min = MIN(min,src[i]) ; max = MAX(max,src[i]); } }

#endif
