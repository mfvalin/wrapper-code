#include <stdio.h>

void fortran_constructor() ;

void __attribute__ ((constructor)) FortranConstructor(void) {
   printf("FortranConstructor calling fortran_constructor\n");
   fortran_constructor() ;  // call user supplied Fortran constructor for library
}
