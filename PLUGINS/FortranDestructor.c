#include <stdio.h>

void fortran_destructor() ;

void __attribute__ ((destructor)) PluginDestructor(void) {
   printf("plugin destructor calling fortran_destructor\n");
   fortran_destructor() ;  // call user supplied Fortran constructor for library
}
