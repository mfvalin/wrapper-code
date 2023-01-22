#include <stdint.h>
#include <stdio.h>

#include <misc_operators.h>

int main(int argc, char **argv){
  int32_t x ;
  printf(" /2 /4 /8       truncate             round             divide\n") ;
  for (x = -17 ; x < 18 ; x++){
    printf("x = %4d | %4d, %4d, %4d | %4d, %4d, %4d | %4d, %4d, %4d\n", x, IDIV2T(x), IDIV4T(x), IDIV8T(x), IDIV2R(x), IDIV4R(x), IDIV8R(x), x/2, x/4, x/8) ;
  }
}
