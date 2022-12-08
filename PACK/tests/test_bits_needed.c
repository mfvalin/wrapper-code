#include <stdint.h>
#include <stdio.h>
#include <immintrin.h>
#include <misc_operators.h>

int main(int argc, char **argv){
  int i, n=0 ;
  for(i=1 ; i<33 ; i++){
    printf(" %11d %8.8x (%2d) %11d %8.8x (%2d) %11d %8.8x (%2d) %11d %8.8x (%2d) %11d %8.8x (%2d) %11d %8.8x (%2d)\n", 
           n   , n   , BitsNeeded_32(n),
           n+1 , n+1 , BitsNeeded_32(n+1),
           -n+1, -n+1, BitsNeeded_32(-n+1),
           -n  , -n  , BitsNeeded_32(-n),
           -n-1, -n-1, BitsNeeded_32(-n-1),
           -n-2, -n-2, BitsNeeded_32(-n-2)
          );
    n = (n<<1) | 1 ;
  }
#if 0
  n = 0 ;
  for(i=1 ; i<33 ; i++){
    printf(" %11d %8.8x (%2d) %11d %8.8x (%2d) %11d %8.8x (%2d) %11d %8.8x (%2d) %11d %8.8x (%2d) %11d %8.8x (%2d)\n", 
           n   , n   , BitsNeeded32(n),
           n+1 , n+1 , BitsNeeded32(n+1),
           -n+1, -n+1, BitsNeeded32(-n+1),
           -n  , -n  , BitsNeeded32(-n),
           -n-1, -n-1, BitsNeeded32(-n-1),
           -n-2, -n-2, BitsNeeded32(-n-2)
          );
    n = (n<<1) | 1 ;
  }
#endif
}
