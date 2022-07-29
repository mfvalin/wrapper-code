
#include <stdint.h>

// a[8] , b[8] -> c[4] , row length of a = 8
extern inline void average_rows_8x2_8(float *a, float *c){
  int i, j ;
  float t[8] ;
  for(i = 0 ; i < 8 ; i++){           // add rows
    t[i] = a[i] + a[i+8] ;
  }
  for(i =0, j=0 ; j<4 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[8] , b[8] -> c[4] , row length of a = lrow
extern inline void average_rows_8x2(float *a, float *c, int lrow){
  int i, j ;
  float t[8] ;
  for(i = 0 ; i < 8 ; i++){           // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<4 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[16,2] -> c[8] , row length of a = 16
extern inline void average_rows_16x2_16(float *a, float *c){
  int i, j ;
  float t[16] ;
  for(i = 0 ; i < 16 ; i++){           // add rows
    t[i] = a[i] + a[i+16] ;
  }
  for(i =0, j=0 ; j<8 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[16,2] -> c[8] , row length of a = lrow
extern inline void average_rows_16x2(float *a, float *c, int lrow){
  int i, j ;
  float t[16] ;
  for(i = 0 ; i < 16 ; i++){           // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<8 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[32,2] -> c[16] , row length of a = 32
extern inline void average_rows_32x2_32(float *a, float *c){
  int i, j ;
  float t[32] ;
  for(i = 0 ; i < 32 ; i++){           // add rows
    t[i] = a[i] + a[i+32] ;
  }
  for(i =0, j=0 ; j<16 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[32,2] -> c[16] , row length of a = lrowlength of a = lrow
extern inline void average_rows_32x2(float *a, float *c, int lrow){
  int i, j ;
  float t[32] ;
  for(i = 0 ; i < 32 ; i++){           // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<16 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[64,2] -> c[32] , row length of a = 64
extern inline void average_rows_64x2_64(float *a, float *c){
  int i, j ;
  float t[64] ;
  for(i = 0 ; i < 64 ; i++){           // add rows
    t[i] = a[i] + a[i+64] ;
  }
  for(i =0, j=0 ; j<32 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[64,2] -> c[32] , row length of a = lrow
extern inline void average_rows_64x2(float *a, float *c, int lrow){
  int i, j ;
  float t[64] ;
  for(i = 0 ; i < 64 ; i++){           // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<32 ; j++, i+=2){   // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[n,2] -> c[n] , row length of a = lrow
extern inline void average_rows_nx2(float *a, float *c, int n, int lrow){
  int i, j ;
  float t[8] ;
  while(n >= 64){                        // slices of [64,2] -> [32]
    average_rows_64x2(a, c, lrow) ;
    n -= 64 ; a+= 64 ; c += 32 ;         // next slice
  }
  if(n >= 32){                           // slices of [32,2] -> [16]
    average_rows_32x2(a, c, lrow) ;
    n -= 32 ; a+= 32 ; c += 16 ;         // next slice
  }
  if(n >= 16){                           // slices of [16,2] -> [8]
    average_rows_32x2(a, c, lrow) ;
    n -= 16 ; a+= 16 ; c += 8 ;          // next slice
  }
  if(n >= 8){                            // slices of [8,2] -> [4]
    average_rows_32x2(a, c, lrow) ;
    n -= 8 ; a+= 8 ; c += 4 ;            // next slice
  }
  // process last (<8) points
  for(i = 0 ; i < n ; i++){              // add rows
    t[i] = a[i] + a[i+lrow] ;
  }
  for(i =0, j=0 ; j<n/2 ; j++, i+=2){    // add even terms with odd terms in row
    c[j] = t[i] + t[i+1] ;
    c[j] *= .25f ;
  }
}

// a[64,8] -> a2[32,4] -> a4[16,2] -> a8[8] , row length of a = 64
void averages_64x8_64(float *a, float *a2, float *a4, float *a8){
  // 8 lines of a -> 4 lines of a2
  average_rows_64x2_64(a       , a2     ) ;  // [64,2] -> [32]
  average_rows_64x2_64(a  + 128, a2 + 32) ;  // [64,2] -> [32]
  average_rows_64x2_64(a  + 256, a2 + 64) ;  // [64,2] -> [32]
  average_rows_64x2_64(a  + 384, a2 + 96) ;  // [64,2] -> [32]
  // 4 lines of a2 -> 2 lines of a4
  average_rows_32x2_32(a2      , a4     ) ;  // [32,2] -> [16]
  average_rows_32x2_32(a2 +  64, a4 + 16) ;  // [32,2] -> [16]
  // 2 lines of a4 -> 1 line of a8
  average_rows_16x2_16(a4      , a8     ) ;  // [16,2] -> [8]
}

// a[64,8] -> a2[32,4] -> a4[16,2] -> a8[8] , row length of a = ltow
void averages_64x8(float *a, float *a2, float *a4, float *a8, int lrow){
  // 8 lines of a -> 4 lines of a2
  average_rows_64x2(a       , a2     , lrow) ;  // [64,2] -> [32]
  average_rows_64x2(a  + 128, a2 + 32, lrow) ;  // [64,2] -> [32]
  average_rows_64x2(a  + 256, a2 + 64, lrow) ;  // [64,2] -> [32]
  average_rows_64x2(a  + 384, a2 + 96, lrow) ;  // [64,2] -> [32]
  // 4 lines of a2 -> 2 lines of a4
  average_rows_32x2(a2      , a4     , lrow) ;  // [32,2] -> [16]
  average_rows_32x2(a2 +  64, a4 + 16, lrow) ;  // [32,2] -> [16]
  // 2 lines of a4 -> 1 line of a8
  average_rows_16x2(a4      , a8     , lrow) ;  // [16,2] -> [8]
}

// a[64,64] a2[32,32]  a4[16,16]  a8[8,8] , row length of a = 64
void averages_64x64_64(float *a, float *a2, float *a4, float *a8){
  int i ;
  for(i = 0 ; i < 7 ; i++){   // 8 blocks of a[64,8]  a2[32,4]  a4[16,2]  a8[8]
    averages_64x8_64(a + i*64*8, a2 + i*32*4, a4 + i*16*2, a8 + i*8) ;
  }
}

// a[64,64] a2[32,32]  a4[16,16]  a8[8,8] , row length of a = lrow
void averages_64x64(float *a, float *a2, float *a4, float *a8, int lrow){
  int i ;
  for(i = 0 ; i < 7 ; i++){   // 8 blocks of a[64,8]  a2[32,4]  a4[16,2]  a8[8]
    averages_64x8(a + i*64*8, a2 + i*32*4, a4 + i*16*2, a8 + i*8, lrow) ;
  }
}

