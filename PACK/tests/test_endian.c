#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#if ! defined(USE_INLINE)
#define PROTOTYPES_ONLY
#endif

#include <fast_endian.h>
#include <fast_endian.h>

#include <misc_timers.h>
#define NTIMES 100000

// make sure that NPTS is (multiple of 32) - 1 to exercize all code branches
#define NPTS (127)
#define NPTS16 (NPTS/2)
#define NPTS32 (NPTS/4)
#define NPTS64 (NPTS/8)

// xor indexing using ^1 makes 0 1 2 3 4 5 6 7 become 1 0 3 2 5 4 7 6
// xor indexing using ^3 makes 0 1 2 3 4 5 6 7 become 3 2 1 0 7 6 5 4
// xor indexing using ^7 makes 0 1 2 3 4 5 6 7 become 7 6 5 4 3 2 1 0
int main(int argc, char **argv){
  uint8_t c[NPTS+8];
  int64_t i;
  int errors, total_errors ;
  uint16_t *u16 = (uint16_t *) c;
  uint32_t *u32 = (uint32_t *) c;
  uint64_t *u64 = (uint64_t *) c;
  uint64_t t0, t1, t[NTIMES], freq, tmin, tmax ;
  double nano ;
  uint64_t in_cache[4096] ;
  uint64_t in_memory[4096*4096*2] ;
  int nelem ;

  freq = cycles_counter_freq() ;
  nano = 1000000000 ;
  nano /= freq ;
#if defined(__AVX2__) && defined(__x86_64__) && defined(WITH_SIMD)
  printf("using SIMD version\n");
#endif

  total_errors = 0 ;
  // 8 bits in 16 bits (2 tokens to swap)
  for(i=0;i<NPTS;i++) c[i] = i & 0xFF;  u16[NPTS16] = 0x1234;  // add verifyable guard at end of stream
  Swap_8_in_16(u16,u16,NPTS16);
  errors = 0;
  for(i=0;i<NPTS16*2;i++) { if( (i & 0xFF) != c[i^1] ) errors++; }  // xor indexing to check order 2 swap
  total_errors += errors ;
  printf("errors in Swap_8_in_16 fwd = %d, last = %4.4x\n",errors,u16[NPTS16]);
  Swap_8_in_16(u16,u16,NPTS16);
  errors = 0;
  for(i=0;i<NPTS16*2;i++) { if( (i & 0xFF) != c[i] ) errors++; }
  total_errors += errors ;
  printf("errors in Swap_8_in_16 bak = %d, last = %4.4x\n",errors,u16[NPTS16]);

  // 8 bits in 32 bits (4 tokens to swap)
  for(i=0;i<NPTS;i++) c[i] = i & 0xFF;  u32[NPTS32] = 0x12345678;  // add verifyable guard at end of stream
  Swap_8_in_32(u32,u32,NPTS32);
  errors = 0;
  for(i=0;i<NPTS32*4;i++) { if( (i & 0xFF) != c[i^3] ) errors++; }  // xor indexing to check order 4 swap
  total_errors += errors ;
  printf("errors in Swap_8_in_32 fwd = %d, last = %8.8x\n",errors,u32[NPTS32]);
  Swap_8_in_32(u32,u32,NPTS32);
  errors = 0;
  for(i=0;i<NPTS32*4;i++) { if( (i & 0xFF) != c[i] ) errors++; }
  total_errors += errors ;
  printf("errors in Swap_8_in_32 bak = %d, last = %8.8x\n",errors,u32[NPTS32]);

  // 16 bits in 32 bits (2 tokens to swap)
  for(i=0;i<NPTS16;i++) u16[i] = i & 0xFFFF;  u32[NPTS32] = 0x12345678;  // add verifyable guard at end of stream
  Swap_16_in_32(u32,u32,NPTS32);
  errors = 0;
  for(i=0;i<NPTS32*2;i++) { if( (i & 0xFFFF) != u16[i^1] ) errors++; }  // xor indexing to check order 2 swap
  total_errors += errors ;
  printf("errors in Swap_16_in_32 fwd = %d, last = %8.8x\n",errors,u32[NPTS32]);
  Swap_16_in_32(u32,u32,NPTS32);
  errors = 0;
  for(i=0;i<NPTS32*2;i++) { if( (i & 0xFF) != u16[i] ) errors++; }
  total_errors += errors ;
  printf("errors in Swap_16_in_32 bak = %d, last = %8.8x\n",errors,u32[NPTS32]);

  // 8 bits in 64 bits (8 tokens to swap)
  for(i=0;i<NPTS;i++) c[i] = i & 0xFF;  u64[NPTS64] = 0x12345678ABCDEF01;  // add verifyable guard at end of stream
  Swap_8_in_64(u64,u64,NPTS64);
  errors = 0;
  for(i=0;i<NPTS64*8;i++) { if( (i & 0xFF) != c[i^7] ) errors++; }  // xor indexing to check order 8 swap
  total_errors += errors ;
  printf("errors in Swap_8_in_64 fwd = %d, last = %16.16lx\n",errors,u64[NPTS64]);
  Swap_8_in_64(u64,u64,NPTS64);
  errors = 0;
  for(i=0;i<NPTS64*8;i++) { if( (i & 0xFF) != c[i] ) errors++; }
  total_errors += errors ;
  printf("errors in Swap_8_in_64 bak = %d, last = %16.16lx\n",errors,u64[NPTS64]);

  // 16 bits in 64 bits (4 tokens to swap)
  for(i=0;i<NPTS16;i++) u16[i] = i & 0xFFFF;  u64[NPTS64] = 0x12345678ABCDEF01;  // add verifyable guard at end of stream
  Swap_16_in_64(u64,u64,NPTS64);
  errors = 0;
  for(i=0;i<NPTS64*4;i++) { if( (i & 0xFFFF) != u16[i^3] ) errors++; }  // xor indexing to check order 4 swap
  total_errors += errors ;
  printf("errors in Swap_16_in_64 fwd = %d, last = %16.16lx\n",errors,u64[NPTS64]);
  Swap_16_in_64(u64,u64,NPTS64);
  errors = 0;
  for(i=0;i<NPTS64*4;i++) { if( (i & 0xFFFF) != u16[i] ) errors++; }
  total_errors += errors ;
  printf("errors in Swap_16_in_64 bak = %d, last = %16.16lx\n",errors,u64[NPTS64]);

  // 32 bits in 64 bits (2 tokens to swap)
  for(i=0;i<NPTS64;i++) u64[i] = (i << 32) + (i+1);  u64[NPTS64] = 0x12345678ABCDEF01;  // add verifyable guard at end of stream
  Swap_32_in_64(u64,u64,NPTS64);
  errors = 0;
  for(i=0;i<NPTS64;i++) { if( u64[i] != ( ((i+1)<<32) + i) )  errors++; }
  total_errors += errors ;
  printf("errors in Swap_32_in_64 fwd = %d, last = %16.16lx\n",errors,u64[NPTS64]);
  Swap_32_in_64(u64,u64,NPTS64);
  errors = 0;
  for(i=0;i<NPTS64;i++) { if( u64[i] != ( (i << 32) + (i+1) ) ) errors++; }
  total_errors += errors ;
  printf("errors in Swap_32_in_64 bak = %d, last = %16.16lx\n",errors,u64[NPTS64]);

  if(total_errors > 0) return 1 ;

  tmin = 1000000000;
  tmax = 0 ;
  nelem = sizeof(in_cache)/(sizeof(int)) ;
  Swap_8_in_32(in_cache, in_cache, nelem);
  Swap_8_in_32(in_cache, in_cache, nelem);
  Swap_8_in_32(in_cache, in_cache, nelem);
  for(i=0 ; i<NTIMES ; i++){
    t0 = elapsed_cycles() ;
    Swap_8_in_32(in_cache, in_cache, nelem);
    t1 = elapsed_cycles() ;
    t[i] = t1 - t0 ;
    tmin = (t1 - t0 < tmin) ? t1 - t0 : tmin ;
    tmax = (t1 - t0 > tmax) ? t1 - t0 : tmax ;
  }
  printf("cache(min) Swap_8_in_32 = %ld cycles, %8.1f ns, %8.2f ns/word\n", tmin, tmin*nano, tmin*nano/nelem) ;
  printf("cache(max) Swap_8_in_32 = %ld cycles, %8.1f ns, %8.2f ns/word\n", tmax, tmax*nano, tmax*nano/nelem) ;

  nelem = sizeof(in_memory)/(sizeof(int)) ;
  t0 = elapsed_cycles() ;
  Swap_8_in_32(in_memory, in_memory, nelem);
  t1 = elapsed_cycles() ;
  tmin = t1 - t0 ;
  printf("memory(first) Swap_8_in_32 = %ld cycles, %8.1f ns, %8.2f ns/word\n", tmin, tmin*nano, tmin*nano/nelem) ;
  tmin = 1000000000;
  tmax = 0 ;
  for(i=0 ; i<NTIMES/1000 ; i++){
    t0 = elapsed_cycles() ;
    Swap_8_in_32(in_memory, in_memory, nelem);
    t1 = elapsed_cycles() ;
    t[i] = t1 - t0 ;
    tmin = (t1 - t0 < tmin) ? t1 - t0 : tmin ;
    tmax = (t1 - t0 > tmax) ? t1 - t0 : tmax ;
  }
  printf("memory(min) Swap_8_in_32 = %ld cycles, %8.1f ns, %8.2f ns/word\n", tmin, tmin*nano, tmin*nano/nelem) ;
  printf("memory(max) Swap_8_in_32 = %ld cycles, %8.1f ns, %8.2f ns/word\n", tmax, tmax*nano, tmax*nano/nelem) ;

  return 0;
}
