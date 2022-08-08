#include <stdio.h>
#include <misc_helpers.h>

int main(int argc, char **argv){
  uint64_t t1, t2, tc1, tc2, freq ;
  double cpi ;
  double cycle, ns ;
  int i=0;

// wall clock timer test
  t2 = t1 = elapsed_us() ;
  while(t2-t1 < 1000){
    t2 = elapsed_us() ;
    i++ ;
  }
  printf("timeofday iter = %d, %ld microseconds, ns/call = %ld\n", i, t2-t1, 1000 * (t2-t1) / i) ;

// try to determine the cycles counter frequency
  freq = cycles_counter_freq() ;
  cycle = 1.0 / freq ;
  ns = cycle * 1000000000 ;
  printf("\nestimated cycles counter frequency = %4.0f MHz\n\n", freq/1000000.0) ;

// cycle counters tests
  i = 0 ;
  t2 = t1 = elapsed_cycles_nofence() ;
  while(t2-t1 < 1000000){  // prime the pump before tests
    t2 = elapsed_cycles() ;
    i++ ;
  }
  cpi = (t2-t1) ; cpi /= i ;
//   printf("rdtscp        elapsed cycles iter = %6d, %9ld cycles, cycles/iter = %8.2f, ns/iter = %8.2f\n", i, t2-t1, cpi, ns*cpi) ;

  i = 0 ;
  t2 = t1 = elapsed_cycles_nofence() ;
  while(t2-t1 < 1000000){
    t2 = elapsed_cycles() ;
    i++ ;
  }
  cpi = (t2-t1) ; cpi /= i ;
  printf("rdtscp        elapsed cycles iter = %6d, %9ld cycles, cycles/iter = %8.2f, ns/iter = %8.2f\n", i, t2-t1, cpi, ns*cpi) ;

  i = 0 ;
  t2 = t1 = elapsed_cycles_fenced() ;
  while(t2-t1 < 1000000){
    t2 = elapsed_cycles_fenced() ;
    i++ ;
  }
  cpi = (t2-t1) ; cpi /= i ;
  printf("rdtscp+mfence elapsed cycles iter = %6d, %9ld cycles, cycles/iter = %8.2f, ns/iter = %8.2f\n", i, t2-t1, cpi, ns*cpi) ;

  i = 0 ;
  t2 = t1 = elapsed_cycles() ;
  while(t2-t1 < 1000000){
    t2 = elapsed_cycles() ;
    i++ ;
  }
  cpi = (t2-t1) ; cpi /= i ;
  printf("lfence+rtdtsc elapsed cycles iter = %6d, %9ld cycles, cycles/iter = %8.2f, ns/iter = %8.2f\n", i, t2-t1, cpi, ns*cpi) ;

  i = 0 ;
  t2 = t1 = elapsed_cycles_fast() ;
  while(t2-t1 < 1000000){
    t2 = elapsed_cycles_fast() ;
    i++ ;
  }
  cpi = (t2-t1) ; cpi /= i ;
  printf("rdtsc         elapsed cycles iter = %6d, %9ld cycles, cycles/iter = %8.2f, ns/iter = %8.2f\n", i, t2-t1, cpi, ns*cpi) ;

  return 0 ;
}
