/*
 * Hopefully useful code for C and Fortran
 * Copyright (C) 2022  Recherche en Prevision Numerique
 *
 * This code is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This code is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * in the ARM V8 case, all the cycle counters use the same instruction
 * in the X86_64 case, various cycle counters and fences are used
 *
 * useful references :
 * https://sites.utexas.edu/jdm4372/2018/07/   John McCalpin's blog
 * https://github.com/jdmccalpin/low-overhead-timers
 */
#if defined(IN_FORTRAN_CODE)
interface
  function elapsed_us() result(t) bind(C,name='ElapsedUs')
    import C_INT64_T
    implicit none
    integer(C_INT64_T) :: t
  end function elapsed_us
  function elapsed_cycles() result(t) bind(C,name='ElapsedCycles')
    import C_INT64_T
    implicit none
    integer(C_INT64_T) :: t
  end function elapsed_cycles
  function cycles_counter_freq() result(t) bind(C,name='CyclesCounterFreq')
    import C_INT64_T
    implicit none
    integer(C_INT64_T) :: t
  end function cycles_counter_freq
end interface
#else

// tmin, tmax, tavg can be floats or doubles, niter is expected to be an integer
//
// timing loop top part
// niter : number of iterations
#define TIME_LOOP_TOP(niter) \
{ uint64_t t, to , mint = ~0, maxt = 0, avgt = 0.0 ; int iter = niter ; \
  for(j=0 ; j < iter ; j++) { to = elapsed_cycles() ; t = elapsed_cycles() ; to = t - to ;
//
// code to be timed in a loop goes between TIME_LOOP_TOP and TIME_LOOP_BOT
//
// timing loop bottom part
// tmin  : best timing in cycles
// tmax  : worst timing in cycles
// tavg  : average timing in cycles
// npts  : number of points processed
// buf   : buffer to receive diagnostic text
// bufsiz: size of buf
#define TIME_LOOP_BOT(tmin, tmax, tavg, npts, buf, bufsiz) \
    t = elapsed_cycles() -t -to ; \
    avgt += t ; mint = (t < mint) ? t : mint ; maxt = (t > maxt) ? t : maxt ; \
  } \
  tmin = mint ; tmax = maxt ; tavg = avgt/iter ; \
  if(npts > 0 && buf != NULL) \
    snprintf(buf, (size_t)bufsiz, " npts = %d, niter = %d,  ns= %6.0f (%6.0f), %6.2f ns/pt\n", \
             npts, niter, tmin*nano, tavg*nano, tavg*nano/npts) ; \
}

// time a piece of code in a loop, get min, max, average time
// tmin  : best timing in cycles
// tmax  : worst timing in cycles
// tavg  : average timing in cycles
// niter : number of iterations
// npts  : number of points processed
// TimedCode : code to be timed n (>=1) statement(s)
// buf   : buffer to receive diagnostic text
// bufsiz: size of buf
#define TIME_LOOP(tmin, tmax, tavg, niter, npts, buf, bufsiz, TimedCode) \
  TIME_LOOP_TOP(niter) ; \
  TimedCode ; \
  TIME_LOOP_BOT(tmin, tmax, tavg, npts, buf, bufsiz) ;

#if ! defined(MISC_TIMERS)
#define MISC_TIMERS

#include <stdint.h>
#include <sys/time.h>
#include <stddef.h>

static inline uint64_t elapsed_us(void){
  struct timeval t ;
  uint64_t elapsed ;
  gettimeofday(&t, NULL) ;
  elapsed = t.tv_sec ;
  elapsed *= 1000000 ;
  elapsed += t.tv_usec ;
  return elapsed ;
}

// NO serializing, NO fencing
static inline uint64_t elapsed_cycles_fast(void) {
#if defined(__x86_64__)
  uint64_t lo, hi ;
  __asm__ volatile ("rdtsc" : /* outputs   */ "=a" (lo), "=d" (hi) );
  return lo | (hi << 32);
#elif defined(__aarch64__)
  uint64_t time0 ;
  asm volatile ("isb ; mrs %0, cntvct_el0" : "=r" (time0));
  return time0;
#else
  return elapsed_us() * 1000 ;  // nanoseconds
#endif
}

// WITH serializing, NO fencing
static inline uint64_t elapsed_cycles_nofence(void) {
#if defined(__x86_64__)
  uint64_t lo, hi, misc ;
  __asm__ volatile ("rdtscp": /* outputs   */ "=a" (lo), "=d" (hi), "=c" (misc) );
  return lo | (hi << 32);
#elif defined(__aarch64__)
  uint64_t time0 ;
  asm volatile ("isb ; mrs %0, cntvct_el0" : "=r" (time0));
  return time0;
#else
  return elapsed_us() * 1000 ;  // nanoseconds
#endif
}

// WITH serializing, WITH memory fencing after
static inline uint64_t elapsed_cycles_fenced(void) {
#if defined(__x86_64__)
  uint64_t lo, hi, misc ;
  __asm__ volatile ("rdtscp": /* outputs   */ "=a" (lo), "=d" (hi), "=c" (misc) );
  __asm__ volatile ("mfence");
  return lo | (hi << 32);
#elif defined(__aarch64__)
  uint64_t time0 ;
  asm volatile ("isb ; mrs %0, cntvct_el0" : "=r" (time0));
  return time0;
#else
  return elapsed_us() * 1000 ;  // nanoseconds
#endif
}

// NO serializing, LOCAL fencing before
// default form
static inline uint64_t elapsed_cycles(void) {
#if defined(__x86_64__)
  uint64_t lo, hi ;
  __asm__ volatile ("lfence");
  __asm__ volatile ("rdtsc" : /* outputs   */ "=a" (lo), "=d" (hi) );
  return lo | (hi << 32);
#elif defined(__aarch64__)
  uint64_t time0 ;
  asm volatile ("isb ; mrs %0, cntvct_el0" : "=r" (time0));
  return time0;
#else
  return elapsed_us() * 1000 ;  // nanoseconds
#endif
}

static inline uint64_t cycles_counter_freq(void){
  static uint64_t timerfreq = 0;
  uint64_t t1, t2, tc1, tc2 ;

  if(timerfreq != 0) return timerfreq ;

#if defined(__aarch64__)
  asm volatile("isb ; mrs %0, cntfrq_el0" : "=r"(timerfreq));
  return timerfreq ;
#endif

  t2 = t1 = elapsed_us() ;
  tc1 = elapsed_cycles() ;
  while(t2-t1 < 1000){    // 1 millisecond
    t2 = elapsed_us() ;
    tc2 = elapsed_cycles() ;
  }
  timerfreq = 1000000 * (tc2-tc1) / (t2-t1);
  return timerfreq ;
}

uint64_t ElapsedUs(void);
uint64_t ElapsedCycles(void);
uint64_t CyclesCounterFreq(void);

#endif

#endif
