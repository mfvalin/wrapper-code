#include <stdio.h>
#include <stdlib.h>

#include <misc_pack.h>
#include <misc_operators.h>
#include <misc_timers.h>
#include <misc_timers.h>

int main(int argc, char **argv){
  float *zi ;
  float *zo ;
  int *iw ;
  QuantizeHeader p;
  int i, j;
  double toler;
  double avg, avgi, avgo ;
  int error=0;
  float minval, maxval ;
  double delta, maxdelta ;
  double avgdelta = 0.0 ;
  int NBITS = 22 ;
  float quantum ;
  float epsilon = .00000006f;
  float eps_plus = 1.0f + epsilon ;
  float eps_minus = 1.0f - epsilon ;
  int NPTS = 32800 ;
  int nloops = 1000 ;
  uint64_t freq ;
  TIME_LOOP_DATA ;

  freq = cycles_counter_freq() ;
  printf("time base : %6.2G ns / tick, %6.2G GHz\n", 1.0E9/freq, freq/1.0E9) ;

  if(argc > 1){
    NPTS= atoi(argv[1]) * 64 ;
    if(NPTS<19) NPTS = 19 ;
  }
  zi = (float *) malloc(NPTS*sizeof(float)) ;
  zo = (float *) malloc(NPTS*sizeof(float)) ;
  iw = (int *) malloc(NPTS*sizeof(float)) ;

  printf("=============== rounding test with 19 points =================\n");
  quantum = .25f ;
  printf("eps = %15.9f, eps_plus = %15.9f, eps_minus = %15.9f, quantum =%15.9f \n",
         epsilon, eps_plus, eps_minus, quantum) ;
  printf("(i - 9) * .125f\n");
  for(i=0 ; i<19 ; i++){ zi[i] = (i - 9) * .125f ; /* zi[i] *= 8.0f */ ; printf("%8.3f", zi[i]) ; } printf("\n") ;
  printf("quantizing/restoring (i - 9) * .125f\n");
  float_quantize_prep(NBITS, &p, zi[18], zi[0], quantum);
  float_quantize(iw, zi, 19, 19, 19, 1, &p ) ;
  float_unquantize(iw, zo, 19, 19, 19, 1, &p ) ;
  for(i=0 ; i<19 ; i++) printf("%8.3f", zo[i]) ; printf("\n") ;
  for(i=0 ; i<19 ; i++) printf("%8d", iw[i]) ; printf("\n") ;

  printf("quantizing/restoring (i - 9) * .125f * (1 - epsilon)\n");
  for(i=0 ; i<19 ; i++){ zi[i] = (i - 9) * .125f ; zi[i] *= eps_minus ; }
  float_quantize_prep(NBITS, &p, 1.125f, -1.125f, quantum);
  float_quantize(iw, zi, 19, 19, 19, 1, &p ) ;
  float_unquantize(iw, zo, 19, 19, 19, 1, &p ) ;
  for(i=0 ; i<19 ; i++) printf("%8.3f", zo[i]) ; printf("\n") ;

  printf("\n=================================================\n") ;
  for(i=0 ; i<NPTS ; i++) { zi[i] = .00001 + .012345 * (i - NPTS/3 ) ;  zo[i] = -99999.0 ; }
  avgi = 0.0 ;
  minval = maxval = zi[0] ;
  for(i=0 ; i<NPTS ; i++) { 
    maxval = (zi[i] > maxval) ? zi[i] : maxval ; 
    minval = (zi[i] < minval) ? zi[i] : minval ;
    avgi += zi[i] ;
  }
  avgi /= NPTS ;
  printf("npts = %d, minval = %g, maxval = %g, range = %g, average = %g", NPTS, minval, maxval, maxval - minval, avgi);

  quantum = 1.0f / 8.0f ;
  printf(", quantum = %g", quantum) ;

  float_quantize_prep(NBITS, &p, maxval, minval, quantum);
  printf(", nbits = %d\n", p.nbits) ;
  TIME_ONCE_EZ(1, float_quantize_prep(NBITS, &p, maxval, minval, quantum) )
  printf("float_quantize_prep : %s\n",timer_msg);

  TIME_ONCE_TOP ;
  float_quantize_prep(NBITS, &p, maxval, minval, quantum) ;
  float_quantize_prep(NBITS, &p, maxval, minval, quantum) ;
  float_quantize_prep(NBITS, &p, maxval, minval, quantum) ;
  float_quantize_prep(NBITS, &p, maxval, minval, quantum) ;
  float_quantize_prep(NBITS, &p, maxval, minval, quantum) ;
  TIME_ONCE_BOT_EZ(5) ;
  printf("float_quantize_prep : %s\n",timer_msg);

  float_quantize(iw, zi, NPTS/2, NPTS/2, NPTS/2, 2, &p ) ;
  TIME_LOOP_EZ(nloops, NPTS, float_quantize(iw, zi, NPTS/2, NPTS/2, NPTS/2, 2, &p ) )
  printf("float_quantize      : %s\n",timer_msg);

  TIME_LOOP_TOP(nloops)
  float_quantize(iw, zi, NPTS/2, NPTS/2, NPTS/2, 2, &p ) ;
  TIME_LOOP_BOT_EZ(NPTS)
  printf("float_quantize      : %s\n",timer_msg);

  float_unquantize(iw, zo, NPTS/2, NPTS/2, NPTS/2, 2, &p);
  TIME_LOOP_EZ(nloops, NPTS, float_unquantize(iw, zo, NPTS/2, NPTS/2, NPTS/2, 2, &p) )
  printf("float_unquantize    : %s\n",timer_msg);

  TIME_LOOP_TOP(nloops)
  float_unquantize(iw, zo, NPTS/2, NPTS/2, NPTS/2, 2, &p) ;
  TIME_LOOP_BOT_EZ(NPTS)
  printf("float_unquantize    : %s\n",timer_msg);

  printf("\n=================================================\n") ;
  avgo = 0.0 ;
  for(i=0 ; i<NPTS ; i++) avgo += zi[i] ;
  avgo /= NPTS ;

  for(j = 2 ; j < 25 ; j+=2) {
    NBITS = j ;
    i = 1 << (NBITS);
    toler = zi[NPTS-1]; toler -= zi[0] ;
    toler /= i; //  toler *= 2;
    if(j == 24) toler *= 2 ;

    quantum = toler * 2.0 ;
    float_quantize_prep(NBITS, &p, maxval, minval, quantum);
    if(p.e > -127 && p.e < 127){
      float_quantize(iw, zi, NPTS/2, NPTS/2, NPTS/2, 2, &p ) ;
      float_unquantize(iw, zo, NPTS/2, NPTS/2, NPTS/2, 2, &p);
    }else{
      printf(" [64] ");
      float_quantize(iw, zi, NPTS/2, NPTS/2, NPTS/2, 2, &p ) ;
      float_unquantize(iw, zo, NPTS/2, NPTS/2, NPTS/2, 2, &p);
    }
    printf("%2d bits: ", NBITS);
    avgo = 0.0 ;
    for(i=0 ; i<NPTS ; i++) avgo += zo[i] ;
    avgo /= NPTS ;
    printf("avgo(%9.4g), ", avgo) ;
    error = 0 ; maxdelta = 0.0 ; avgdelta = 0.0 ; avg = 0.0 ;
    for(i=0 ; i<NPTS ; i++) {
      delta = zo[i] ; delta -= zi[i] ;
      avg += delta ;
      delta = ABS(delta) ;
      avgdelta += delta ;
      maxdelta = (delta > maxdelta) ? delta : maxdelta ;
      if(delta > toler) error++;
    }
    avg /= NPTS;
    printf("maxd= %9.4g (%9.4g), ", maxdelta, avgdelta/NPTS);
    printf("bias= %9.3g, toler= %8.3g, bias/avg = %8.2g, maxd/t = %6.4f, exceed=%d\n",
          avg, toler, avg/avgi, maxdelta/toler, error);
  }

  IntPair int_extrema ;
  FloatPair float_extrema ;
  int nbits ;
  float factor ;

  printf("\n===============       simple quantizer       =================\n");
  printf("=============== rounding test with 19 points =================\n");
  quantum = quantum_adjust(.35f) ;
  printf("eps = %15.9f, eps_plus = %15.9f, eps_minus = %15.9f, quantum = %15.9f\n", 
         epsilon, eps_plus, eps_minus, quantum) ;

  printf("quantizing/restoring (9 - i) * .125f\n");
  for(i=0 ; i<19 ; i++) zi[i] = (9 - i) * .125f ;
  for(i=0 ; i<19 ; i++) printf("%8.3f ", zi[i]) ; printf("\n") ;
  nbits   = float_quantize_simple_1D(zi, iw, 19, quantum, &int_extrema) ;
  float_unquantize_simple_1D(zo, iw, 19, quantum, &float_extrema) ;
  for(i=0 ; i<19 ; i++) printf("%8.3f ", zo[i]) ; printf("\n") ;
  printf("float_extrema = %f %f\n", float_extrema.t[0], float_extrema.t[1]) ;
  for(i=0 ; i<19 ; i++) printf("%8d ", iw[i]) ; printf("\n") ;
  printf("int_extrema = %d (%d bits) %d (%d bits), nbits = %d\n",
         int_extrema.t[0], NeedBits(int_extrema.t[0]), int_extrema.t[1], NeedBits(int_extrema.t[1]), nbits) ;

  printf("quantizing/restoring (i - 9) * .125f * (1 - epsilon)\n");
  for(i=0 ; i<19 ; i++){ zi[i] = (i - 9) * .125f ; zi[i] *= eps_minus ; }
  nbits = float_quantize_simple_1D(zi, iw, 19, quantum, &int_extrema) ;
  float_unquantize_simple_1D(zo, iw, 19, quantum, &float_extrema) ;
  for(i=0 ; i<19 ; i++) printf("%8.3f ", zi[i]) ; printf("\n") ;
  for(i=0 ; i<19 ; i++) printf("%8.3f ", zo[i]) ; printf("\n") ;
  for(i=0 ; i<19 ; i++) printf("%8d ", iw[i]) ; printf("\n") ;
  printf("int_extrema = %d (%d bits) %d (%d bits), nbits = %d\n",
         int_extrema.t[0], NeedBits(int_extrema.t[0]), int_extrema.t[1], NeedBits(int_extrema.t[1]), nbits) ;

  printf("\n=================================================\n") ;

  for(factor = 4.0f ; factor <= 16384.0f ; factor *= 8.0f ){
    for(i=0 ; i<NPTS ; i++) { zi[i] = .00001 + .012345 * (i - NPTS/3 ) ;  zo[i] = -99999.0 ; }
    for(i=0 ; i<NPTS ; i++) zi[i] *= (32767.0f/(NPTS-1)) ;
    avgi = 0.0 ;
    minval = maxval = zi[0] ;
    for(i=0 ; i<NPTS ; i++) { 
      maxval = (zi[i] > maxval) ? zi[i] : maxval ; 
      minval = (zi[i] < minval) ? zi[i] : minval ;
      avgi += zi[i] ;
    }
    avgi /= NPTS ;
    quantum = quantum_adjust(1.00001f / factor) ;
    printf("npts = %dx64, minval = %g, maxval = %g, range = %g, average = %g", NPTS/64, minval, maxval, maxval - minval, avgi);
    printf(", quantum = %g\n", quantum) ;

    nbits = float_quantize_simple_1D(zi, iw, NPTS, quantum, &int_extrema ) ;
    printf("range = %d, min = %d, max = %d",
          int_extrema.t[1] - int_extrema.t[0], int_extrema.t[0], int_extrema.t[1]) ;
    printf(", int_extrema = %d (%d bits) %d (%d bits), nbits = %d\n",
          int_extrema.t[0], NeedBits(int_extrema.t[0]), int_extrema.t[1], NeedBits(int_extrema.t[1]), nbits) ;
    TIME_LOOP_EZ(nloops, NPTS, nbits = float_quantize_simple_1D(zi, iw, NPTS, quantum, &int_extrema ) ) ;
    printf("float_quantize_simple_1D   : %s\n",timer_msg);
    TIME_LOOP_EZ(nloops, NPTS, nbits = float_quantize_simple(zi, iw, NPTS/64, NPTS/64, NPTS/64, 64, quantum, &int_extrema ) ) ;
    printf("float_quantize_simple      : %s\n",timer_msg);

    float_unquantize_simple_1D(zo, iw, NPTS, quantum, &float_extrema) ;
    TIME_LOOP_EZ(nloops, NPTS, float_unquantize_simple_1D(zo, iw, NPTS, quantum, &float_extrema) ) ;
    printf("float_unquantize_simple_1D : %s\n",timer_msg);
    TIME_LOOP_EZ(nloops, NPTS, float_unquantize_simple(zo, iw, NPTS/64, NPTS/64, NPTS/64, 64, quantum, &float_extrema) ) ;
    printf("float_unquantize_simple    : %s\n",timer_msg);

    for(i=0 ; i<NPTS ; i++) { zo[i] = -999999999.0 ; iw[i] = -999999999 ; }
    nbits = float_quantize_simple(zi, iw, NPTS/64, NPTS/64, NPTS/64, 64, quantum, &int_extrema ) ;
    float_unquantize_simple(zo, iw, NPTS/64, NPTS/64, NPTS/64, 64, quantum, &float_extrema) ;
    printf("nbits = %d, ", nbits);

    toler = quantum * .5 ;
    avgo = 0.0 ;
    for(i=0 ; i<NPTS ; i++) avgo += zo[i] ;
    avgo /= NPTS ;
    printf("avgo(%9.4g), ", avgo) ;
    error = 0 ; maxdelta = 0.0 ; avgdelta = 0.0 ; avg = 0.0 ;
    for(i=0 ; i<NPTS ; i++) {
      delta = zo[i] ; delta -= zi[i] ;
      avg += delta ;
      delta = ABS(delta) ;
      avgdelta += delta ;
      maxdelta = (delta > maxdelta) ? delta : maxdelta ;
      if(delta > toler) error++;
    }
    avg /= NPTS;
    printf("maxd= %9.4g (%9.4g), ", maxdelta, avgdelta/NPTS);
    printf("bias= %9.3g, toler= %9.4g, bias/avg = %7.2g, maxd/t = %4.2f, exceed=%d\n",
          avg, toler, avg/avgi, maxdelta/toler, error);
    printf("\n") ;
  }
}
