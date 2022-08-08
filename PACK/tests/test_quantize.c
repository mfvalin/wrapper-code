#include <stdio.h>

#include <misc_pack.h>
#include <misc_operators.h>

#define NPTS 32800

int main(){
  float zi[NPTS] ;
  float zo[NPTS] ;
  int iw[NPTS] ;
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

  for(i=0 ; i<NPTS ; i++) { zi[i] = .00001 + .012345 * (i - NPTS/2 ) ;  zo[i] = -99999.0 ; }
  avgi = 0.0 ;
  for(i=0 ; i<NPTS ; i++) { zi[i] *= 1.6E+30 ; avgi += zi[i] ;}
  avgi /= NPTS ;
  printf("avgi = %15.8g\n\n", avgi) ;
  minval = maxval = zi[0] ;
  for(i=1 ; i<NPTS ; i++) { 
    maxval = (zi[i] > maxval) ? zi[i] : maxval ; minval = (zi[i] < minval) ? zi[i] : minval ;
  }

  for(j = 1 ; j < 25 ; j++) {
    NBITS = j ;
    i = 1 << (NBITS);
    toler = zi[NPTS-1]; toler -= zi[0] ;
    toler /= i; //  toler *= 2;

    quantum = toler * 2.0 ;
    float_quantize_prep(NBITS, &p, maxval, minval, quantum);
//     printf("offset = %10d, exp = %d, nbits = %d, \n", p.o, p.e, NBITS);
    if(p.e > -127 && p.e < 127){
      printf(" [32] ");
      float_quantize(iw, zi, NPTS/2, NPTS/2, NPTS/2, 2, &p ) ;
      float_unquantize(iw, zo, NPTS/2, NPTS/2, NPTS/2, 2, &p);
    }else{
      printf(" [64] ");
      float_quantize(iw, zi, NPTS/2, NPTS/2, NPTS/2, 2, &p ) ;
      float_unquantize(iw, zo, NPTS/2, NPTS/2, NPTS/2, 2, &p);
    }
    printf("nbits = %2d, max, min, rng = %9.5g %9.5g %9.5g, ", NBITS, maxval, minval, maxval-minval);
    avgo = 0.0 ;
    for(i=0 ; i<NPTS ; i++) avgo += zo[i] ;
    avgo /= NPTS ;
    printf("avgo = %15.8g, ", avgo) ;
    error = 0 ; maxdelta = 0.0 ; avgdelta = 0.0 ; avg = 0.0 ;
    for(i=0 ; i<NPTS ; i++) {
//       delta = ABS(zo[i]-zi[i]) ;
      delta = zo[i] ; delta -= zi[i] ;
      avg += delta ;
      delta = ABS(delta) ;
      avgdelta += delta ;
      maxdelta = (delta > maxdelta) ? delta : maxdelta ;
      if(delta > toler) error++;
    }
    avg /= NPTS;
//     printf("%d points from %15.5f to %15.5f, maxdelta = %15.8f, avgdelta = %15.8f\n",
//            NPTS, zi[0], zi[NPTS-1], maxdelta, avgdelta/NPTS);
    printf("maxdelta = %10.4g, avgdelta = %10.4g, ",
           maxdelta, avgdelta/NPTS);
    printf("bias = %15.8g, toler = %9.5g, bias/toler = %8.4f, maxdelta/toler = %6.4f, exceed=%d\n\n",
          avg, toler, avg/toler, maxdelta/toler, error);
  }
}
