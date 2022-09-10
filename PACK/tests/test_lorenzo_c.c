#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>

#include <misc_timers.h>

#if !defined(NPTS)
#define NPTS 64
#endif

#define NTIMES 100000

void LorenzoPredict2D(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);
void LorenzoUnpredict2D(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);
void LorenzoPredict2D_IJ(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);
void LorenzoPredictShort(int32_t * restrict orig, int32_t * restrict diff, int ni, int lnio, int lnid, int nj);

int main(int argc, char **argv){
  int32_t data[NPTS+1][NPTS] ;
  int32_t data2[NPTS+1][NPTS] ;
  int32_t pred[NPTS+1][NPTS] ;
  int32_t pred2[NPTS+1][NPTS] ;
  int i, j, errors ;
  uint64_t tmin, tmax, freq ;
  double nano, tavg ;
  int niter = NTIMES ;
  char buf[1024] ;
  size_t bufsiz = sizeof(buf) ;

  freq = cycles_counter_freq() ;
  nano = 1000000000.0 ;
  printf("nano = %8.2G\n", nano) ;
  nano = nano / freq ;
  printf("nano = %8.2G\n", nano) ;

  for(j=0 ; j<NPTS+1 ; j++){
    for(i=0 ; i<NPTS ; i++){
      data[j][i]  = (2*i + 3*j + 5) ;
      data2[j][i] = 999999 ;
      pred[j][i]  = 999999 ;
      pred2[j][i] = 999999 ;
    }
  }

  LorenzoPredict2D(&data[0][0], &pred[0][0], NPTS, NPTS, NPTS, NPTS) ;
  for(i=0 ; i<8 ; i++) printf("%d ",pred[NPTS][i]); printf("\n");
  LorenzoUnpredict2D(&data2[0][0], &pred[0][0], NPTS, NPTS, NPTS, NPTS) ;
  for(i=0 ; i<8 ; i++) printf("%d ",data2[NPTS][i]); printf("\n");
  errors = 0 ;
  for(j=0 ; j<NPTS ; j++){
    for(i=0 ; i<NPTS ; i++){
      if(data2[j][i] != (2*i + 3*j + 5) ) errors++;
    }
  }
  printf("LorenzoPredict2D    : errors = %d\n\n",errors);

  LorenzoPredictShort(&data[0][0], &pred2[0][0], NPTS, NPTS, NPTS, NPTS) ;
  for(i=0 ; i<8 ; i++) printf("%d ",pred2[NPTS][i]); printf("\n");
  for(j=0 ; j<NPTS+1 ; j++) for(i=0 ; i<NPTS ; i++) data2[j][i] = 999999 ;
  LorenzoUnpredict2D(&data2[0][0], &pred2[0][0], NPTS, NPTS, NPTS, NPTS) ;
  for(i=0 ; i<8 ; i++) printf("%d ",data2[NPTS][i]); printf("\n");
  errors = 0 ;
  for(j=0 ; j<NPTS ; j++){
    for(i=0 ; i<NPTS ; i++){
      if(data2[j][i] != (2*i + 3*j + 5) ) errors++;
    }
  }
  printf("LorenzoPredict2D_S  : errors = %d\n\n",errors);

  TIME_LOOP(tmin, tmax, tavg, NTIMES, (NPTS*NPTS), buf, bufsiz, LorenzoPredict2D(&data[0][0], &pred[0][0], NPTS, NPTS, NPTS, NPTS) )
  printf("LorenzoPredict2D    : %s\n",buf);

  TIME_LOOP(tmin, tmax, tavg, NTIMES, (NPTS*NPTS), buf, bufsiz, LorenzoPredictShort(&data[0][0], &pred[0][0], NPTS, NPTS, NPTS, NPTS) )
  printf("LorenzoPredict2D_S  : %s\n",buf);

  TIME_LOOP(tmin, tmax, tavg, NTIMES, (NPTS*NPTS), buf, bufsiz, LorenzoUnpredict2D(&data[0][0], &pred[0][0], NPTS, NPTS, NPTS, NPTS) )
  printf("LorenzoUnpredict2D  : %s\n",buf);
}
