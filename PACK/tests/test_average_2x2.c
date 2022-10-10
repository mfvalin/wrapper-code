#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <misc_timers.h>
#include <misc_operators.h>
#include <average_2x2.h>

#define NTIMES 10000

int main(int argc, char **argv){
  int ni, nj, i, j, ij, nia, nja ;
  int32_t *isrc, *iavg, *idst, *idst2, *ires ;
  float *fsrc, *favg, *fdst, *fdst2, *fres ;
  uint64_t tmin, tmax, freq ;
  double nano, tavg ;
  int niter = NTIMES ;
  char buf[1024] ;
  size_t bufsiz = sizeof(buf) ;
  int errors ;

  freq = cycles_counter_freq() ;
  nano = 1000000000.0 ;
  printf("nano = %8.2G\n", nano) ;
  nano = nano / freq ;
  printf("nano = %8.2G\n", nano) ;

  if(argc<3) {
    printf("usage : %s ni nj\n",argv[0]) ;
    return 1 ;
  }
  printf("SOURCE :") ; for(i=-12 ; i<13 ; i++) printf("%4d", i       ) ; printf("\n") ;
  printf("IDIV2R :") ; for(i=-12 ; i<13 ; i++) printf("%4d", IDIV2R(i)) ; printf("\n") ;
  printf("IDIV4R :") ; for(i=-12 ; i<13 ; i++) printf("%4d", IDIV4R(i)) ; printf("\n") ;
  printf("IDIV8R :") ; for(i=-12 ; i<13 ; i++) printf("%4d", IDIV8R(i)) ; printf("\n") ;
  printf("================================================================\n") ;
  ni = atoi(argv[1]) ; nia = (ni+1) / 2 ;
  nj = atoi(argv[2]) ; nja = (nj+1) / 2 ;
  printf("ni,nja = (%d,%d), nj,nja = (%d,%d)\n", ni, nia, nj, nja) ;
  isrc = (int32_t *) malloc(sizeof(int32_t) * ni  * nj ) ;
  fsrc = (float *) malloc(sizeof(float) * ni  * nj ) ;
  ires = (int32_t *) malloc(sizeof(int32_t) * ni  * nj ) ;
  fres = (float *) malloc(sizeof(float) * ni  * nj ) ;
  idst = (int32_t *) malloc(sizeof(int32_t) * ni  * nj ) ;
  fdst = (float *) malloc(sizeof(float) * ni  * nj ) ;
  idst2 = (int32_t *) malloc(sizeof(int32_t) * ni * 2 ) ;
  fdst2 = (float *) malloc(sizeof(float) * ni * 2 ) ;
  iavg = (int32_t *) malloc(sizeof(int32_t) * nia * nja) ;
  favg = (float *) malloc(sizeof(float) * nia * nja) ;
  for(j=0 ; j<nj ; j++){
    ij = j * ni ;
    for(i=0 ; i<ni ; i++){
//       isrc[ij+i] = ij+i ;
      isrc[ij+i] = j+i ;
      isrc[ij+i] = 4 * (isrc[ij+i] - (ni+nj)/2 + 1) ;
      fsrc[ij+i] = isrc[ij+i] ;
//       isrc[ij+i] = (-isrc[i]) ;
      idst[ij+i] = 999 ; fdst[ij+i] = 999.0f ;
      idst2[i]   = 888 ; fdst2[i]   = 888.0f ;
      ires[ij+i] = 777 ; fres[ij+i] = 777.0f ;
    }
  }
  if(ni*nj<400)
  for(j=nj-1 ; j>=0 ; j--){
    printf("%3d -", j);
    ij = j * ni ;
    for(i=0 ; i<ni ; i++)printf("%4d", isrc[ij+i]) ;
    printf("\n") ;
  }
  printf("================================================================\n") ;

  if(ni*nj<400)
    for(i=0 ; i<ni ; i++)printf("%4d", isrc[i]) ; printf("\n") ;
  for(i=0 ; i<nia*nja ; i++) iavg[i] = 666 ;
  average_2x1_I32(isrc, iavg, ni);
  if(ni*nj<400)
    for(i=0 ; i<nia+1 ; i++)printf("%8d", iavg[i]) ; printf("\n") ;
//   restore_2x1(idst2, iavg, ni) ;
  for(i=0 ; i<ni*2 ; i++) idst2[i] = 666 ;
  expand_2x1_row_along_x_I32(idst2, iavg, ni) ;
  if(ni*nj<400)
    for(i=0 ; i<ni+1 ; i++)printf("%4d", idst2[i]) ; printf("\n") ;
  errors = 0 ;
  for(i=0 ; i<ni ; i++) if(isrc[i] != idst2[i]) errors++ ;
    printf("expand_2x1_row_along_x_I32    : errors = %d\n", errors);
  TIME_LOOP(tmin, tmax, tavg, NTIMES, (ni), buf, bufsiz, expand_2x1_row_along_x_I32(idst2, iavg, ni) )
  printf("expand_2x1_row_along_x_I32    : %s\n",buf);
//   TIME_LOOP(tmin, tmax, tavg, NTIMES, (1), buf, bufsiz, i=0 )
//   printf("overhead    : %s\n",buf);
  printf("================================================================\n") ;

  if(ni*nj<400)
    for(i=0 ; i<ni ; i++)printf("%4.0f", fsrc[i]) ; printf("\n") ;
  average_2x1_F32(fsrc, favg, ni);
  if(ni*nj<400)
    for(i=0 ; i<nia ; i++)printf("%8.0f", favg[i]) ; printf("\n") ;
  expand_2x1_row_along_x_F32(fdst2, favg, ni) ;
  if(ni*nj<400)
    for(i=0 ; i<ni ; i++)printf("%4.0f", fdst2[i]) ; printf("\n") ;
  errors = 0 ;
  for(i=0 ; i<ni ; i++) if(fsrc[i] != fdst2[i]) errors++ ;
    printf("expand_2x1_row_along_x_F32    : errors = %d\n", errors);
  TIME_LOOP(tmin, tmax, tavg, NTIMES, (ni), buf, bufsiz, expand_2x1_row_along_x_F32(fdst2, favg, ni) )
  printf("expand_2x1_row_along_x_F32    : %s\n",buf);
  printf("================================================================\n") ;

  average_2x2_2D_I32(isrc, iavg, ni, ni, nj) ;
  if(ni*nj<400)
  for(j=nja-1 ; j>=0 ; j--){
    printf("%3d -", j);
    ij = j * nia ;
    for(i=0 ; i<nia ; i++)printf("%4d", iavg[ij+i]) ;
    printf("\n") ;
  }
  TIME_LOOP(tmin, tmax, tavg, NTIMES, (ni*nj), buf, bufsiz, average_2x2_2D_I32(isrc, iavg, ni, ni, nj) )
  printf("average_2x2_2D_I32    : %s\n",buf);
  printf("================================================================\n") ;
  average_2x2_2D_F32(fsrc, favg, ni, ni, nj) ;
  if(ni*nj<400)
  for(j=nja-1 ; j>=0 ; j--){
    printf("%3d -", j);
    ij = j * nia ;
    for(i=0 ; i<nia ; i++)printf("%4.0f", favg[ij+i]) ;
    printf("\n") ;
  }
  TIME_LOOP(tmin, tmax, tavg, NTIMES, (ni*nj), buf, bufsiz, average_2x2_2D_F32(fsrc, favg, ni, ni, nj) )
  printf("average_2x2_2D_F32    : %s\n",buf);
  printf("================================================================\n") ;
  TIME_LOOP(tmin, tmax, tavg, NTIMES, (ni*nj), buf, bufsiz, avgres_2x2_2D_I32(isrc, iavg, ires, ni, ni, nj) )
  printf("avgres_2x2_2D_I32    : %s\n",buf);
  printf("================================================================\n") ;
  expand_2x2_2D_I32(idst, iavg, ni, ni, nj) ;
  if(ni*nj<400)
  for(j=nj-1 ; j>=0 ; j--){
    printf("%3d -", j);
    ij = j * ni ;
    for(i=0 ; i<ni ; i++)printf("%4d", idst[ij+i]) ;
    printf("\n") ;
  }
  TIME_LOOP(tmin, tmax, tavg, NTIMES, (ni*nj), buf, bufsiz, expand_2x2_2D_I32(idst, iavg, ni, ni, nj) )
  errors = 0 ;
  for(i=0 ; i<ni*nj ; i++) if(isrc[i] != idst[i]) errors++ ;
  printf("expand_2x2_2D_I32    : %s, errors = %d\n", buf, errors);
  printf("================================================================\n") ;
  expand_2x2_2D_F32(fdst, favg, ni, ni, nj) ;
  if(ni*nj<400)
  for(j=nj-1 ; j>=0 ; j--){
    printf("%3d -", j);
    ij = j * ni ;
    for(i=0 ; i<ni ; i++)printf("%4.0f", fdst[ij+i]) ;
    printf("\n") ;
  }
  TIME_LOOP(tmin, tmax, tavg, NTIMES, (ni*nj), buf, bufsiz, expand_2x2_2D_F32(fdst, favg, ni, ni, nj) )
  errors = 0 ;
  for(i=0 ; i<ni*nj ; i++) if(fsrc[i] != fdst[i]) errors++ ;
  printf("expand_2x2_2D_F32    : %s, errors = %d\n", buf, errors);
}
