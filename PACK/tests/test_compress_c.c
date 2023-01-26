#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

#include <idwt53.h>
#include <rmn/misc_pack.h>
#include <rmn/lorenzo.h>
#include <rmn/misc_operators.h>

static inline int indx3(int i, int j, int k, int ni, int nj,  int nk){
  return (i + ni*j + nj*k) ;
}

void *read_raw_data(char *filename, int fdi, int *ni, int *nj, int *nk, int *ndim);

void print_8x8_i(int32_t *what, int32_t lni){
  int i, j ;
  for(j=0 ; j<8 ; j++){
    for(i=0 ; i<8 ; i++) fprintf(stderr, "%8d", what[i]) ;
    fprintf(stderr, "\n") ;
    what += lni ;
  }
}

void to_zigzag_8x8(int32_t *what, int32_t lni){
  int i, j ;
  uint32_t *zigzag ;
  for(j=0 ; j<8 ; j++){
    zigzag = (uint32_t *) what;
    for(i=0 ; i<8 ; i++) zigzag[i] = to_zigzag_32(what[i]) ;
    what += lni ;
  }
}

void from_zigzag_8x8(int32_t *what, int32_t lni){
  int i, j ;
  uint32_t *zigzag ;
  for(j=0 ; j<8 ; j++){
    zigzag = (uint32_t *) what;
    for(i=0 ; i<8 ; i++) what[i] = from_zigzag_32(zigzag[i]) ;
    what += lni ;
  }
}

int compare_integers(int *f1, int *f2, int n){
  int errors = 0 ;
  int i ;
  for(i=0 ; i<n ; i++) errors += ( (f1[i] != f2[i]) ? 1 : 0 ) ;
  return errors ;
}

float bias_err(float *f1, float *f2, int n){
  double sumerr = 0.0 ;
  float bias ;
  int i ;
  for(i=0 ; i<n ; i++) sumerr += (f1[i] - f2[i]) ;
  return bias = sumerr / n ;
}

float max_abs_err(float *f1, float *f2, int n){
  float errabs, maxerr = 0.0 ;
  int i ;
  for(i=0 ; i<n ; i++){
    errabs = f1[i] - f2[i] ;
    errabs = (errabs <0) ? -errabs : errabs ;
    maxerr = (errabs > maxerr) ? errabs : maxerr ;
  }
  return maxerr ;
}

static int blocks = 0 ;
static int btab[33] ;

int code_8x8_u32(void *iwhat, int ni){
  uint32_t *what = (uint32_t *) iwhat ;
  int nbits = 0 ;
  int nbmax = 0 ;
  int i, j ;
  for(j=0 ; j<8 ; j++){
    for(i=0 ; i<8 ; i++){
      nbits = BitsNeeded_u32((uint32_t) what[i]) ;
      nbmax = (nbits > nbmax) ? nbits : nbmax ;
    }
    what += ni ;
  }
  blocks++ ;
  btab[nbmax]++ ;
// if( (blocks&0xFFF) == 0) fprintf(stderr,"DEBUG, nbmax = %d\n", nbmax) ;
  return nbmax * 64 ;
}

int main(int argc, char **argv){
  int i, j, ni, nj, nk, ndim, fd, nerrors, gi, gj, bki, bkj ;
  float *data = NULL ;
  float *d = NULL ;
  int32_t *q = NULL ;
  int32_t *ql = NULL ;
  int32_t *qw = NULL ;
  float quantum = 0.0f ;
  IntPair t ;
  FloatPair tf ;
  char vname ;
  float totbits ;

  if(argc < 2){
    fprintf(stderr,"USAGE: %s file_name\n", argv[0]) ;
    exit(1) ;
  }
  ndim = 2 ; ni = nj = nk = 1 ;
  if(argc >2) ndim = atoi(argv[2]) ;

  data = read_raw_data(argv[1], fd = 0, &ni, &nj, &nk, &ndim) ;
  vname = argv[1][4] ;
  if(data == NULL){
    fprintf(stderr,"ERROR reading %dD %s\n", ndim, argv[1]) ;
  }else{
    bki = ni/8 ;
    bkj = nj/8 ;
    fprintf(stderr,"INFO: ndim = %d, ni = %d, nj = %d, nk = %d, , %x x %d blocks, data[0,last] = %11f %11f (%s)\n", 
            ndim, ni, nj, nk, bki, bkj, data[0], data[ni*nj*nk-1], argv[1]) ;
  }

  d  = malloc(ni*nj*nk*sizeof(float)) ;
  q  = malloc(ni*nj*nk*sizeof(int32_t)) ;
  ql = malloc(ni*nj*nk*sizeof(int32_t)) ;
  qw = malloc(ni*nj*nk*sizeof(int32_t)) ;
//  uint32_t float_quantize_simple(float * restrict z, int32_t * restrict q, int ni, int lniz, int lniq, int nj, float quantum, IntPair *t);
  if(vname == 'T') quantum = 0.001 ; // degree C
  if(vname == 'E') quantum = 0.001 ;
  if(vname == 'G') quantum = 0.05  ; // dam
  if(vname == 'Z') quantum = 0.01  ; // m
  if(vname == 'U') quantum = 0.01  ;
  if(vname == 'V') quantum = 0.01  ;
  if(vname == 'W') quantum = 0.001 ;

  int base = indx3(128, 256, 1, ni, nj, nk) ;
  base = 0 ;
  base = indx3(112, 955, 1, ni, nj, nk) ;

  // quantize and check unquantize errors, expected error is ~ quantum/2
  uint32_t nbits = float_quantize_simple(data, q, ni, ni, ni, nj, quantum, &t) ;
  fprintf(stderr,"INFO: quantized %c using at most %d bits, quantum = %f, t = %d %d\n", vname, nbits, quantum, t.t[0], t.t[1]) ;
  print_8x8_i(q+base, ni) ;
  float_unquantize_simple(d, q, ni, ni, ni, nj, quantum, &tf) ;
  fprintf(stderr,"INFO: restored %c, tf = %f %f, maxerr = %f, bias/quantum = %8.1e\n", 
          vname, tf.t[0], tf.t[1], max_abs_err(data, d, ni*nj), bias_err(data, d, ni*nj)/quantum) ;

  // check Lorenzo predict/unpredict errors, 0 is expected
  fprintf(stderr,"INFO: Lorenzo predict/unprecict\n") ;
  LorenzoPredict_c(q, ql, ni, ni, ni, nj);
//   print_8x8_i(ql+base, ni) ;
  fprintf(stderr, "zigzag\n") ;
  to_zigzag_8x8(ql+base, ni) ;
  print_8x8_i(ql+base, ni) ;
  fprintf(stderr, "\n") ;
  from_zigzag_8x8(ql+base, ni) ;
  LorenzoUnpredictInplace_c(ql, ni, ni, nj);
  print_8x8_i(ql+base, ni) ;
  nerrors = compare_integers(q, ql, ni*nj) ;
  fprintf(stderr, "errors = %d\n", nerrors) ;

  LorenzoPredict_c(q, ql, ni, ni, ni, nj);
  int base0 ;
  nbits = 0 ;
  for(gj=0 ; gj < bkj ; gj++){
    for(gi=0 ; gi < bki ; gi++){
      base0 = indx3(8*gi, 8*gj, 1, ni, nj, nk) ;
      to_zigzag_8x8(ql+base0, ni) ;
      nbits += code_8x8_u32(ql+base0, ni) ;
    }
  }
  fprintf(stderr,"INFO: btab =") ; for(i=0 ; i<33 ; i++) fprintf(stderr," %d", btab[i]) ; fprintf(stderr,"\n") ;
  totbits = nbits ;
  totbits = totbits/(ni*nj) ;
  fprintf(stderr,"INFO: nbits = %d (%5.2f/point), blocks = %d\n", nbits, totbits, blocks) ;

  // check 53 wavelet transform/restore and errors, 0 is expected
  fprintf(stderr,"INFO: DWT52 forward/inverse\n") ;
  FDWT53i_2D_split_inplace_n(q+base, 8, ni, 8, 3) ;
//   print_8x8_i(q+base, ni) ;
  fprintf(stderr, "zigzag\n") ;
  to_zigzag_8x8(q+base, ni) ;
  print_8x8_i(q+base, ni) ;
  fprintf(stderr, "\n") ;
  from_zigzag_8x8(q+base, ni) ;

  IDWT53i_2D_split_inplace_n(q+base, 8, ni, 8, 3) ;
  print_8x8_i(q+base, ni) ;
  nerrors = compare_integers(q, ql, ni*nj) ;
  fprintf(stderr, "errors = %d\n", nerrors) ;
}

