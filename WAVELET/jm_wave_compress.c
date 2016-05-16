#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* These scale factors should give the HF basis functions roughly the same
   energy as the LF basis functions. They are SNR-optimal for an orthogonal
   transform. Since the 5/3 wavelet is biorthogonal, they should probably be
   tuned a bit. */
#define LF_SCALE 1.41421356237f
#define HF_SCALE 0.7071067812f

#define LF_SCALE_1 (1./LF_SCALE);
#define HF_SCALE_1 (1./HF_SCALE);

/* Number of transform levels */
#define LEVELS 5

/* Predict odd samples from even samples (and subtract prediction). */
static void predict(float *x, int N) {
  int i;
  for (i=1;i<N-1;i+=2) {
    x[i] -= .5*(x[i-1] + x[i+1]);
  }
  if ((N&1)==0) x[i] -= x[i-1];
}

/* Update even samples from odd samples so that the even samples
   end up being a low-pass-filtered version of the signal. */
static void update(float *x, int N) {
  int i;
  x[0] += .5*x[1];
  for (i=2; i < N-1; i+=2) {
    x[i] += .25*(x[i-1] + x[i+1]);
  }
  if (N&1) x[N-1] += .5*x[N-2];
}

/* Undo predict step. */
static void unpredict(float *x, int N) {
  int i;
  for (i=1;i<N-1;i+=2) {
    x[i] += .5*(x[i-1] + x[i+1]);
  }
  if ((N&1)==0) x[i] += x[i-1];
}

/* Undo update step. */
static void unupdate(float *x, int N) {
  int i;
  x[0] -= .5*x[1];
  for (i=2; i < N-1; i+=2) {
    x[i] -= .25*(x[i-1] + x[i+1]);
  }
  if (N&1) x[N-1] -= .5*x[N-2];
}

/* Implements a 5/3 lifting wavelet transform */
static void forward_transform(float *A, int M, int N, int stride) {
  int i, j;
  /* Horizontal transform */
  for (i=0;i<M;i++) {
    float tmp[N];
    for (j=0;j<N;j++) tmp[j] = A[i*stride + j];
    predict(tmp, N);
    update(tmp, N);
    /* Reorder LF first */
    for (j=0;j<(N+1)/2;j++) A[i*stride + j] = tmp[2*j]*LF_SCALE;
    for (j=0;j<N/2;j++) A[i*stride + j + (N+1)/2] = tmp[2*j + 1]*HF_SCALE;
  }
  /* Vertical transform */
  for (j=0;j<N;j++) {
    float tmp[M];
    for (i=0;i<M;i++) tmp[i] = A[i*stride + j];
    predict(tmp, M);
    update(tmp, M);
    /* Reorder LF first */
    for (i=0;i<(M+1)/2;i++) A[i*stride + j] = tmp[2*i]*LF_SCALE;
    for (i=0;i<M/2;i++) A[(i+(M+1)/2)*stride + j] = tmp[2*i + 1]*HF_SCALE;
  }
}

static void inverse_transform(float *A, int M, int N, int stride) {
  int i, j;
  /* Vertical inverse transform */
  for (j=0;j<N;j++) {
    float tmp[M];
    /* Reorder from LF first */
    for (i=0;i<(M+1)/2;i++) tmp[2*i] = A[i*stride + j]*LF_SCALE_1;
    for (i=0;i<M/2;i++) tmp[2*i + 1] = A[(i+(M+1)/2)*stride + j]*HF_SCALE_1;
    unupdate(tmp, M);
    unpredict(tmp, M);
    for (i=0;i<M;i++) A[i*stride + j] = tmp[i];
  }
  /* Horizontal inverse transform */
  for (i=0;i<M;i++) {
    float tmp[N];
    /* Reorder LF first */
    for (j=0;j<(N+1)/2;j++) tmp[2*j] = A[i*stride + j]*LF_SCALE_1;
    for (j=0;j<N/2;j++) tmp[2*j + 1] = A[i*stride + j + (N+1)/2]*HF_SCALE_1;
    unupdate(tmp, N);
    unpredict(tmp, N);
    for (j=0;j<N;j++) A[i*stride + j] = tmp[j];
  }
}

/* inverse integer 2D transform, linear wavelet 5/3, ni and nj MUST be even */
int InverseDwt_53(int *z, int ni, int nj, int np){
  return -1;
}

/* forward integer 2D transform, linear wavelet 5/3, ni and nj MUST be even */
int ForwardDwt_53(int *z, int ni, int nj, int np){
  int todd[ni*nj/2] ;   /* size of uppre half of array z, will contain odd rows */
  int *even, *odd, *zeven, *zodd, *te, *to, *ze, *zo ;
  int i, j, i2, j2 ;
  int ni2, nj2 ;

  /* ======================================================================== */
  /*                                PASS 1 along I                            */
  /* ======================================================================== */
  ni2 = ni >> 1;
  te = todd;
  ze = z;
  for( j=0 ; j<nj ; j+=2 ) {    /* loop over pairs of even+odd rows */
    zo = ze + ni2;
    to = te + ni2;

    /*       even row        */
    for( i=0,i2=0 ; i<ni2 ; i++) { /* split even row into ze/to */
      ze[i] = z[i2++];             /* even points */
      to[i] = z[i2++];             /* odd points  */
    }
    for( i=0 ; i<ni2-1 ; i++) {    /* predict step, odd columns of even row */
      zo[i] = to[i] - ( (ze[i] + ze[i+1]) >> 1 ) ; /* odd[i] = odd[i] - (even[i] + even[i+1]) / 2 */
    }
    zo[i] = to[i] - ze[i];         /* last point for predict step */
    ze[0] = ze[0] + (to[0] >> 1);  /* first point for update step */
    for( i=1 ; i<ni2 ; i++) {      /* update step, even columns of even row */
      ze[i] = ze[i] + ( (to[i-1] + to[i]) >> 2 );  /* even[i] = even[i] + (odd[i-1] + odd[i]) / 4 */
    }

    /*       odd row        */
    for( i=0,i2=0 ; i<ni2 ; i++) { /* split odd row into te/to */
      te[i] = z[i2++];             /* even points */
      to[i] = z[i2++];             /* odd points  */
    }
    for( i=0 ; i<ni2-1 ; i++) {    /* predict step, odd row */
      to[i] = to[i] - ( (te[i] + te[i+1]) >> 1 ) ; /* odd[i] = odd[i] - (even[i] + even[i+1]) / 2 */
    }
    to[i] = to[i] - ze[i];         /* last point for predict step */
    te[0] = te[0] + (to[0] >> 1);  /* first point for update step */
    for( i=1 ; i<ni2 ; i++) {      /* update step, odd row */
      te[i] = te[i] + ( (to[i-1] + to[i]) >> 2 );  /* even[i] = even[i] + (odd[i-1] + odd[i]) / 4 */
    }

    te += ni;      /* ni is storage first dimension of te/todd */
    ze += np;      /* np is storage first dimension of z/ze */
  }
  /* we now have the lower part of z filled with the even rows, todd contains the odd rows */

  /* ======================================================================== */
  /*                                PASS 2 along J                            */
  /* ======================================================================== */
  /* even rows are in zeven, odd rows are in todd */
  nj2 = nj >> 1;
  zeven = z;
  zodd = z + np*nj2;
  te = todd;
  for( j=0 ; j<nj2-1 ; j++ ) {    /* predict step, predict odd rows from even rows */
    for( i=0 ; i<ni ; i++ ) {
      zodd[i] = te[i] - ( (zeven[i] + zeven[i+np]) >> 1 );   /* odd row k = odd row k - (even row k + even row [k+1]) / 2 */
    }
    te += ni;      /* ni is storage first dimension of te/todd */
    zeven += np;   /* np is storage first dimension of z/zeven/zodd */
    zodd += np;
  }
  for( i=0 ; i<ni ; i++ ) {       /* last row of predict step */
    zodd[i] = te[i] - zeven[i] ;
  }
  zeven = z;
  zodd = z + np*nj2;
  for( i=0 ; i<ni ; i++ ) {       /* first row of update step */
    zeven[i] = zeven[i] + ( zodd[i] >> 1);
  }
  for( j=1 ; j<nj2 ; j++ ) {    /* update step */
    zeven += np;
    zodd += np;
    for( i=0 ; i<ni ; i++ ) {
      zeven[i] = zeven[i] + ( (zodd[i-np] + zodd[i])  >> 2);  /* even row k = even row k + (odd row [k-1] + odd row k) / 4 */
    }
  }
  return 0;
}

#if defined (SELF_TEST)
int main() {
  int n, m;
  int i, j;
  float *A, *Q;
  int N, M;
  int level;
  float quant;
  scanf("%d%d", &n, &m);
  /* Pad a a multiple of 1<<LEVELS. While it's easy to apply the transform on
     any size, the tree-based coefficients coding will be a pain if we don't
     have complete trees.*/
  N = (n + (1<<LEVELS) - 1) >> LEVELS << LEVELS;
  M = (m + (1<<LEVELS) - 1) >> LEVELS << LEVELS;
  A = calloc(N*M, sizeof(*A));
  Q = calloc(N*M, sizeof(*A));
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) {
      scanf("%f", &A[i*N + j]);
    }
  }
  /* TODO: Do proper padding on the edges rather than leave zeros. */
  fprintf(stderr, "original size: %d %d\npadded size: %d %d\n", n, m, N, M);
  for (level=0;level<LEVELS;level++)
  {
    forward_transform(A, M>>level, N>>level, N);
  }
  /* This controls the quantization resolution. */
  quant = .0025;
  for (i=0;i<M;i++) {
    for (j=0;j<N;j++) {
      Q[i*N+j] = (int)floor(.5 + A[i*N+j]/quant);
    }
  }
  /* Do the actual encoding here */
  for (i=0;i<M;i++) {
    for (j=0;j<N;j++) {
      A[i*N+j] = Q[i*N+j]*quant;
    }
  }
  /* Apply inverse transform stages in reverse order. */
  for (level=LEVELS-1;level>=0;level--)
  {
    inverse_transform(A, M>>level, N>>level, N);
  }
#if 1
  /* Print reconstructed data. */
  for (i=0;i<m;i++) {
    for (j=0;j<n;j++) {
      printf("%f ", A[i*N + j]);
    }
    printf("\n");
  }
#else
  /* Print quantized coefficients. */
  for (i=0;i<M;i++) {
    for (j=0;j<N;j++) {
      printf("%f ", Q[i*N + j]);
    }
    printf("\n");
  }
#endif
  return 0;
}
#endif
