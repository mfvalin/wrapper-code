#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(CDF97)
#include <cdf97.h>
#else
#include <cdf53.h>
#endif

#if ! defined(NPTS)
#define NPTS 16
#endif

void quadrants(float *z, int n, int ln, float *min, float *max){
  int i, j ;
  min[0] = min[1] = min[2] = min[3] =  999999999.0 ;
  max[0] = max[1] = max[2] = max[3] = -999999999.0 ;
  for (j=0;j<(n+1)/2;j++) {
    for (i=0;i<(n+1)/2;i++) {          // LL quadrant, level 1
      min[0] = (z[i+ln*j] < min[0]) ? z[i+ln*j] : min[0] ;
      max[0] = (z[i+ln*j] > max[0]) ? z[i+ln*j] : max[0] ;
    }
    for (i=(n+1)/2;i<n;i++) {       // HL quadrant, level 1
      min[1] = (z[i+ln*j] < min[1]) ? z[i+ln*j] : min[1] ;
      max[1] = (z[i+ln*j] > max[1]) ? z[i+ln*j] : max[1] ;
    }
  }
  for (j=(n+1)/2;j<n;j++) {
    for (i=0;i<(n+1)/2;i++) {          // LH quadrant, level 1
      min[2] = (z[i+ln*j] < min[2]) ? z[i+ln*j] : min[2] ;
      max[2] = (z[i+ln*j] > max[2]) ? z[i+ln*j] : max[2] ;
    }
    for (i=(n+1)/2;i<n;i++) {       // HH quadrant, level 1
      min[3] = (z[i+ln*j] < min[3]) ? z[i+ln*j] : min[3] ;
      max[3] = (z[i+ln*j] > max[3]) ? z[i+ln*j] : max[3] ;
    }
  }
}

int main(int argc, char **argv) {
  float x[NPTS+1], y[NPTS+1], e[NPTS+1], o[NPTS+1], z[NPTS+1], d[NPTS+1];
  float xy[NPTS][NPTS], xy0[NPTS][NPTS] ;
  int i, j, k;
  double sum2;
  float quantum, quantum2, quantum4, expected ;
  float maxll, maxlh, maxhl, maxhh, minll, minlh, minhl, minhh ;
  float mins[4], maxs[4] ;
  int spanll, spanhl, spanlh, spanhh ;
  int npts2 = (NPTS+1)/2;
  int npts4 = (npts2+1)/2;
#if ! defined(NOQUANTUM)
  int use_quantum = 1 ;
#else
  int use_quantum = 0 ;
#endif
  int levels = 3 ;

  if(argc > 1) use_quantum = atoi(argv[1]);  // default is normally 1
  quantum = quantum2 = quantum4 = .05f ;
  if(argc > 2) quantum = quantum2 = quantum4 = atof(argv[2]) ;
  if(argc > 3) levels = atoi(argv[3]) ;

  for (i=0;i<NPTS;i++) x[i]=5+i+0.4*i*i-0.02*i*i*i;
  for (i=0;i<NPTS;i++) y[i]=x[i];
  for (i=0;i<NPTS;i++) d[i]=x[i];
  for (i=0;i<NPTS;i++) { e[i] = 0 ; o[i] = 0 ; }
  for (j=0;j<NPTS;j++) {
  // Makes a fancy cubic signal
//     for (i=0;i<NPTS;i++) xy[j][i] = (3+i+0.4f*i*i-0.02f*i*i*i) * (3+j+0.4f*j*j-0.02f*j*j*j);
    for (i=0;i<NPTS;i++) {
      xy[j][i] = (3+i+0.4*i*i-0.02*i*i*i) * (3+i+0.4*j*j-0.02*j*j*j);
//       xy[j][i] = 3*i + 4*j ;
//     xy[j][i] = sqrt((i-7.45)*(i-7.45) + (j-7.55)*(j-7.55));
      xy0[j][i] = xy[j][i] ;
    }
  }
  
#if defined(CDF97)
  printf("CDF97 test : Original 2D signal:\n");
#else
  printf("CDF53 test : Original 2D signal:\n");
#endif
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
    printf("\n");
  }
#if defined(CDF97)
  F_CDF97_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, levels);
//   F_CDF97_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
//   F_CDF97_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   F_CDF97_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
#else
  F_CDF53_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, levels);
//   F_CDF53_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
//   F_CDF53_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   F_CDF53_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
#endif

//   printf("quantum used = %8.2f\n",quantum);
//   for (j=0;j<NPTS;j++) {               // quantification pass
//     for (i=0;i<NPTS;i++) { k = xy[j][i] / quantum + .5f ; xy[j][i] = k * quantum ; }
//   }
  minll = minhl = minlh = minhh =  999999999.0 ;
  maxll = maxhl = maxlh = maxhh = -999999999.0 ;
  for (j=0;j<NPTS/2;j++) {
    for (i=0;i<NPTS/2;i++) {          // LL quadrant, level 1
      k = xy[j][i] / quantum + .5f ;
      if(use_quantum) xy[j][i] = k * quantum ;
      minll = (xy[j][i] < minll) ? xy[j][i] : minll ;
      maxll = (xy[j][i] > maxll) ? xy[j][i] : maxll ;
    }
    for (i=NPTS/2;i<NPTS;i++) {       // HL quadrant, level 1
      k = xy[j][i] / quantum2 + .5f ;
      if(use_quantum) xy[j][i] = k * quantum2 ;
      minhl = (xy[j][i] < minhl) ? xy[j][i] : minhl ;
      maxhl = (xy[j][i] > maxhl) ? xy[j][i] : maxhl ;
    }
  }
  for (j=NPTS/2;j<NPTS;j++) {
    for (i=0;i<NPTS/2;i++) {          // LH quadrant, level 1
      k = xy[j][i] / quantum2 + .5f ;
      if(use_quantum) xy[j][i] = k * quantum2 ;
      minlh = (xy[j][i] < minlh) ? xy[j][i] : minlh ;
      maxlh = (xy[j][i] > maxlh) ? xy[j][i] : maxlh ;
    }
    for (i=NPTS/2;i<NPTS;i++) {       // HH quadrant, level 1
      k = xy[j][i] / quantum4 + .5f ;
      if(use_quantum) xy[j][i] = k * quantum4 ;
      minhh = (xy[j][i] < minhh) ? xy[j][i] : minhh ;
      maxhh = (xy[j][i] > maxhh) ? xy[j][i] : maxhh ;
    }
  }

  if(use_quantum){
    printf("Transformed 2D signal (after quantification by %8.2f):\n",quantum) ;
  }else{
    printf("Transformed 2D signal\n");
  }
  for (j=NPTS-1;j>=0;j--) {
    for (i=0;i<NPTS;i++) printf(" %8.2f",xy[j][i]);
    printf("\n");
  }
  quadrants(&xy[0][0], NPTS, NPTS, mins, maxs) ;
//   printf("max(1) = %8.2f, %8.2f, %8.2f, %8.2f\n", maxs[0], maxs[1], maxs[2], maxs[3]) ;
//   printf("min(1) = %8.2f, %8.2f, %8.2f, %8.2f\n", mins[0], mins[1], mins[2], mins[3]) ;
  printf("level 1 (%8.2f, %8.2f) (%8.2f, %8.2f)\n", mins[2], maxs[2], mins[3], maxs[3]) ;
  if(levels == 1) printf("        (%8.2f, %8.2f) (%8.2f, %8.2f)\n", mins[0], maxs[0], mins[1], maxs[1]) ;
  if(levels >  1) printf("                             (%8.2f, %8.2f)\n", mins[1], maxs[1]) ;
  if(levels > 1){
    quadrants(&xy[0][0], NPTS/2, NPTS, mins, maxs) ;
//     printf("max(2) = %8.2f, %8.2f, %8.2f, %8.2f\n", maxs[0], maxs[1], maxs[2], maxs[3]) ;
//     printf("min(2) = %8.2f, %8.2f, %8.2f, %8.2f\n", mins[0], mins[1], mins[2], mins[3]) ;
    printf("level 2 (%8.2f, %8.2f) (%8.2f, %8.2f)\n", mins[2], maxs[2], mins[3], maxs[3]) ;
    if(levels == 2) printf("        (%8.2f, %8.2f) (%8.2f, %8.2f)\n", mins[0], maxs[0], mins[1], maxs[1]) ;
    if(levels >  2) printf("                             (%8.2f, %8.2f)\n", mins[1], maxs[1]) ;
  }
  if(levels > 2){
    quadrants(&xy[0][0], NPTS/2, NPTS, mins, maxs) ;
//     printf("max(3) = %8.2f, %8.2f, %8.2f, %8.2f\n", maxs[0], maxs[1], maxs[2], maxs[3]) ;
//     printf("min(3) = %8.2f, %8.2f, %8.2f, %8.2f\n", mins[0], mins[1], mins[2], mins[3]) ;
    printf("level 3 (%8.2f, %8.2f) (%8.2f, %8.2f)\n", mins[2], maxs[2], mins[3], maxs[3]) ;
    printf("        (%8.2f, %8.2f) (%8.2f, %8.2f)\n", mins[0], maxs[0], mins[1], maxs[1]) ;
  }
//   printf("max = %8.2f, %8.2f, %8.2f, %8.2f\n", maxll, maxhl, maxlh, maxhh) ;
//   printf("min = %8.2f, %8.2f, %8.2f, %8.2f\n", minll, minhl, minlh, minhh) ;
  if(use_quantum) printf("quanta = %8.3f, %8.3f, %8.3f, %8.3f\n", quantum, quantum2, quantum2, quantum4) ;
  spanll = (maxll-minll)/quantum ;
  spanhl = (maxhl-minhl)/quantum2 ;
  spanlh = (maxlh-minlh)/quantum2 ;
  spanhh = (maxhh-minhh)/quantum4 ;
  if(use_quantum) printf("spanll = %8d, spanhl = %8d, spanlh = %8d, spanhh = %8d\n", spanll, spanhl, spanlh, spanhh) ;

#if defined(CDF97)
  I_CDF97_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, levels);
//   I_CDF97_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
//   I_CDF97_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   I_CDF97_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
//   printf("Restored 2D signal:\n");
//   for (j=NPTS-1;j>=0;j--) {
//     for (i=0;i<NPTS;i++) printf(" %8.3f",xy[j][i]);
//     printf("\n");
//   }
#else
  I_CDF53_2D_split_inplace_n((float *)xy, NPTS,  NPTS, NPTS, levels);
//   I_CDF53_2D_split_inplace((float *)xy, npts4, NPTS, npts4);
//   I_CDF53_2D_split_inplace((float *)xy, npts2, NPTS, npts2);
//   I_CDF53_2D_split_inplace((float *)xy, NPTS,  NPTS, NPTS);
#endif
  if(use_quantum) {
    printf("Restored 2D signal error in units of quantum (%5.2f):\n", quantum);
  }else{
    printf("Restored 2D signal error :\n");
    quantum = 1.0f ;
  }
  for (j=NPTS-1;j>=0;j--) {
//     for (i=0;i<NPTS;i++) printf(" %8.2f",(xy[j][i] - (3+i+0.4*i*i-0.02*i*i*i) * (3+i+0.4*j*j-0.02*j*j*j)) / quantum);
    for (i=0;i<NPTS;i++) printf(" %8.1e",(xy[j][i] - xy0[j][i]) / quantum);
    printf("\n");
  }

}
