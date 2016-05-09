#include <stdio.h>

typedef union{
  float f;
  int i;
  } intfloat;

void v_fake_log(float *X, float min0, int n) {  /* |x.f| will be adjusted  >= 1.0 (inverse function of fake_exp) */
  int exp;
  int sign;
  intfloat x;
  int j;
  int lowexp, offsetexp;;

  x.f = min0;
  lowexp = 0xFF & (x.i >> 23);                /* exponent of smallest non zero absolute value */
  offsetexp = 127 - lowexp;                   /* will be added to exponents to make sure that "X" >= 1.0 */

  for(j = 0 ; j < n ; j++) {
    x.f = X[j];
    sign = x.i & (1 << 31);                   /* extract sign */
    exp = 0xFF & (x.i >> 23);                 /* IEEE exponent */
    exp = exp - 128 + offsetexp;              /* get log2(x) -1  from IEEE exponent */
    x.i = (x.i & 0x3FFFFFFF) | 0x3f800000;    /* 1.0 <= X < 2.0 (force exponent to 127) */
    x.f = x.f + exp;                          /* pseudo log  */
    x.i = x.i | sign;                         /* restore sign */
    X[j] = x.f;
  }
}

static float fake_log(intfloat x, float min0) {  /* |x.f| must be  >= 1.0 (inverse function of fake_exp) */
  int exp;
  int sign;

  sign = x.i & (1 << 31);                   /* extract sign */
  exp = 0xFF & (x.i >> 23);                 /* IEEE exponent */
  exp = exp - 128;                          /* get log2(x) -1  from IEEE exponent */
  x.i = (x.i & 0x3FFFFFFF) | 0x3f800000;    /* 1.0 <= X < 2.0 (force exponent to 127) */
  x.f = x.f + exp;                          /* pseudo log  */
  x.i = x.i | sign;                         /* restore sign */
  return (x.f);
}

void v_fake_exp(float *f, float min0, int n){   /* |result| will be >= 1.0 (inverse function of fake_log) */
  intfloat x;
  int exp;
  int sign;
  int j;
  int lowexp;

  x.f = min0;
  lowexp = 0xFF & (x.i >> 23);                      /* exponent of smallest non zero absolute value */

  for(j = 0 ; j < n ; j++) {
    x.f = f[j];
    sign = x.i & (1 << 31);                         /* save sign */
    x.i = x.i & 0x7FFFFFFF;                         /* strip sign */
    exp = x.f;                                      /* get log2-1  */
    x.f = x.f - exp + 1;                            /*  1.0 <= X < 2.0 */
    exp = exp + lowexp;                             /* IEEE exponent with proper offset */
    x.i = (x.i & 0x7FFFFF ) | (exp << 23) | sign;   /* rebuild float and restore sign */
    f[j] = x.f;
  }
}

static float fake_exp(float f, float min0){   /* |result| will be >= 1.0 (inverse function of fake_log) */
  intfloat x;
  int exp;
  int sign;

  x.f = f;
  sign = x.i & (1 << 31);                         /* save sign */
  x.i = x.i & 0x7FFFFFFF;                         /* strip sign */
  exp = x.f;                                      /* get log2-1  */
  x.f = x.f - exp + 1;                            /*  1.0 <= X < 2.0 */
  exp += 127;                                     /* IEEE exponent with proper offset */
  x.i = (x.i & 0x7FFFFF ) | (exp << 23) | sign;   /* rebuild float and restore sign */
  return (x.f);
}

#if defined(SELF_TEST)
main() {
  union{
  float f;
  int i;
  } x;
  float z, zz;
  float zero = 0.0;
  float z0 = 2.0000001;

  x.i = 127 << 23;
  x.i += 1<<22;
  printf("%f %f %f\n",x.f,fake_exp(zero,z0),fake_log((intfloat)zero,z0));
  for( z = z0 ; z < 5000.0 ; z = z*-1.17) { 

    x.f = (z/fake_exp(fake_log((intfloat)z,z0),z0)) ;
    printf("S: %12.6f %12.6f %13.8f %8d\n", z, fake_log((intfloat)z,z0), x.f, (x.i << 9) >> 9 ); 

    zz = z ; v_fake_log(&zz,z0,1) ; x.f = zz ; v_fake_exp(&x.f,z0,1);
    x.f = z / x.f;
    printf("V: %12.6f %12.6f %13.8f %8d\n", z, zz, x.f, (x.i << 9) >> 9 ); 
  }

  x.i = 0x3f800005;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
  x.i -= 1;
  printf("%13.8f  %8.8x\n",x.f,x.i);
}
#endif
