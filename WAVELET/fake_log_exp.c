#include <stdio.h>

typedef union{
  float f;
  int i;
  } intfloat;

float fake_log(intfloat x) {  /* |x.f| must be  >= 1.0 (inverse function of fake_exp) */
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

float fake_exp(float f){   /* |result| will be >= 1.0 (inverse function of fake_log) */
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

  x.i = 127 << 23;
  x.i += 1<<22;
  printf("%f %f %f\n",x.f,fake_exp(zero),fake_log((intfloat)zero));
  for( z = 1.0000001 ; z < 5000.0 ; z = z*-1.17) { 
    x.f = (z/fake_exp(fake_log((intfloat)z))) ;
    printf("%12.6f %12.6f %13.8f %8d\n", z, fake_log((intfloat)z), x.f, (x.i << 9) >> 9 ); 
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
