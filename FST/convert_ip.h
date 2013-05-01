#ifndef CONVERT_IP_DEFS
#define CONVERT_IP_DEFS
#define COMVERT_IP_OK 0
#define COMVERT_IP_ERROR -1

typedef struct {
  float v1;
  float v2;
  int kind;
} float_ip;

void ConvIp(int *ip, float *p, int *kind, int mode);

int EncodeIp( int *ip1, int *ip2, int *ip3, float_ip *p1, float_ip *p2, float_ip *p3);
int DecodeIp(int *ip1, int *ip2, int *ip3, float_ip *p1, float_ip *p2, float_ip *p3); /* ip as pointers for better FORTRAN compatibility */

int ConvertPKtoIP(int *ip1, int *ip2, int *ip3, float *p1, float *p2, float *p3, int *kind1, int *kind2, int *kind3);
int ConvertIPtoPKP(float *p1, float *p2, float *p3, int *kind1, int *kind2, int *kind3, int *ip1, int *ip2, int *ip3);

#endif
