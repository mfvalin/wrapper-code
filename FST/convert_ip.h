#ifndef CONVERT_IP_DEFS
#define CONVERT_IP_DEFS

#define TO_IP 1
#define TO_RP -1
#define CONVERT_OK 0
#define CONVERT_GUESS 14
#define CONVERT_GOOD_GUESS 2
#define CONVERT_BAD_GUESS 4
#define CONVERT_TERRIBLE_GUESS 8
#define CONVERT_WARNING 32
#define CONVERT_ERROR 64

typedef struct { /* if v1 == v2, it is not a range but a single value */
  float v1;      /* first value of range */
  float v2;      /* second value of range */
  int kind;
} float_ip;

static float_ip invalid_float_ip={0.0,0.0,-1};
#define NULL_float_ip &invalid_float_ip

void ConvIp(int *ip, float *p, int *kind, int mode);

int EncodeIp( int *ip1, int *ip2, int *ip3, float_ip *p1, float_ip *p2, float_ip *p3);
int DecodeIp(float_ip *p1, float_ip *p2, float_ip *p3, int ip1, int ip2, int ip3);

int EncodeIp_v(int ip[3],float_ip p[3]);
int DecodeIp_v(float_ip p[3],int ip[3]);

int ConvertPKtoIP(int *ip1, int *ip2, int *ip3, float p1, float p2, float p3, int kind1, int kind2, int kind3);
int ConvertIPtoPK(float *p1, float *p2, float *p3, int *kind1, int *kind2, int *kind3, int ip1, int ip2, int ip3);

int ConvertPKtoIP_v(int ip[3],float p[3],int kind[3]);
int ConvertIPtoPK_v(float p[3],int kind[3],int ip[3]);

#endif
