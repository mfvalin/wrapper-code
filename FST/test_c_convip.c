#include <stdio.h>
#include <stdlib.h>
#include <convert_ip.h>
#include <convert_ip.h>


int main()
{
  float p;
  int ip =850;
  int kind=-1;
  int mode;

  ConvIp(&ip,&p,&kind,-1);
  fprintf(stderr," p=%g, kind=%d, ip=%d \n",p,kind,ip);
  p=850.0;
  kind=0;
  mode = 2;
  ConvIp(&ip,&p,&kind,mode);
  fprintf(stderr," p=%g, kind=%d, ip=%d \n",p,kind,ip);
  p=-1.0;
  kind=-1;
  ConvIp(&ip,&p,&kind,-1);
  fprintf(stderr," p=%g, kind=%d, ip=%d \n",p,kind,ip);
  return 0;
}
