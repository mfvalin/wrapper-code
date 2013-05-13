#include <stdio.h>
#include <stdlib.h>
#include <convert_ip.h>
#include <convert_ip.h>


int Cmain()  /* called by Fortran test */
{
  float p;
  int ip;
  int kind;
  int mode;
  int status;
  int ip1,ip2,ip3,kind1,kind2,kind3;
  int ip123[3], kind123[3];
  float p1,p2,p3;
  float p123[3];
  float_ip fp1, fp2, fp3;
  float_ip fp123[3];

  fprintf(stderr,"old style ip -> p,kind single value conversion\n");
  ip=850;
  p=0.0;
  kind=-1;
  ConvertIp(&ip,&p,&kind,-1);
  fprintf(stderr,"ip=%d, p=%g, kind=%d\n\n",ip,p,kind);

  ip=-1;
  fprintf(stderr,"p,kind -> ip single value conversion\n");
  ConvertIp(&ip,&p,&kind,2);
  fprintf(stderr,"p=%g, kind=%d, ip=%d\n\n",p,kind,ip);

  fprintf(stderr,"p,kind -> ip single value conversion\n");
  p=850.0;
  kind=0;
  ip=-1;
  mode = 2;
  ConvertIp(&ip,&p,&kind,mode);
  fprintf(stderr,"p=%g, kind=%d, ip=%d\n\n",p,kind,ip);

  fprintf(stderr,"ip -> p,kind single value conversion\n");
  p=-1.0;
  kind=-1;
  ConvertIp(&ip,&p,&kind,-1);
  fprintf(stderr,"ip=%d, p=%g, kind=%d\n\n",ip,p,kind);
  
  fprintf(stderr,"==============================================================\n\n");
  
  fprintf(stderr,"multi value conversion, 1 level, 2 times in the wrong order\n");
  
  p1=950.0 ; kind1 = 2 ; fp1.v1=950.0 ; fp1.v2=950.0 ; fp1.kind=2 ;
  p2=48.0  ; kind2 = 10; fp2.v1= 48.0 ; fp2.v2= 24.0 ; fp2.kind=10;
  p3=24.0  ; kind3 = 10; fp3.v1=  0.0 ; fp3.v2=  0.0 ; fp3.kind=-1;
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d\n",p1,kind1,p2,kind2,p3,kind3);
  status=ConvertPKtoIP(&ip1,&ip2,&ip3,p1,p2,p3,kind1,kind2,kind3);
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);

  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind);
  status=EncodeIp(&ip1,&ip2,&ip3,&fp1,&fp2,&fp3);
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);
  
  p1=0.0 ; kind1 = -1; fp1=(float_ip){0.0,0.0,-1};
  p2=0.0 ; kind2 = -1; fp2=(float_ip){0.0,0.0,-1};
  p3=0.0 ; kind3 = -1; fp3=(float_ip){0.0,0.0,-1};
  status=ConvertIPtoPK(&p1,&p2,&p3,&kind1,&kind2,&kind3,ip1,ip2,ip3);
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);
  status=DecodeIp(&fp1,&fp2,&fp3,ip1,ip2,ip3);
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  fprintf(stderr,"==============================================================\n\n");
  
  fprintf(stderr,"multi value conversion, 2 pressure levels in the wrong order, 1 time\n");
  
  p1=850.0 ; kind1 = 2 ; fp1.v1=850.0 ; fp1.v2=950.0 ; fp1.kind=2 ;
  p2=950.0 ; kind2 = 2 ; fp2.v1= 24.0 ; fp2.v2= 24.0 ; fp2.kind=10 ;
  p3=24.0  ; kind3 = 10; fp3.v1=  0.0 ; fp3.v2=  0.0 ; fp3.kind=-1;
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d\n",p1,kind1,p2,kind2,p3,kind3);
  status=ConvertPKtoIP(&ip1,&ip2,&ip3,p1,p2,p3,kind1,kind2,kind3);
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);

  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind);
  status=EncodeIp(&ip1,&ip2,&ip3,&fp1,&fp2,&fp3);
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);
  
  p1=0.0 ; kind1 = -1; fp1=(float_ip){0.0,0.0,-1};
  p2=0.0 ; kind2 = -1; fp2=(float_ip){0.0,0.0,-1};
  p3=0.0 ; kind3 = -1; fp3=(float_ip){0.0,0.0,-1};
  status=ConvertIPtoPK(&p1,&p2,&p3,&kind1,&kind2,&kind3,ip1,ip2,ip3);
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);
  status=DecodeIp(&fp1,&fp2,&fp3,ip1,ip2,ip3);
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  fprintf(stderr,"==============================================================\n\n");
  
  fprintf(stderr,"multi value conversion, 2 heights in the wrong order, 1 time\n");
  
  p1=950.0 ; kind1 = 0 ; fp1=(float_ip){950.0,850.0,0};
  p2=850.0 ; kind2 = 0 ; fp2=(float_ip){24.0,24.0,10};
  p3=24.0  ; kind3 = 10; fp3=(float_ip){0.0,0.0,-1};
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d\n",p1,kind1,p2,kind2,p3,kind3);
  status=ConvertPKtoIP(&ip1,&ip2,&ip3,p1,p2,p3,kind1,kind2,kind3);
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);

  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind);
  status=EncodeIp(&ip1,&ip2,&ip3,&fp1,&fp2,&fp3);
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);
  
  p1=0.0 ; kind1 = -1; fp1=(float_ip){0.0,0.0,-1};
  p2=0.0 ; kind2 = -1; fp2=(float_ip){0.0,0.0,-1};
  p3=0.0 ; kind3 = -1; fp3=(float_ip){0.0,0.0,-1};
  status=ConvertIPtoPK(&p1,&p2,&p3,&kind1,&kind2,&kind3,ip1,ip2,ip3);
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);
  status=DecodeIp(&fp1,&fp2,&fp3,ip1,ip2,ip3);
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  return 0;
}
