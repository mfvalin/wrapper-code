#include <stdio.h>
#include <stdlib.h>
#include <convert_ip.h>
#include <convert_ip.h>

#define EPSILON 0.00001
#define TRUE    1
#define FALSE   0


int Cmain1()  /* called by Fortran test */
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
  ip_info fp1, fp2, fp3;
  ip_info fp123[3];

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
  p2=24.0  ; kind2 = 10; fp2.v1= 48.0 ; fp2.v2= 24.0 ; fp2.kind=10;
  p3=48.0  ; kind3 = 10; fp3.v1=  0.0 ; fp3.v2=  0.0 ; fp3.kind=-1;

  fprintf(stderr,"FROM     pk1=%f %d, pk2=%f,%d, pk3=%f %d\n",p1,kind1,p2,kind2,p3,kind3);
  status=ConvertPKtoIP(&ip1,&ip2,&ip3,p1,kind1,p2,kind2,p3,kind3);
  if(ip2>ip3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);

  fprintf(stderr,"FROM     fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind);
  status=EncodeIp(&ip1,&ip2,&ip3,&fp1,&fp2,NULL_ip_info);  /* fp3 is not used , use NULL_ip_info */
  if(ip2>ip3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);
  
  p1=0.0 ; kind1 = -1; fp1=(ip_info){0.0,0.0,-1};
  p2=0.0 ; kind2 = -1; fp2=(ip_info){0.0,0.0,-1};
  p3=0.0 ; kind3 = -1; fp3=(ip_info){0.0,0.0,-1};

  fprintf(stderr,"FROM     ip1=%d, ip2=%d, ip3=%d\n",ip1,ip3,ip2);
  status=ConvertIPtoPK(&p1,&kind1,&p2,&kind2,&p3,&kind3,ip1,ip3,ip2);  /* ip2 and ip3 deliberately reversed */
  if(p2>p3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);

  status=DecodeIp(&fp1,&fp2,&fp3,ip1,ip3,ip2);  /* ip2 and ip3 deliberately reversed */
  if(fp2.v2>fp2.v1) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  fprintf(stderr,"FROM     ip1=%d, ip2=%d, ip3=%d\n",ip1,ip2,ip3);
  status=ConvertIPtoPK(&p1,&kind1,&p2,&kind2,&p3,&kind3,ip1,ip2,ip3);
  if(p2>p3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);

  status=DecodeIp(&fp1,&fp2,&fp3,ip1,ip2,ip3);
  if(fp2.v2>fp2.v1) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  fprintf(stderr,"==============================================================\n\n");
  
  fprintf(stderr,"multi value conversion, 2 pressure levels in the wrong order, 1 time\n");
  
  p1=850.0 ; kind1 = 2 ; fp1.v1=850.0 ; fp1.v2=950.0 ; fp1.kind=2 ;
  p2=950.0 ; kind2 = 2 ; fp2.v1= 24.0 ; fp2.v2= 24.0 ; fp2.kind=10 ;
  p3=24.0  ; kind3 = 10; INIT_ip_info(fp3);

  status=ConvertPKtoIP(&ip1,&ip2,&ip3,p1,kind1,p2,kind2,p3,kind3);
  fprintf(stderr,"FROM     pk1=%f %d, pk2=%f,%d, pk3=%f %d\n",p1,kind1,p2,kind2,p3,kind3);
  if(ip1>ip3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);

  fprintf(stderr,"FROM     fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind);

  status=EncodeIp(&ip1,&ip2,&ip3,&fp1,&fp2,NULL_ip_info);  /* fp3 is not used , use NULL_ip_info */
  if(ip1>ip3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);
  
  p1=0.0 ; kind1 = -1; fp1=(ip_info){0.0,0.0,-1};
  p2=0.0 ; kind2 = -1; fp2=(ip_info){0.0,0.0,-1};
  p3=0.0 ; kind3 = -1; fp3=(ip_info){0.0,0.0,-1};

  fprintf(stderr,"FROM     ip1=%d, ip2=%d, ip3=%d\n",ip3,ip2,ip1);
  status=ConvertIPtoPK(&p1,&kind1,&p2,&kind2,&p3,&kind3,ip3,ip2,ip1);  /* ip1 and ip3 deliberately reversed */
  if(p1>p3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);

  status=DecodeIp(&fp1,&fp2,&fp3,ip3,ip2,ip1);  /* ip1 and ip3 deliberately reversed */
  if(fp1.v1>fp1.v2) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  fprintf(stderr,"FROM     ip1=%d, ip2=%d, ip3=%d\n",ip1,ip2,ip3);
  status=ConvertIPtoPK(&p1,&kind1,&p2,&kind2,&p3,&kind3,ip1,ip2,ip3);
  if(p1>p3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);

  status=DecodeIp(&fp1,&fp2,&fp3,ip1,ip2,ip3);
  if(fp1.v1>fp1.v2) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  fprintf(stderr,"==============================================================\n\n");
  
  fprintf(stderr,"multi value conversion, 2 heights in the wrong order, 1 time\n");
  
  p1=950.0 ; kind1 = 0 ; fp1=(ip_info){950.0,850.0,0};
  p2=850.0 ; kind2 = 0 ; fp2=(ip_info){24.0,24.0,10};
  p3=24.0  ; kind3 = 10; fp3=(ip_info){0.0,0.0,-1};
  fprintf(stderr,"FROM     pk1=%f %d, pk2=%f,%d, pk3=%f %d\n",p1,kind1,p2,kind2,p3,kind3);
  status=ConvertPKtoIP(&ip1,&ip2,&ip3,p1,kind1,p2,kind2,p3,kind3);
  if(ip1<ip3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);

  fprintf(stderr,"FROM     fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind);
  status=EncodeIp(&ip1,&ip2,&ip3,&fp1,&fp2,NULL_ip_info);  /* fp3 is not used , put NULL_ip_info */
  if(ip1<ip3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," ip1=%d, ip2=%d, ip3=%d, status=%d\n\n",ip1,ip2,ip3,status);
  
  p1=0.0 ; kind1 = -1; fp1=(ip_info){0.0,0.0,-1};
  p2=0.0 ; kind2 = -1; fp2=(ip_info){0.0,0.0,-1};
  p3=0.0 ; kind3 = -1; fp3=(ip_info){0.0,0.0,-1};
  fprintf(stderr,"FROM     ip1=%d, ip2=%d, ip3=%d, status=%d\n",ip3,ip2,ip1,status);
  status=ConvertIPtoPK(&p1,&kind1,&p2,&kind2,&p3,&kind3,ip3,ip2,ip1);  /* ip1 and ip3 deliberately reversed */
  if(p1<p3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);
  
  status=DecodeIp(&fp1,&fp2,&fp3,ip3,ip2,ip1);  /* ip1 and ip3 deliberately reversed */
  if(fp1.v1<fp1.v2) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  fprintf(stderr,"FROM     ip1=%d, ip2=%d, ip3=%d, status=%d\n",ip1,ip2,ip3,status);
  status=ConvertIPtoPK(&p1,&kind1,&p2,&kind2,&p3,&kind3,ip1,ip2,ip3);
  if(p1<p3) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," pk1=%f %d, pk2=%f,%d, pk3=%f %d, status=%d\n",p1,kind1,p2,kind2,p3,kind3,status);
  
  status=DecodeIp(&fp1,&fp2,&fp3,ip1,ip2,ip3);
  if(fp1.v1<fp1.v2) fprintf(stderr,"SUCCESS "); else fprintf(stderr,"FAILURE ");
  fprintf(stderr," fp1=%f %f %d, fp2=%f %f %d, pk3=%f %f %d, status=%d\n\n",
	    fp1.v1,fp1.v2,fp1.kind,fp2.v1,fp2.v2,fp2.kind,fp3.v1,fp3.v2,fp3.kind,status);
  
  return 0;
}

void Test1(void)
{
   float p    =   0.0;
   int   ip   = 850;
   int   kind =  -1;
 
   fprintf(stderr, "Test #1 : Decode un IP de %d.\n", ip);

   ConvertIp(&ip, &p, &kind, -1);

   if ((ip != 850) || ((p - 850.0) > EPSILON) || (kind != 2))
      fprintf(stderr, "\033[31;1mErreur dans le test #1\033[0m\n");
 
   fprintf(stderr,"ip= %d, p= %f, kind= %d\n\n", ip, p, kind);
}


void Test2(void)
{
   float p    = 850.0;
   int   ip   = -1;
   int   kind = 2;
 
   fprintf(stderr, "Test #2 : Encode %f avec kind=%d dans le NEWSTYLE\n", p, kind);

   ConvertIp(&ip, &p, &kind, 2);

   if ((ip != 41744464) || ((p - 850.0) > EPSILON) || (kind != 2))
      fprintf(stderr, "\033[31;1mErreur dans le test #2\033[0m\n");
 
   fprintf(stderr,"ip= %d, p= %f, kind= %d\n\n", ip, p, kind);
}


void Test3(void)
{
   float p    = 0.3;
   int   ip   = -1;
   int   kind = 1;
 
   fprintf(stderr, "Test #3 : Encode %f avec kind=%d dans le NEWSTYLE\n", p, kind);

   ConvertIp(&ip, &p, &kind, 2);

   if ((ip != 27562976) || ((p - 0.3) > EPSILON) || (kind != 1))
      fprintf(stderr, "\033[31;1mErreur dans le test #3\033[0m\n");
 
   fprintf(stderr,"ip= %d, p= %f, kind= %d\n\n", ip, p, kind);
}


void Test4(void)
{
   float p    = 0.0;
   int   ip   = 107189252;
   int   kind = -1;
 
   fprintf(stderr, "Test #4 : Decode un IP de %d \n", ip);

   ConvertIp(&ip, &p, &kind, -1);

   if ((ip != 107189252) || ((p - 2345.0) > EPSILON) || (kind != 6))
      fprintf(stderr, "\033[31;1mErreur dans le test #4\033[0m\n");
 
   fprintf(stderr,"ip= %d, p= %f, kind= %d\n\n", ip, p, kind);
}


void Test100(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){2345.0, 2345.0,  6};
   ip_info fp2 = (ip_info){  12.0,   24.0, 10};
   int      ip1toCmp = 107189252;  /* 2345.0 */
   int      ip2toCmp = 176400768;  /*   24.0 */
   int      ip3toCmp = 176280768;  /*   12.0 */
 
   fprintf(stderr, "Test #100 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, NULL_ip_info);  

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #100\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test101(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){2345.0, 2345.0,  6};
   ip_info fp2 = (ip_info){  24.0,   12.0, 10};
   int      ip1toCmp = 107189252;  /* 2345.0 */
   int      ip2toCmp = 176400768;  /*   24.0 */
   int      ip3toCmp = 176280768;  /*   12.0 */
 
   fprintf(stderr, "Test #101 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, NULL_ip_info);  

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #101\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test102(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){48.0, 12.0, 10};
   ip_info fp2 = (ip_info){ 0.3,  0.3,  1};
   int      ip1toCmp = 27562976;   /* 0.3 */
   int      ip2toCmp = 176640768;  /* 48.0 */
   int      ip3toCmp = 176280768;  /* 12.0 */
 
   fprintf(stderr, "Test #102 : Erreur en essayant d'encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, NULL_ip_info);  

   if (status == CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #102\033[0m\n");
      fprintf(stderr, "\033[31;1mif (status != CONVERT_OK)\033[0m\n");
   }
 
   fprintf(stderr, "status=%d\n\n", status);
}


void Test103(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){ 0.3,  0.7,  1};
   ip_info fp2 = (ip_info){12.0, 12.0, 10};
   int      ip1toCmp = 27962976;  /* 0.7  */
   int      ip2toCmp = 176280768; /* 12.0 */
   int      ip3toCmp = 27562976;  /* 0.3  */
 
   fprintf(stderr, "Test #103 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, NULL_ip_info);  

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #103\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test104(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){ 0.7,  0.3,  1};
   ip_info fp2 = (ip_info){12.0, 12.0, 10};
   int      ip1toCmp = 27962976;  /* 0.7  */
   int      ip2toCmp = 176280768; /* 12.0 */
   int      ip3toCmp = 27562976;  /* 0.3  */
 
   fprintf(stderr, "Test #104 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, NULL_ip_info);  

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #104\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test105(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){100.0, 560.0,  4};
   ip_info fp2 = (ip_info){ 12.0,  12.0, 10};
   int      ip1toCmp = 74548896;  /* 100.0  */
   int      ip2toCmp = 176280768; /*  12.0  */
   int      ip3toCmp = 75008896;  /* 560.0  */
 
   fprintf(stderr, "Test #105 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, NULL_ip_info);  

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #105\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test106(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){100.0, 560.0,  2};
   ip_info fp2 = (ip_info){ 12.0,  12.0, 10};
   int      ip1toCmp = 41454464;  /* 560.0  */
   int      ip2toCmp = 176280768; /*  12.0  */
   int      ip3toCmp = 40994464;  /* 100.0  */
 
   fprintf(stderr, "Test #106 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, NULL_ip_info);  

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #106\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test107(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){100.0, 560.0,  4};
   ip_info fp2 = (ip_info){ 12.0,  24.0, 10};
 
   fprintf(stderr, "Test #107 : Erreur en essayant d'encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, NULL_ip_info);  

   if (status == CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #107\033[0m\n");
      fprintf(stderr, "\033[31;1mif (status != CONVERT_OK)\033[0m\n");
   }
   fprintf(stderr,"status=%d\n\n", status);
}


void Test108(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){100.0, 100.0,  2};
   ip_info fp2 = (ip_info){ 12.0,  12.0, 10};
   ip_info fp3 = (ip_info){ 36.0,  36.0,  3};
   int      ip1toCmp = 40994464;   /* 100.0  */
   int      ip2toCmp = 176280768;  /*  12.0  */
   int      ip3toCmp = 59080256;   /*  36.0  */
 
   fprintf(stderr, "Test #108 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d}, {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, &fp3);  

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #108\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test109(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){ 188.0,  188.0,  2};
   ip_info fp2 = (ip_info){  12.0,   12.0, 10};
   ip_info fp3 = (ip_info){1013.0, 1013.0,  2};
   int      ip1toCmp = 41082464;   /*  188.0  */
   int      ip2toCmp = 176280768;  /*   12.0  */
   int      ip3toCmp = 39947188;   /* 1013.0  */
 
   fprintf(stderr, "Test #109 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d}, {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, &fp3);  

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #109\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test110(void)
{
   int      status;
   int      ip1 = -1; 
   int      ip2 = -1;
   int      ip3 = -1;
   ip_info fp1 = (ip_info){ 188.0,  188.0,  2};
   ip_info fp2 = (ip_info){1013.0, 1013.0,  2};
   ip_info fp3 = (ip_info){  12.0,   12.0, 10};
   int      ip1toCmp = 41082464;   /*  188.0  */
   int      ip2toCmp = 176280768;  /*   12.0  */
   int      ip3toCmp = 39947188;   /* 1013.0  */
 
   fprintf(stderr, "Test #110 : Encode IP1, IP2 et IP3 d'après les structures ip_info {%f, %f, %d}, {%f, %f, %d} et {%f, %f, %d}\n",
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind);

   status = EncodeIp(&ip1, &ip2, &ip3, &fp1, &fp2, &fp3);  

   if (status == CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #110\033[0m\n");
      fprintf(stderr, "\033[31;1mif (status != CONVERT_OK)\033[0m\n");
   }
   fprintf(stderr,"status=%d\n\n", status);
}


int ip_info_equal(ip_info fp1, ip_info fp2)
{
   if (fp1.v1 == fp2.v1 && fp1.v2 == fp2.v2 && fp1.kind == fp2.kind)
      return TRUE;
   else
      return FALSE;
}

void Test200(void)
{
   int   status;
   int   ip1 = 40994464;   /* 100.0  2 */
   int   ip2 = 176280768;  /*  12.0 10 */
   int   ip3 = 59080256;   /*  36.0  3 */
   ip_info fp1 = (ip_info){0.0,0.0,-1};
   ip_info fp2 = (ip_info){0.0,0.0,-1};
   ip_info fp3 = (ip_info){0.0,0.0,-1};
   ip_info fp1toCmp = (ip_info){100.00, 100.0, 2};
   ip_info fp2toCmp = (ip_info){ 12.00,  12.0, 10};
   ip_info fp3toCmp = (ip_info){ 36.00,  36.0, 3};
 
   fprintf(stderr, "Test #200 : DecodeIP d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = DecodeIp(&fp1, &fp2, &fp3, ip1, ip2, ip3);  

   if (ip_info_equal(fp1, fp1toCmp) != TRUE ||
       ip_info_equal(fp2, fp2toCmp) != TRUE ||
       ip_info_equal(fp3, fp3toCmp) != TRUE ||
       status                       != CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #200\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || status != %d)\033[0m\n", fp1.v1, fp1.v2, fp1.kind, fp1toCmp.v1, fp1toCmp.v2, fp1toCmp.kind, fp2.v1, fp2.v2, fp2.kind, fp2toCmp.v1, fp2toCmp.v2, fp2toCmp.kind, fp3.v1, fp3.v2, fp3.kind, fp3toCmp.v1, fp3toCmp.v2, fp3toCmp.kind, CONVERT_OK);
   }

   fprintf(stderr,"fp1={%f, %f, %d}, fp2={%f, %f, %d}, fp3={%f, %f, %d}, status=%d\n\n", 
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind, status);
}


void Test201(void)
{
   int   status;
   int   ip1 = 40994464;   /* 100.0  2 */
   int   ip2 = 176280768;  /*  12.0 10 */
   int   ip3 = 176520768;  /*  36.0 10 */
   ip_info fp1 = (ip_info){0.0,0.0,-1};
   ip_info fp2 = (ip_info){0.0,0.0,-1};
   ip_info fp3 = (ip_info){0.0,0.0,-1};
   ip_info fp1toCmp = (ip_info){100.0, 100.0,  2};
   ip_info fp2toCmp = (ip_info){ 12.0,  36.0, 10};
   ip_info fp3toCmp = (ip_info){  0.0,   0.0, -1};
 
   fprintf(stderr, "Test #201 : DecodeIP d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = DecodeIp(&fp1, &fp2, &fp3, ip1, ip2, ip3);  

   if (ip_info_equal(fp1, fp1toCmp) != TRUE ||
       ip_info_equal(fp2, fp2toCmp) != TRUE ||
       ip_info_equal(fp3, fp3toCmp) != TRUE ||
       status                       != CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #201\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || status != %d)\033[0m\n", fp1.v1, fp1.v2, fp1.kind, fp1toCmp.v1, fp1toCmp.v2, fp1toCmp.kind, fp2.v1, fp2.v2, fp2.kind, fp2toCmp.v1, fp2toCmp.v2, fp2toCmp.kind, fp3.v1, fp3.v2, fp3.kind, fp3toCmp.v1, fp3toCmp.v2, fp3toCmp.kind, CONVERT_OK);
   }

   fprintf(stderr,"fp1={%f, %f, %d}, fp2={%f, %f, %d}, fp3={%f, %f, %d}, status=%d\n\n", 
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind, status);
}


void Test202(void)
{
   int   status;
   int   ip1 = 40994464;   /* 100.0  2 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 176280768;  /*  12.0 10 */
   ip_info fp1 = (ip_info){0.0,0.0,-1};
   ip_info fp2 = (ip_info){0.0,0.0,-1};
   ip_info fp3 = (ip_info){0.0,0.0,-1};
   ip_info fp1toCmp = (ip_info){100.0, 100.0,  2};
   ip_info fp2toCmp = (ip_info){ 12.0,  36.0, 10};
   ip_info fp3toCmp = (ip_info){  0.0,   0.0, -1};
 
   fprintf(stderr, "Test #202 : DecodeIP d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = DecodeIp(&fp1, &fp2, &fp3, ip1, ip2, ip3);  

   if (ip_info_equal(fp1, fp1toCmp) != TRUE ||
       ip_info_equal(fp2, fp2toCmp) != TRUE ||
       ip_info_equal(fp3, fp3toCmp) != TRUE ||
       status                       != CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #202\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || status != %d)\033[0m\n", fp1.v1, fp1.v2, fp1.kind, fp1toCmp.v1, fp1toCmp.v2, fp1toCmp.kind, fp2.v1, fp2.v2, fp2.kind, fp2toCmp.v1, fp2toCmp.v2, fp2toCmp.kind, fp3.v1, fp3.v2, fp3.kind, fp3toCmp.v1, fp3toCmp.v2, fp3toCmp.kind, CONVERT_OK);
   }

   fprintf(stderr,"fp1={%f, %f, %d}, fp2={%f, %f, %d}, fp3={%f, %f, %d}, status=%d\n\n", 
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind, status);
}


void Test203(void)
{
   int   status;
   int   ip1 = 176520768;  /*  36.0 10 */
   int   ip2 = 40994464;   /* 100.0  2 */
   int   ip3 = 176280768;  /*  12.0 10 */
   ip_info fp1 = (ip_info){0.0,0.0,-1};
   ip_info fp2 = (ip_info){0.0,0.0,-1};
   ip_info fp3 = (ip_info){0.0,0.0,-1};
   ip_info fp1toCmp = (ip_info){ 36.0,  36.0, 10};
   ip_info fp2toCmp = (ip_info){100.0, 100.0,  2};
   ip_info fp3toCmp = (ip_info){ 12.0,  12.0, 10};
 
   fprintf(stderr, "Test #203 : DecodeIP d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = DecodeIp(&fp1, &fp2, &fp3, ip1, ip2, ip3);  

   if (status == CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #203\033[0m\n");
      fprintf(stderr, "\033[31;1mif (status != CONVERT_OK)\033[0m\n");
   }
   fprintf(stderr,"status=%d\n\n", status);
}


void Test204(void)
{
   int   status;
   int   ip1 = 40994464;   /* 100.0  2 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 41219464;   /* 325.0  2 */
   ip_info fp1 = (ip_info){0.0,0.0,-1};
   ip_info fp2 = (ip_info){0.0,0.0,-1};
   ip_info fp3 = (ip_info){0.0,0.0,-1};
   ip_info fp1toCmp = (ip_info){325.0, 100.0,  2};
   ip_info fp2toCmp = (ip_info){ 36.0,  36.0, 10};
   ip_info fp3toCmp = (ip_info){  0.0,   0.0, -1};
 
   fprintf(stderr, "Test #204 : DecodeIP d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = DecodeIp(&fp1, &fp2, &fp3, ip1, ip2, ip3);  

   if (ip_info_equal(fp1, fp1toCmp) != TRUE ||
       ip_info_equal(fp2, fp2toCmp) != TRUE ||
       ip_info_equal(fp3, fp3toCmp) != TRUE ||
       status                       != CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #204\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || status != %d)\033[0m\n", fp1.v1, fp1.v2, fp1.kind, fp1toCmp.v1, fp1toCmp.v2, fp1toCmp.kind, fp2.v1, fp2.v2, fp2.kind, fp2toCmp.v1, fp2toCmp.v2, fp2toCmp.kind, fp3.v1, fp3.v2, fp3.kind, fp3toCmp.v1, fp3toCmp.v2, fp3toCmp.kind, CONVERT_OK);
   }

   fprintf(stderr,"fp1={%f, %f, %d}, fp2={%f, %f, %d}, fp3={%f, %f, %d}, status=%d\n\n", 
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind, status);
}


void Test205(void)
{
   int   status;
   int   ip1 = 41219464;   /* 325.0  2 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 40994464;   /* 100.0  2 */
   ip_info fp1 = (ip_info){0.0,0.0,-1};
   ip_info fp2 = (ip_info){0.0,0.0,-1};
   ip_info fp3 = (ip_info){0.0,0.0,-1};
   ip_info fp1toCmp = (ip_info){325.0, 100.0,  2};
   ip_info fp2toCmp = (ip_info){ 36.0,  36.0, 10};
   ip_info fp3toCmp = (ip_info){  0.0,   0.0, -1};
 
   fprintf(stderr, "Test #205 : DecodeIP d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = DecodeIp(&fp1, &fp2, &fp3, ip1, ip2, ip3);  

   if (ip_info_equal(fp1, fp1toCmp) != TRUE ||
       ip_info_equal(fp2, fp2toCmp) != TRUE ||
       ip_info_equal(fp3, fp3toCmp) != TRUE ||
       status                       != CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #205\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || status != %d)\033[0m\n", fp1.v1, fp1.v2, fp1.kind, fp1toCmp.v1, fp1toCmp.v2, fp1toCmp.kind, fp2.v1, fp2.v2, fp2.kind, fp2toCmp.v1, fp2toCmp.v2, fp2toCmp.kind, fp3.v1, fp3.v2, fp3.kind, fp3toCmp.v1, fp3toCmp.v2, fp3toCmp.kind, CONVERT_OK);
   }

   fprintf(stderr,"fp1={%f, %f, %d}, fp2={%f, %f, %d}, fp3={%f, %f, %d}, status=%d\n\n", 
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind, status);
}


void Test206(void)
{
   int   status;
   int   ip1 = 7665032;    /* 325.0  0 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 8124032;    /* 784.0  0 */
   ip_info fp1 = (ip_info){0.0,0.0,-1};
   ip_info fp2 = (ip_info){0.0,0.0,-1};
   ip_info fp3 = (ip_info){0.0,0.0,-1};
   ip_info fp1toCmp = (ip_info){325.0, 784.0,  0};
   ip_info fp2toCmp = (ip_info){ 36.0,  36.0, 10};
   ip_info fp3toCmp = (ip_info){  0.0,   0.0, -1};
 
   fprintf(stderr, "Test #206 : DecodeIP d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = DecodeIp(&fp1, &fp2, &fp3, ip1, ip2, ip3);  

   if (ip_info_equal(fp1, fp1toCmp) != TRUE ||
       ip_info_equal(fp2, fp2toCmp) != TRUE ||
       ip_info_equal(fp3, fp3toCmp) != TRUE ||
       status                       != CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #206\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || status != %d)\033[0m\n", fp1.v1, fp1.v2, fp1.kind, fp1toCmp.v1, fp1toCmp.v2, fp1toCmp.kind, fp2.v1, fp2.v2, fp2.kind, fp2toCmp.v1, fp2toCmp.v2, fp2toCmp.kind, fp3.v1, fp3.v2, fp3.kind, fp3toCmp.v1, fp3toCmp.v2, fp3toCmp.kind, CONVERT_OK);
   }

   fprintf(stderr,"fp1={%f, %f, %d}, fp2={%f, %f, %d}, fp3={%f, %f, %d}, status=%d\n\n", 
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind, status);
}


void Test207(void)
{
   int   status;
   int   ip1 = 8124032;    /* 784.0  0 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 7665032;    /* 325.0  0 */
   ip_info fp1 = (ip_info){0.0,0.0,-1};
   ip_info fp2 = (ip_info){0.0,0.0,-1};
   ip_info fp3 = (ip_info){0.0,0.0,-1};
   ip_info fp1toCmp = (ip_info){325.0, 784.0,  0};
   ip_info fp2toCmp = (ip_info){ 36.0,  36.0, 10};
   ip_info fp3toCmp = (ip_info){  0.0,   0.0, -1};
 
   fprintf(stderr, "Test #207 : DecodeIP d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = DecodeIp(&fp1, &fp2, &fp3, ip1, ip2, ip3);  

   if (ip_info_equal(fp1, fp1toCmp) != TRUE ||
       ip_info_equal(fp2, fp2toCmp) != TRUE ||
       ip_info_equal(fp3, fp3toCmp) != TRUE ||
       status                       != CONVERT_OK)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #207\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || (ip_info_equal({%f, %f, %d}, {%f, %f, %d}) != TRUE || status != %d)\033[0m\n", fp1.v1, fp1.v2, fp1.kind, fp1toCmp.v1, fp1toCmp.v2, fp1toCmp.kind, fp2.v1, fp2.v2, fp2.kind, fp2toCmp.v1, fp2toCmp.v2, fp2toCmp.kind, fp3.v1, fp3.v2, fp3.kind, fp3toCmp.v1, fp3toCmp.v2, fp3toCmp.kind, CONVERT_OK);
   }

   fprintf(stderr,"fp1={%f, %f, %d}, fp2={%f, %f, %d}, fp3={%f, %f, %d}, status=%d\n\n", 
           fp1.v1, fp1.v2, fp1.kind, fp2.v1, fp2.v2, fp2.kind, fp3.v1, fp3.v2, fp3.kind, status);
}


void Test300(void)
{
   int      status;
   int      ip[3] = {-1, -1, -1};
   ip_info  fp[3] = {(ip_info){2345.0, 2345.0,  6},
                     (ip_info){  12.0,   24.0, 10},
                     (ip_info){   0.0,    0.0, -1}};
   int      ip1toCmp = 107189252;  /* 2345.0 */
   int      ip2toCmp = 176400768;  /*   24.0 */
   int      ip3toCmp = 176280768;  /*   12.0 */
 
   fprintf(stderr, "Test #300 : Encode IP1, IP2 et IP3 d'après le vecteur de structures ip_info {{%f, %f, %d}, {%f, %f, %d}, {%f, %f, %d}}\n",
           fp[0].v1, fp[0].v2, fp[0].kind, fp[1].v1, fp[1].v2, fp[1].kind, fp[2].v1, fp[2].v2, fp[2].kind);

   status = EncodeIp_v(&ip[0], &fp[0]);  

   if ((ip[0] != ip1toCmp) || (ip[1] != ip2toCmp) || (ip[2] != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #300\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip[0] != %d) || (ip[1] != %d) || (ip[2] != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip[0]= %d, ip[1]= %d, ip[2]= %d, status=%d\n\n", ip[0], ip[1], ip[2], status);
}


void Test400(void)
{
}


void Test500(void)
{
   int   status;
   int   ip1 = -1; 
   int   ip2 = -1;
   int   ip3 = -1;
   float p1 = 2345.0; int kind1 =  6;
   float p2 =   12.0; int kind2 = 10;
   float p3 =   24.0; int kind3 = 10;
   int   ip1toCmp = 107189252;  /* 2345.0  */
   int   ip2toCmp = 176400768;  /*   24.0  */
   int   ip3toCmp = 176280768;  /*   12.0  */
 
   fprintf(stderr, "Test #500 : Encode IP1, IP2 et IP3 d'après p1=%f, kind1=%d, p2=%f, kind2=%d, p3=%f, kind3=%d\n",
           p1, kind1, p2, kind2, p3, kind3);

   status = ConvertPKtoIP(&ip1, &ip2, &ip3, p1, kind1, p2, kind2, p3, kind3);

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #500\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test501(void)
{
   int   status;
   int   ip1 = -1; 
   int   ip2 = -1;
   int   ip3 = -1;
   float p1 = 2345.0; int kind1 =  6;
   float p2 =   24.0; int kind2 = 10;
   float p3 =   12.0; int kind3 = 10;
   int   ip1toCmp = 107189252;  /* 2345.0  */
   int   ip2toCmp = 176400768;  /*   24.0  */
   int   ip3toCmp = 176280768;  /*   12.0  */
 
   fprintf(stderr, "Test #501 : Encode IP1, IP2 et IP3 d'après p1=%f, kind1=%d, p2=%f, kind2=%d, p3=%f, kind3=%d\n",
           p1, kind1, p2, kind2, p3, kind3);

   status = ConvertPKtoIP(&ip1, &ip2, &ip3, p1, kind1, p2, kind2, p3, kind3);

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #501\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test502(void)
{
   int   status;
   int   ip1 = -1; 
   int   ip2 = -1;
   int   ip3 = -1;
   float p1 = 1013.0; int kind1 =  2;
   float p2 =  188.0; int kind2 =  2;
   float p3 =   12.0; int kind3 = 10;
   int   ip1toCmp = 39947188;   /* 1013.0  */
   int   ip2toCmp = 176280768;  /*   12.0  */
   int   ip3toCmp = 41082464;   /*  188.0  */
 
   fprintf(stderr, "Test #502 : Encode IP1, IP2 et IP3 d'après p1=%f, kind1=%d, p2=%f, kind2=%d, p3=%f, kind3=%d\n",
           p1, kind1, p2, kind2, p3, kind3);

   status = ConvertPKtoIP(&ip1, &ip2, &ip3, p1, kind1, p2, kind2, p3, kind3);

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || status != CONVERT_WARNING)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #502\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || status != %d)\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_WARNING);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test503(void)
{
   int   status;
   int   ip1 = -1; 
   int   ip2 = -1;
   int   ip3 = -1;
   float p1 =  188.0; int kind1 =  2;
   float p2 = 1013.0; int kind2 =  2;
   float p3 =   12.0; int kind3 = 10;
   int   ip1toCmp = 39947188;   /* 1013.0  */
   int   ip2toCmp = 176280768;  /*   12.0  */
   int   ip3toCmp = 41082464;   /*  188.0  */
 
   fprintf(stderr, "Test #503 : Encode IP1, IP2 et IP3 d'après p1=%f, kind1=%d, p2=%f, kind2=%d, p3=%f, kind3=%d\n",
           p1, kind1, p2, kind2, p3, kind3);

   status = ConvertPKtoIP(&ip1, &ip2, &ip3, p1, kind1, p2, kind2, p3, kind3);

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || status != CONVERT_WARNING)
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #503\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || status != %d)\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_WARNING);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test504(void)
{
   int   status;
   int   ip1 = -1; 
   int   ip2 = -1;
   int   ip3 = -1;
   float p1 =  188.0; int kind1 =  2;
   float p2 =   12.0; int kind2 = 10;
   float p3 = 1013.0; int kind3 =  2;
   int   ip1toCmp = 41082464;   /*  188.0  */
   int   ip2toCmp = 176280768;  /*   12.0  */
   int   ip3toCmp = 39947188;   /* 1013.0  */
 
   fprintf(stderr, "Test #504 : Encode IP1, IP2 et IP3 d'après p1=%f, kind1=%d, p2=%f, kind2=%d, p3=%f, kind3=%d\n",
           p1, kind1, p2, kind2, p3, kind3);

   status = ConvertPKtoIP(&ip1, &ip2, &ip3, p1, kind1, p2, kind2, p3, kind3);

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #504\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test505(void)
{
   int   status;
   int   ip1 = -1; 
   int   ip2 = -1;
   int   ip3 = -1;
   float p1 = 1013.0; int kind3 =  2;
   float p2 =   12.0; int kind2 = 10;
   float p3 =  188.0; int kind1 =  2;
   int   ip1toCmp = 39947188;   /* 1013.0  */
   int   ip2toCmp = 176280768;  /*   12.0  */
   int   ip3toCmp = 41082464;   /*  188.0  */
 
   fprintf(stderr, "Test #505 : Encode IP1, IP2 et IP3 d'après p1=%f, kind1=%d, p2=%f, kind2=%d, p3=%f, kind3=%d\n",
           p1, kind1, p2, kind2, p3, kind3);

   status = ConvertPKtoIP(&ip1, &ip2, &ip3, p1, kind1, p2, kind2, p3, kind3);

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #505\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test506(void)
{
   int   status;
   int   ip1 = -1; 
   int   ip2 = -1;
   int   ip3 = -1;
   float p1 = 124.0; int kind1 =  4;
   float p2 =  12.0; int kind2 = 10;
   float p3 = 836.0; int kind3 =  4;
   int   ip1toCmp = 74572896;   /* 124.0  */
   int   ip2toCmp = 176280768;  /*  12.0  */
   int   ip3toCmp = 75284896;   /* 836.0  */
 
   fprintf(stderr, "Test #506 : Encode IP1, IP2 et IP3 d'après p1=%f, kind1=%d, p2=%f, kind2=%d, p3=%f, kind3=%d\n",
           p1, kind1, p2, kind2, p3, kind3);

   status = ConvertPKtoIP(&ip1, &ip2, &ip3, p1, kind1, p2, kind2, p3, kind3);

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #506\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test507(void)
{
   int   status;
   int   ip1 = -1; 
   int   ip2 = -1;
   int   ip3 = -1;
   float p1 = 836.0; int kind3 =  4;
   float p2 =  12.0; int kind2 = 10;
   float p3 = 124.0; int kind1 =  4;
   int   ip1toCmp = 75284896;   /* 836.0  */
   int   ip2toCmp = 176280768;  /*  12.0  */
   int   ip3toCmp = 74572896;   /* 124.0  */
 
   fprintf(stderr, "Test #507 : Encode IP1, IP2 et IP3 d'après p1=%f, kind1=%d, p2=%f, kind2=%d, p3=%f, kind3=%d\n",
           p1, kind1, p2, kind2, p3, kind3);

   status = ConvertPKtoIP(&ip1, &ip2, &ip3, p1, kind1, p2, kind2, p3, kind3);

   if ((ip1 != ip1toCmp) || (ip2 != ip2toCmp) || (ip3 != ip3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #507\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((ip1 != %d) || (ip2 != %d) || (ip3 != %d) || (status != %d))\033[0m\n", ip1toCmp, ip2toCmp, ip3toCmp, CONVERT_OK);
   }
 
   fprintf(stderr,"ip1= %d, ip2= %d, ip3= %d, status=%d\n\n", ip1, ip2, ip3, status);
}


void Test600(void)
{
   int   status;
   int   ip1 = 40994464;   /* 100.0  2 */
   int   ip2 = 176280768;  /*  12.0 10 */
   int   ip3 = 176520768;  /*  36.0 10 */
   float p1 = 0.0;
   float p2 = 0.0;
   float p3 = 0.0;
   int   kind1 = -1;
   int   kind2 = -1;
   int   kind3 = -1;
   float p1toCmp = 100.0;
   float p2toCmp = 36.0;
   float p3toCmp = 12.0;
   int   kind1toCmp = 2;
   int   kind2toCmp = 10;
   int   kind3toCmp = 10;
 
   fprintf(stderr, "Test #600 : ConvertIPtoPK d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = ConvertIPtoPK(&p1, &kind1, &p2, &kind2, &p3, &kind3, ip1, ip2, ip3);  

   if ((p1 != p1toCmp) || (kind1 != kind1toCmp) || (p2 != p2toCmp) || (kind2 != kind2toCmp) || (p3 != p3toCmp) || (kind3 != kind3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #600\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (status != %d))\033[0m\n", p1, p1toCmp, kind1, kind1toCmp, p2, p2toCmp, kind2, kind2toCmp, p3, p3toCmp, kind3, kind3toCmp, CONVERT_OK);
   }

   fprintf(stderr,"p1= %f %d, p2= %f %d, p3= %f %d, status=%d\n\n", p1, kind1, p2, kind2, p3, kind3, status);
}


void Test601(void)
{
   int   status;
   int   ip1 = 40994464;   /* 100.0  2 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 176280768;  /*  12.0 10 */
   float p1 = 0.0;
   float p2 = 0.0;
   float p3 = 0.0;
   int   kind1 = -1;
   int   kind2 = -1;
   int   kind3 = -1;
   float p1toCmp = 100.0;
   float p2toCmp = 36.0;
   float p3toCmp = 12.0;
   int   kind1toCmp = 2;
   int   kind2toCmp = 10;
   int   kind3toCmp = 10;
 
   fprintf(stderr, "Test #601 : ConvertIPtoPK d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = ConvertIPtoPK(&p1, &kind1, &p2, &kind2, &p3, &kind3, ip1, ip2, ip3);  

   if ((p1 != p1toCmp) || (kind1 != kind1toCmp) || (p2 != p2toCmp) || (kind2 != kind2toCmp) || (p3 != p3toCmp) || (kind3 != kind3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #601\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (status != %d))\033[0m\n", p1, p1toCmp, kind1, kind1toCmp, p2, p2toCmp, kind2, kind2toCmp, p3, p3toCmp, kind3, kind3toCmp, CONVERT_OK);
   }

   fprintf(stderr,"p1= %f %d, p2= %f %d, p3= %f %d, status=%d\n\n", p1, kind1, p2, kind2, p3, kind3, status);
}


void Test602(void)
{
   int   status;
   int   ip1 = 7665032;    /* 325.0  0 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 8124032 ;   /* 784.0  0 */
   float p1 = 0.0;
   float p2 = 0.0;
   float p3 = 0.0;
   int   kind1 = -1;
   int   kind2 = -1;
   int   kind3 = -1;
   float p1toCmp = 325.0;
   float p2toCmp = 36.0;
   float p3toCmp = 784.0;
   int   kind1toCmp = 0;
   int   kind2toCmp = 10;
   int   kind3toCmp = 0;
 
   fprintf(stderr, "Test #602 : ConvertIPtoPK d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = ConvertIPtoPK(&p1, &kind1, &p2, &kind2, &p3, &kind3, ip1, ip2, ip3);  

   if ((p1 != p1toCmp) || (kind1 != kind1toCmp) || (p2 != p2toCmp) || (kind2 != kind2toCmp) || (p3 != p3toCmp) || (kind3 != kind3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #602\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (status != %d)\033[0m\n", p1, p1toCmp, kind1, kind1toCmp, p2, p2toCmp, kind2, kind2toCmp, p3, p3toCmp, kind3, kind3toCmp, CONVERT_OK);
   }

   fprintf(stderr,"p1= %f %d, p2= %f %d, p3= %f %d, status=%d\n\n", p1, kind1, p2, kind2, p3, kind3, status);
}


void Test603(void)
{
   int   status;
   int   ip1 = 8124032 ;   /* 784.0  0 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 7665032;    /* 325.0  0 */
   float p1 = 0.0;
   float p2 = 0.0;
   float p3 = 0.0;
   int   kind1 = -1;
   int   kind2 = -1;
   int   kind3 = -1;
   float p1toCmp = 325.0;
   float p2toCmp = 36.0;
   float p3toCmp = 784.0;
   int   kind1toCmp = 0;
   int   kind2toCmp = 10;
   int   kind3toCmp = 0;
 
   fprintf(stderr, "Test #603 : ConvertIPtoPK d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = ConvertIPtoPK(&p1, &kind1, &p2, &kind2, &p3, &kind3, ip1, ip2, ip3);  

   if ((p1 != p1toCmp) || (kind1 != kind1toCmp) || (p2 != p2toCmp) || (kind2 != kind2toCmp) || (p3 != p3toCmp) || (kind3 != kind3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #603\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (status != %d)\033[0m\n", p1, p1toCmp, kind1, kind1toCmp, p2, p2toCmp, kind2, kind2toCmp, p3, p3toCmp, kind3, kind3toCmp, CONVERT_OK);
   }

   fprintf(stderr,"p1= %f %d, p2= %f %d, p3= %f %d, status=%d\n\n", p1, kind1, p2, kind2, p3, kind3, status);
}


void Test604(void)
{
   int   status;
   int   ip1 = 27962976;   /*   0.7  1 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 27362976;   /*   0.1  1 */
   float p1 = 0.0;
   float p2 = 0.0;
   float p3 = 0.0;
   int   kind1 = -1;
   int   kind2 = -1;
   int   kind3 = -1;
   float p1toCmp = 0.7;
   float p2toCmp = 36.0;
   float p3toCmp = 0.1;
   int   kind1toCmp = 1;
   int   kind2toCmp = 10;
   int   kind3toCmp = 1;
 
   fprintf(stderr, "Test #604 : ConvertIPtoPK d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = ConvertIPtoPK(&p1, &kind1, &p2, &kind2, &p3, &kind3, ip1, ip2, ip3);  

   if ((p1 != p1toCmp) || (kind1 != kind1toCmp) || (p2 != p2toCmp) || (kind2 != kind2toCmp) || (p3 != p3toCmp) || (kind3 != kind3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #604\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (status != %d))\033[0m\n", p1, p1toCmp, kind1, kind1toCmp, p2, p2toCmp, kind2, kind2toCmp, p3, p3toCmp, kind3, kind3toCmp, CONVERT_OK);
   }

   fprintf(stderr,"p1= %f %d, p2= %f %d, p3= %f %d, status=%d\n\n", p1, kind1, p2, kind2, p3, kind3, status);
}


void Test605(void)
{
   int   status;
   int   ip1 = 27362976;   /*   0.1  1 */
   int   ip2 = 176520768;  /*  36.0 10 */
   int   ip3 = 27962976;   /*   0.7  1 */
   float p1 = 0.0;
   float p2 = 0.0;
   float p3 = 0.0;
   int   kind1 = -1;
   int   kind2 = -1;
   int   kind3 = -1;
   float p1toCmp = 0.7;
   float p2toCmp = 36.0;
   float p3toCmp = 0.1;
   int   kind1toCmp = 1;
   int   kind2toCmp = 10;
   int   kind3toCmp = 1;
 
   fprintf(stderr, "Test #605 : ConvertIPtoPK d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = ConvertIPtoPK(&p1, &kind1, &p2, &kind2, &p3, &kind3, ip1, ip2, ip3);  

   if ((p1 != p1toCmp) || (kind1 != kind1toCmp) || (p2 != p2toCmp) || (kind2 != kind2toCmp) || (p3 != p3toCmp) || (kind3 != kind3toCmp) || (status != CONVERT_OK))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #605\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (status != %d))\033[0m\n", p1, p1toCmp, kind1, kind1toCmp, p2, p2toCmp, kind2, kind2toCmp, p3, p3toCmp, kind3, kind3toCmp, CONVERT_OK);
   }

   fprintf(stderr,"p1= %f %d, p2= %f %d, p3= %f %d, status=%d\n\n", p1, kind1, p2, kind2, p3, kind3, status);
}


void Test606(void) /* time/level/level reordered to level/time/level with warning*/
{
   int   status;
   int   ip1 = 176520768;  /*  36.0 10 */
   int   ip2 = 27962976;   /*   0.7  1 */
   int   ip3 = 27362976;   /*   0.1  1 */
   float p1 = 0.0;
   float p2 = 0.0;
   float p3 = 0.0;
   int   kind1 = -1;
   int   kind2 = -1;
   int   kind3 = -1;
   float p2toCmp = 36.0;
   float p1toCmp = 0.7;
   float p3toCmp = 0.1;
   int   kind2toCmp = 10;
   int   kind1toCmp = 1;
   int   kind3toCmp = 1;
 
   fprintf(stderr, "Test #606 : ConvertIPtoPK d'après ip1=%d, ip2=%d, ip3=%d\n", ip1, ip2, ip3);

   status = ConvertIPtoPK(&p1, &kind1, &p2, &kind2, &p3, &kind3, ip1, ip2, ip3);  

   if ((p1 != p1toCmp) || (kind1 != kind1toCmp) || (p2 != p2toCmp) || (kind2 != kind2toCmp) || (p3 != p3toCmp) || (kind3 != kind3toCmp) || (status != CONVERT_WARNING))
   { 
      fprintf(stderr, "\033[31;1mErreur dans le test #606\033[0m\n");
      fprintf(stderr, "\033[31;1mif ((%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (%f != %f) || (%d != %d) || (status != %d))\033[0m\n", p1, p1toCmp, kind1, kind1toCmp, p2, p2toCmp, kind2, kind2toCmp, p3, p3toCmp, kind3, kind3toCmp, CONVERT_WARNING);
   }

   fprintf(stderr,"p1= %f %d, p2= %f %d, p3= %f %d, status=%d\n\n", p1, kind1, p2, kind2, p3, kind3, status);
}


void Test700(void)
{
}


void Test800(void)
{
}



int Cmain2()  /* called by Fortran test */
{
   fprintf(stderr, "\n================= Test avec la fonction ConvIp =================\n\n");
   Test1();
   Test2();
   Test3();
   Test4();

   fprintf(stderr, "\n================= Test avec la fonction EncodeIp =================\n\n");
   Test100();
   Test101();
   Test102();
   Test103();
   Test104();
   Test105();
   Test106();
   Test107();
   Test108();
   Test109();
   Test110();

   fprintf(stderr, "\n================= Test avec la fonction DecodeIp =================\n\n");  
   Test200();
   Test201();
   Test202();
   Test203();
   Test204();
   Test205();
   Test206();
   Test207();

   fprintf(stderr, "\n================= Test avec la fonction EncodeIp_v =================\n\n");
   Test300();

   fprintf(stderr, "\n================= Test avec la fonction DecodeIp_v =================\n\n");
   Test400();

   fprintf(stderr, "\n================= Test avec la fonction ConvertPKtoIp =================\n\n");
   Test500();
   Test501();
   Test502();
   Test503();
   Test504();
   Test505();
   Test506();
   Test507();

   fprintf(stderr, "\n================= Test avec la fonction ConvertIPtoPKP =================\n\n");
   Test600();
   Test601();
   Test602();
   Test603();
   Test604();
   Test605();
   Test606();

   fprintf(stderr, "\n================= Test avec la fonction ConvertPKtoIP_v =================\n\n");
   Test700();

   fprintf(stderr, "\n================= Test avec la fonction ConvertIPtoPKP_v =================\n\n");
   Test800();

   return 0;
}