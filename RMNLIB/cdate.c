/*
 *
 *  file      :  CDATE.C
 *
 *  author    :  Marc Gagnon  92/03/13
 *
 *  revision  :  V1.0 Michel Grenier; conversion a l'an 2000
 *                                    utilisation de rmnlib.
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  object    :  THIS FILE CONTAINS SOME OF THE MODULES
 *               MANAGING DATES
 *
 *
 */

#define CCDATE

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "cdate.h"
#include "rpnmacros.h"

extern int *GetFmain;

/*
 *  reference date: Jan 1 1950 00:00:00:00
 */

 static int refdate  = 19500101;
 static int reftime  = 0;

/*
 *
 *  module    :  CD_JS_CO
 *
 *  author    :  Marc Gagnon
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  call      :  int js;
 *               js = cd_js_co();
 *
 *  object    :  THIS MODULE RETURNS THE CURRENT DAY OF THE WEEK.
 *
 */

 extern int
 cd_js_co (void)
    {
    struct tm *tm;
    time_t temps;

    temps = time(NULL);
    tm = gmtime(&temps);

    return tm->tm_wday;  /* 0=dimanche */
    }

/*
 *
 *  module    :  CD_HR_CO
 *
 *  author    :  Marc Gagnon
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  call      :  int hr;
 *               hr = cd_hr_co();
 *
 *  object    :  THIS MODULE RETURNS THE CURRENT HOUR.
 *
 */

 extern int
 cd_hr_co (void)
    {
    struct tm *tm;
    time_t temps;

    temps = time(NULL);
    tm = gmtime(&temps);

    return tm->tm_hour;  /* 0-23 */
    }

/*
 *
 *  module    :  CD_JJ_CO
 *
 *  author    :  Marc Gagnon
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  call      :  int jj;
 *               jj = cd_jj_co();
 *
 *  object    :  THIS MODULE RETURNS THE CURRENT JULIAN DAY
 *
 */

 extern int
 cd_jj_co (void)
    {
    struct tm *tm;
    time_t temps;
    int aa,mm,jj;

    temps = time(NULL);
    tm = gmtime(&temps);
    aa = tm->tm_year + 1900;
    mm = tm->tm_mon + 1;
    jj = tm->tm_mday;

    return (cd_gr_jj(aa,mm,jj));
    }

/*
 *
 *  module    :  CD_GR_CO
 *
 *  author    :  Marc Gagnon
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  call      :  int a,m,j;
 *               cd_gr_co(&a,&m,&j);
 *
 *  object    :  THIS MODULE RETURNS THE CURRENT GREGORIAN DATE
 *
 */

 extern void
 cd_gr_co ( int *a, int *m, int *j )
    {
    struct tm *tm;
    time_t temps;

    temps = time(NULL);
    tm = gmtime(&temps);
    *a = tm->tm_year + 1900;
    *m = tm->tm_mon + 1;
    *j = tm->tm_mday;
    }

/*
 *
 *  module    :  CD_GR_JJ
 *
 *  author    :  Marc Gagnon
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  call      :  int a,m,j;
 *               int jj;
 *               jj = cd_gr_jj(aa,mm,jj);
 *
 *  object    :  THIS MODULE CONVERTS A GREGORIAN DATE TO DAYS SINCE REFERENCE
 *
 */

 extern int
 cd_gr_jj ( int aa, int mm, int jj )
    {
    double jhr;

    int dat;
    int mode;
    int refstamp;
    int stamp;
    int tim;
    int jday;

/*
 *  reference date/time to reference stamp
 */
    mode = 3;
    f77name(newdate)( &refstamp, &refdate, &reftime, &mode );

/*
 *  given date/time to stamp
 */
    dat  = (aa*100+mm)*100+jj;
    tim  = 0;
    mode = 3;
    f77name(newdate)( &stamp, &dat, &tim, &mode );

/*
 *  julian hours = difference between date and refdate
 */
    f77name(difdatr)( &stamp, &refstamp, &jhr );

/*
 *  converted to julian days
 */
    jhr += 0.1;
    jday = (int)(jhr/24.0);

    return (jday);
    }

/*
 *
 *  module    :  CD_JJ_GR
 *
 *  author    :  Marc Gagnon
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  call      :  int *a,*m,*j;
 *               int jul;
 *               cd_jj_gr(a,m,j,jul);
 *
 *  object    :  THIS MODULE CONVERTS A DAYS SINCE REFERENCE TO A GREGORIAN DATE
 *
 */

 extern void
 cd_jj_gr ( int *pyyyy, int *pmm, int *pdd, int jul )
    {
    double jhr;
    int dat;
    int mode;
    int refstamp;
    int stamp = 0;
    int tim;

/*
 *  reference date/time to reference stamp
 */
    mode = 3;
    f77name(newdate)( &refstamp, &refdate, &reftime, &mode );

/*
 *  increment reference time by julian hours
 */
    jhr = jul*24.0;
    f77name(incdatr)( &stamp, &refstamp, &jhr );

/*
 *  resulting stamp to printable
 */
    mode = -3;
    f77name(newdate)( &stamp, &dat, &tim, &mode );

/*
 *  parse the resulting date
 */
    *pyyyy =  dat/10000;
    *pmm   = (dat%10000)/100;
    *pdd   =  dat%100;
    }

/*
 *
 *  module    :  LEAP
 *
 *  author    :  Marc Gagnon
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  call      :  int yyyy;
 *               leap(yyyy)
 *
 *  object    :  THIS MODULE RETURNS IF THE GIVEN YEAR IS A LEAP YEAR
 *
 */

 static int
 leap(int yyyy)
    {
    int leapr,a,b,c;

    a = (yyyy%4   == 0);
    b = (yyyy%100 != 0);
    c = (yyyy%400 == 0);

    leapr = ( (a&&b)||c ? 1 : 0 );

    return( leapr );
    }

extern int
#ifdef Mop_linux
main_cdate (argc,argv)
#else
main (argc,argv)
#endif
     int argc;
     char *argv[];
    {
    int a,m,j,h,jul;

/*
 *  test for argument error
 */
    if( argc < 2 ) { usage(); }

/*
 *  argument -jg : conversion GMT date/time from julian hour to gregorian
 */
    if( strcmp(argv[1],"-jg") == 0 )
      {
      if( argc != 3 ) { usage(); }

      jul = atoi(argv[2]);
      cd_jj_gr(&a,&m,&j,jul/24);
      printf("%04d %02d %02d %02d\n",a,m,j,jul % 24);
      exit(0);
      }

/*
 *  argument -gj : conversion GMT date/time from gregorian to julian hour
 */
    if( strcmp(argv[1],"-gj") == 0 )
      {
      if( argc != 6 ) { usage(); }
      a = atoi(argv[2]);
      m = atoi(argv[3]);
      j = atoi(argv[4]);
      h = atoi(argv[5]);
      printf("%d\n",cd_gr_jj(a,m,j)*24+h);
      exit(0);
      }

/*
 *  argument -jc : current GMT date/time in julian hour
 */
    if( strcmp(argv[1],"-jc") == 0 )
      {
      printf("%d\n",cd_jj_co()*24+cd_hr_co());
      exit(0);
      }

/*
 *  argument -gc : current GMT date/time (gregorian)
 */
    if( strcmp(argv[1],"-gc") == 0 )
      {
      h = cd_hr_co();
      cd_gr_co(&a,&m,&j);
      printf("%04d %02d %02d %02d\n",a,m,j,h);
      exit(0);
      }

    usage();
    }

/*
 *
 *  module    :  USAGE
 *
 *  author    :  Marc Gagnon
 *
 *  revision  :  V0.0
 *
 *  status    :  DEVELOPMENT
 *
 *  language  :  C
 *
 *  call      :  usage();
 *
 *  object    :  THIS MODULE PRINTS ON STANDARD OUTPUT THE USAGE OF
 *               PROGRAM CDATE
 *
 */

 extern int
 usage ( void )
    {
    char *p;

    p = getenv("CMCLNG");

    if( p == NULL || p[0] == 'e' || p[0] == 'E' )
      {
      printf("usage: cdate -gj AAAA MM DD HH | -jg NNNNNNNN | -jc | -gc\n");
      printf(" -gj conversion from gregorian date to julian hour.\n");
      printf(" -jg conversion from julian hour to gregorian date.\n");
      printf(" -jc display current julian hour (GMT).\n");
      printf(" -gc display current AAAA MM DD HH (GMT).\n");
      exit(1);
      }

    printf("usage: cdate -gj AAAA MM JJ HH | -jg NNNNNNNN | -jc | -gc\n");
    printf(" -gj conversion date gregorienne en heure julienne.\n");
    printf(" -jg conversion heure julienne en date gregorienne.\n");
    printf(" -jc afficher l'heure julienne courante (GMT).\n");
    printf(" -gc afficher AAAA MM JJ HH courant (GMT).\n");
    exit(1);
    }
