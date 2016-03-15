/*----------------------------------------------------------------------
 *    FICHIER    cdate.h
 *
 *    BUT        Conversions de date.
 *
 *    AUTEUR     Marc Gagnon [92/03/16].
 *
 *    ALGORITHME Revue "Communication of the ACM"
 *
 *    REVISION
 *    $Header: /home/binops/afsi/sio/cvs_backup/programs/cdate/cdate.h,v 1.7 2003/05/28 18:27:37 afsisio Exp $
\*----------------------------------------------------------------------*/

#ifndef HCDATE
#define HCDATE

 int  cd_hr_co(void);
 int  cd_jj_co(void);
 int  cd_js_co(void);

 void cd_gr_co(int *a,int *m,int *j);
 int  cd_gr_jj(int  a,int  m,int  j);
 void cd_jj_gr(int *a,int *m,int *j,int jj);

#endif
