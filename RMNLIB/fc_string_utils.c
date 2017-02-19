#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <malloc.h>

typedef unsigned char uchar ;
typedef struct{
  char *p;
}fc_string ;

/*
  type, BIND(C) :: string                                        !InTf!
    type(C_PTR) :: ptr                                           !InTf!
  end type string                                                !InTf!
  type(string), parameter :: NULL_STRING = string(C_NULL_PTR)    !InTf!
*/
/*
  interface str_addition                                         !InTf!
    subroutine str_plus_str(str3,str1,str2) BIND(C, name='F_StrPlusStr')   !InTf!
      import :: string                                           !InTf!
      type(string), intent(IN)        :: str1,str2               !InTf!
      type(string), intent(OUT) :: str3                          !InTf!
    end subroutine str_plus_str                                  !InTf!
    subroutine chr_plus_str(str2, chr, str) BIND(C, name='F_ChrPlusStr')  !InTf!
      import :: C_CHAR, string                                   !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr         !InTf!
      type(string), intent(IN)        :: str                     !InTf!
      type(string), intent(OUT) :: str2                          !InTf!
    end subroutine chr_plus_str                                  !InTf!
    subroutine str_plus_chr(str2, str, chr) BIND(C, name='F_StrPlusChr')  !InTf!
      import :: C_CHAR, string                                   !InTf!
      type(string), intent(IN)        :: str                     !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr         !InTf!
      type(string), intent(OUT) :: str2                          !InTf!
    end subroutine str_plus_chr                                  !InTf!
    subroutine chr_plus_chr(str1, chr1, chr2) BIND(C, name='F_ChrPlusChr')  !InTf!
      import :: C_CHAR, string                                   !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr1, chr2  !InTf!
      type(string), intent(OUT) :: str1                          !InTf!
    end subroutine chr_plus_chr                                  !InTf!
  end interface                                                  !InTf!
*/
#pragma weak ChrPlusChr=StrPlusStr
#pragma weak ChrPlusStr=StrPlusStr
#pragma weak StrPlusChr=StrPlusStr
fc_string ChrPlusChr(uchar *chr1, uchar *chr2);
fc_string StrPlusChr(fc_string str1, uchar *chr1);
fc_string ChrPlusStr(uchar *chr1, fc_string str1);
fc_string StrPlusStr(fc_string str1, fc_string str2){
  size_t len1, len2 ;
  fc_string result;

  len1 = (str1.p) ? strlen(str1.p) : 0 ;
  len2 = (str2.p) ? strlen(str2.p) : 0 ;

  if(len1 + len2 > 0) {
    result.p = (uchar *) malloc(len1 + len2 + 1) ;
    if(str1.p) strncpy(result.p       , str1.p, len1);    // copy str1 if not null
    if(str2.p) strncpy(result.p + len1, str2.p, len2);    // add str2 if not null
    result.p[len1+len2] = '\0';                           // null terminate
    if(str1.p) free(str1.p);
  }else{
    result.p = NULL ;
  }
  return(result);
}
// Fortran interface has to pass derived types by reference, 
// as pass by value is not reliably portable
#pragma weak F_ChrPlusChr=F_StrPlusStr
#pragma weak F_ChrPlusStr=F_StrPlusStr
#pragma weak F_StrPlusChr=F_StrPlusStr
void F_ChrPlusChr(fc_string *str3, fc_string *str1, fc_string *str2);
void F_StrPlusChr(fc_string *str3, fc_string *str1, fc_string *str2);
void F_ChrPlusStr(fc_string *str3, fc_string *str1, fc_string *str2);
void F_StrPlusStr(fc_string *str3, fc_string *str1, fc_string *str2)
{
  *str3 = StrPlusStr(*str1, *str2);
}

/*
  interface str_to_chr                                           !InTf!
    subroutine str_to_chr(chr, str, n) BIND(C, name='StrToChr')  !InTf!
      import :: C_CHAR, string, C_INT                            !InTf!
      type(string), intent(IN) :: str                            !InTf!
      character(C_CHAR), dimension(*), intent(OUT) :: chr        !InTf!
      integer(C_INT), intent(IN) :: n                            !InTf!
    end subroutine str_to_chr                                    !InTf!
  end interface                                                  !InTf!

  interface assignment(=)                                        !InTf!
    subroutine str_to_str(str, str2) BIND(C, name='StrToStr')    !InTf!
      import :: string                                           !InTf!
      type(string), intent(OUT) :: str                           !InTf!
      type(string), intent(IN) :: str2                           !InTf!
    end subroutine str_to_str                                    !InTf!
    subroutine chr_to_str(str, chr) BIND(C, name='ChrToStr')     !InTf!
      import :: C_CHAR, string                                   !InTf!
      type(string), intent(OUT) :: str                           !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr         !InTf!
    end subroutine chr_to_str                                    !InTf!
  end interface                                                  !InTf!
*/

void StrToChr(uchar *chr, fc_string str, int n){  // copy string into array of bytes
  if(str.p) {
    if(strlen(str.p) < n){
      strncpy(chr, str.p, strlen(str.p));
    }else{
      strncpy(chr, str.p, n-1);
    }
  }
}

#pragma weak ChrToStr=StrToStr
void ChrToStr(void *str, void *str2);
void StrToStr(fc_string *str, fc_string str2){
  size_t len2 = 0;

  if(str->p) free(str->p) ;    // get rid of original string if it exists
  if(str2.p) {                            // string 2 exists
    len2 = strlen(str2.p);
    str->p = (uchar *) malloc(len2 + 1);
    strncpy(str->p, str2.p, len2);        // copy string 2
  }else{                                  // string 2 is null
    str->p = NULL;
  }
}
