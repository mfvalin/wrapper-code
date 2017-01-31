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
  interface operator(+)                                          !InTf!
    function str_plus_str(str1,str2) BIND(C, name='StrPlusStr') result(str3)   !InTf!
      import :: string                                           !InTf!
      type(string), intent(IN), value :: str1,str2               !InTf!
      type(string) :: str3                                       !InTf!
    end function str_plus_str                                    !InTf!
    function chr_plus_str(chr, str)  result(str2) BIND(C, name='ChrPlusStr')  !InTf!
      import :: C_CHAR, string                                   !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr         !InTf!
      type(string), intent(IN), value :: str                     !InTf!
      type(string) :: str2                                       !InTf!
    end function chr_plus_str                                    !InTf!
    function str_plus_chr(str, chr)  result(str2) BIND(C, name='StrPlusChr')  !InTf!
      import :: C_CHAR, string                                   !InTf!
      type(string), intent(IN), value :: str                     !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr         !InTf!
      type(string) :: str2                                       !InTf!
    end function str_plus_chr                                    !InTf!
    function chr_plus_chr(chr1, chr2)  result(str1) BIND(C, name='ChrPlusChr')  !InTf!
      import :: C_CHAR, string                                   !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr1, chr2  !InTf!
      type(string) :: str1                                       !InTf!
    end function chr_plus_chr                                    !InTf!
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
