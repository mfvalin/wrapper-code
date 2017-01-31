
  type, BIND(C) :: string
    type(C_PTR) :: ptr
  end type string
  type(string), parameter :: NULL_STRING = string(C_NULL_PTR)

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

