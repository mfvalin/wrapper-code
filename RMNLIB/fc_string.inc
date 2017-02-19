
  type, BIND(C) :: F_STRNG
    type(C_PTR) :: ptr
  end type F_STRNG
  type(F_STRNG), parameter :: NULL_F_STRNG = F_STRNG(C_NULL_PTR)

  interface operator(+)                                          !InTf!
    function str_plus_str(str1,str2) BIND(C, name='F90StrPlusStr') result(str3)   !InTf!
      import :: F_STRNG                                           !InTf!
      type(F_STRNG), intent(IN)        :: str1,str2               !InTf!
      type(F_STRNG) :: str3                                       !InTf!
    end function str_plus_str                                    !InTf!
    function chr_plus_str(chr, str)  result(str2) BIND(C, name='F90ChrPlusStr')  !InTf!
      import :: C_CHAR, F_STRNG                                   !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr         !InTf!
      type(F_STRNG), intent(IN)        :: str                     !InTf!
      type(F_STRNG) :: str2                                       !InTf!
    end function chr_plus_str                                    !InTf!
    function str_plus_chr(str, chr)  result(str2) BIND(C, name='F90StrPlusChr')  !InTf!
      import :: C_CHAR, F_STRNG                                   !InTf!
      type(F_STRNG), intent(IN)        :: str                     !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr         !InTf!
      type(F_STRNG) :: str2                                       !InTf!
    end function str_plus_chr                                    !InTf!
    function chr_plus_chr(chr1, chr2)  result(str1) BIND(C, name='F90ChrPlusChr')  !InTf!
      import :: C_CHAR, F_STRNG                                   !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr1, chr2  !InTf!
      type(F_STRNG) :: str1                                       !InTf!
    end function chr_plus_chr                                    !InTf!
  end interface                                                  !InTf!

  interface str_to_chr                                           !InTf!
    subroutine str_to_chr(chr, str, n) BIND(C, name='StrToChr')  !InTf!
      import :: C_CHAR, F_STRNG, C_INT                            !InTf!
      type(F_STRNG), intent(IN) :: str                            !InTf!
      character(C_CHAR), dimension(*), intent(OUT) :: chr        !InTf!
      integer(C_INT), intent(IN) :: n                            !InTf!
    end subroutine str_to_chr                                    !InTf!
  end interface                                                  !InTf!

  interface assignment(=)                                        !InTf!
    subroutine str_to_str(str, str2) BIND(C, name='StrToStr')    !InTf!
      import :: F_STRNG                                           !InTf!
      type(F_STRNG), intent(OUT) :: str                           !InTf!
      type(F_STRNG), intent(IN) :: str2                           !InTf!
    end subroutine str_to_str                                    !InTf!
    subroutine chr_to_str(str, chr) BIND(C, name='ChrToStr')     !InTf!
      import :: C_CHAR, F_STRNG                                   !InTf!
      type(F_STRNG), intent(OUT) :: str                           !InTf!
      character(C_CHAR), dimension(*), intent(IN) :: chr         !InTf!
    end subroutine chr_to_str                                    !InTf!
  end interface                                                  !InTf!
