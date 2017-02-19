module fstring_internals
  use ISO_C_BINDING

  type, BIND(C) :: F_STRNG
    type(C_PTR) :: ptr
  end type F_STRNG
  type(F_STRNG), parameter :: NULL_F_STRNG = F_STRNG(C_NULL_PTR)

  interface operator(+)
    module procedure str_plus_str
    module procedure chr_plus_str
    module procedure str_plus_chr
    module procedure chr_plus_chr
  end interface

  interface assignment(=)
    subroutine str_to_str(str, str2) BIND(C, name='StrToStr')  ! F_STRNG = F_STRNG
      import :: F_STRNG
      type(F_STRNG), intent(OUT) :: str
      type(F_STRNG), intent(IN) :: str2
    end subroutine str_to_str
    subroutine chr_to_str(str, chr) BIND(C, name='ChrToStr')   ! F_STRNG = C_STRING
      import :: C_CHAR, F_STRNG
      type(F_STRNG), intent(OUT) :: str
      character(C_CHAR), dimension(*), intent(IN) :: chr
    end subroutine chr_to_str
    module procedure chr_to_fstr                               ! fortran_string = C_STRING
    module procedure str_to_fstr                               ! fortran_string = F_STRNG
  end interface

  interface str_to_chr
    subroutine str_to_chr(chr, str, n) BIND(C, name='StrToChr')
      import :: C_CHAR, F_STRNG, C_INT
      type(F_STRNG), intent(IN) :: str
      character(C_CHAR), dimension(*), intent(OUT) :: chr
      integer(C_INT), intent(IN) :: n
    end subroutine str_to_chr
  end interface

  interface fstr_addition
    subroutine strstr(stro,str1,str2) BIND(C, name='F_StrPlusStr')
      import :: F_STRNG
      type(F_STRNG), intent(OUT) :: stro
      type(F_STRNG), intent(IN)  :: str1,str2
    end subroutine strstr

    subroutine chrstr(stro, chr, str) BIND(C, name='F_ChrPlusStr')
      import :: C_CHAR, F_STRNG
      type(F_STRNG), intent(OUT) :: stro
      character(C_CHAR), dimension(*), intent(IN) :: chr
      type(F_STRNG), intent(IN)  :: str
    end subroutine chrstr

    subroutine chrchr(stro, chr1, chr2) BIND(C, name='F_ChrPlusChr')
      import :: C_CHAR, F_STRNG
      type(F_STRNG), intent(OUT) :: stro
      character(C_CHAR), dimension(*), intent(IN) :: chr1, chr2
    end subroutine chrchr
  end interface
contains
  function str_plus_str(str1,str2) BIND(C, name='F90StrPlusStr') result(stro)
    implicit none
    type(F_STRNG), intent(IN)        :: str1,str2
    type(F_STRNG) :: stro
    call fstr_addition(stro,str1,str2)
  end function str_plus_str

  function chr_plus_str(chr, str)  result(stro) BIND(C, name='F90ChrPlusStr')
    implicit none
    character(C_CHAR), dimension(*), intent(IN) :: chr
    type(F_STRNG), intent(IN)        :: str
    type(F_STRNG) :: stro
    call fstr_addition(stro,chr,str)
  end function chr_plus_str

  function str_plus_chr(str, chr)  result(stro) BIND(C, name='F90StrPlusChr')
    implicit none
    type(F_STRNG), intent(IN)        :: str
    character(C_CHAR), dimension(*), intent(IN) :: chr
    type(F_STRNG) :: stro
    call fstr_addition(stro,chr,str)
  end function str_plus_chr

  function chr_plus_chr(chr1, chr2)  result(stro) BIND(C, name='F90ChrPlusChr')
    implicit none
    character(C_CHAR), dimension(*), intent(IN) :: chr1, chr2
    type(F_STRNG) :: stro
    call fstr_addition(stro,chr1,chr2)
  end function chr_plus_chr

  subroutine chr_to_fstr(fstr, chr)   ! "fortran F_STRNG" = array of char
    implicit none
    character(len=*), intent(OUT) :: fstr
    character(C_CHAR), dimension(*), intent(IN) :: chr
    integer :: i, lfstr

    lfstr = len(fstr)
    fstr = " "
    do i=1,lfstr
      if(ichar(chr(i)) == 0) exit
      fstr(i:i) = chr(i)
    enddo
    return
  end subroutine chr_to_fstr

  subroutine str_to_fstr(fstr, str)   ! "fortran string" = F_STRNG
    implicit none
    character(len=*), intent(OUT) :: fstr
    type(F_STRNG), intent(IN), value :: str
    integer :: i, lfstr
    character(C_CHAR), dimension(:), pointer :: chr

    lfstr = len(fstr)
    fstr = " "
    call c_f_pointer(str%ptr,chr,[lfstr])
    do i=1,lfstr
      if(ichar(chr(i)) == 0) exit
      fstr(i:i) = chr(i)
    enddo
    return
  end subroutine str_to_fstr

  subroutine fstr_to_str(str, fstr, n)    ! F_STRNG = "fortran string"
    implicit none
    type(F_STRNG), intent(OUT) :: str
    character(len=*), intent(IN) :: fstr
    integer(C_INT), intent(IN), value :: n
    integer :: i, lfstr
    character(C_CHAR), dimension(:), pointer :: chr

    lfstr = min(len_trim(fstr),n-1)
    call c_f_pointer(str%ptr,chr,[lfstr])
    do i=1,lfstr
      chr(i) = fstr(i:i)
    enddo
    chr(lfstr+1) = achar(0)
    return
  end subroutine fstr_to_str

end module fstring_internals
