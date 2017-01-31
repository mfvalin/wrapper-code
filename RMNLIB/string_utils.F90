subroutine test_string_utils
  use ISO_C_BINDING
  implicit none
  include 'fc_string.inc'
#include "sizeof.hf"

  interface assignment(=)
    ! "fortran string" = character*1 array (null terminated)
    subroutine chr_to_fstr(fstr, chr)
      import :: C_CHAR
      character(len=*), intent(OUT) :: fstr
      character(C_CHAR), dimension(*), intent(IN) :: chr
    end subroutine chr_to_fstr
    subroutine str_to_fstr(fstr, str)   ! "fortran string" = string
      import :: string
      character(len=*), intent(OUT) :: fstr
      type(string), intent(IN), value :: str
    end subroutine str_to_fstr
  end interface

!   interface
!     subroutine fstr_to_str(str, fstr, n)    ! string = "fortran string"
!       import :: string
!       type(string), intent(OUT) :: str
!       character(len=*), intent(IN) :: fstr
!   integer(C_INT), intent(IN), value :: n
!     end subroutine fstr_to_str
!   end interface

!   interface assignment(=)
!     subroutine chr_to_chr(chr, chr2) BIND(C)  ! character*1 array (null terminated) = character*1 array (null terminated)
!       import :: C_CHAR
!       character(C_CHAR), dimension(*), intent(OUT) :: chr
!       character(C_CHAR), dimension(*), intent(IN) :: chr2
!     end subroutine chr_to_chr
!     subroutine fstr_to_chr(chr, fstr)    ! character*1 array (null terminated) = "fortran string"
!       import :: C_CHAR
!       character(C_CHAR), dimension(*), intent(OUT) :: chr
!       character(len=*), intent(IN) :: fstr
!     end subroutine fstr_to_chr
!   end interface

  character(C_CHAR), dimension(100) :: chr1, chr2
  character(C_CHAR), dimension(:), allocatable :: c
!   character(C_CHAR), dimension(1) :: c
  type(string) :: mystr1, mystr2, mystr3
  character(len=64) :: s64

  mystr1 = chr1
!   chr2 = mystr1
  call str_to_chr(chr2, mystr1, size(chr2))
  mystr2 = chr2
  mystr3 = mystr1 + mystr2 + chr1
  mystr2 = mystr2 + CSTRING("bla bla bla!!")
  mystr3 = mystr3 + transfer("bla bla bla!!"//achar(0),c)
  mystr1 = mystr3
  s64 = transfer(chr2,s64)

end subroutine test_string_utils

subroutine chr_to_fstr(fstr, chr)   ! "fortran string" = array of char
  use ISO_C_BINDING
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

subroutine str_to_fstr(fstr, str)   ! "fortran string" = string
  use ISO_C_BINDING
  implicit none
  include 'fc_string.inc'
  character(len=*), intent(OUT) :: fstr
  type(string), intent(IN), value :: str
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

subroutine fstr_to_str(str, fstr, n)    ! string = "fortran string"
  use ISO_C_BINDING
  implicit none
  include 'fc_string.inc'
  type(string), intent(OUT) :: str
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
