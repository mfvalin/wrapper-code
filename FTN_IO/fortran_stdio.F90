function fortran_fopen(path,mode) result(file)
!
! fortran interface to libc's fopen
! arguments follow man fopen
! needed to convert fortran strings to char *
!
  use ISO_C_BINDING
  implicit none
  character(len=*), intent(IN) :: path,mode
  type(C_PTR) :: file

  interface
    function f__fopen(path,mode) result(file) bind(C,name='fopen')
      import
      implicit none
      type(C_PTR) :: file
      character(C_CHAR), intent(IN), dimension(*) :: path, mode
    end function f__fopen
  end interface

  character(len=1), dimension(:), allocatable :: lpath, lmode

  allocate(lpath(1+len_trim(path)), lmode(1+len_trim(mode)) )
  lpath = transfer(trim(path)//achar(0),lpath)
  lmode = transfer(trim(mode)//achar(0),lmode)

  file = f__fopen(lpath , lmode )
  deallocate(lpath , lmode )

end function fortran_fopen
#if defined(SELF_TEST)
! s.f90 -DSELF_TEST fortran_stdio.F90
! ./a.out
program test_fortran_stdio
  use ISO_C_BINDING
  include 'fortran_stdio.inc'
  type(C_PTR) :: file
  integer, dimension(100), target :: buffer
  integer(C_SIZE_T) :: status, nelem, elemsize
  integer :: i
  integer(C_LONG) :: where

  print *,'calling fortran_fopen'
  file = fortran_fopen('tagada','w+')
  if(.not. C_ASSOCIATED(file) ) then
    print *,'open failed'
    stop
  endif

  nelem = 100
  elemsize = 4
  do i = 1,nelem
    buffer(i) = i
  enddo

  print *,'calling fortran_fwrite'
  where = fortran_ftell(file)
  print *,'Where=',where
  status = fortran_fwrite(c_loc(buffer(1)) , elemsize , nelem , file)
  print *,'status=',status
  where = fortran_ftell(file)
  print *,'Where=',where

  print *,'calling fortran_fclose'
  status = fortran_fclose(file)

  print *,'calling fortran_fopen'
  file = fortran_fopen('tagada','r')
  if(.not. C_ASSOCIATED(file) ) then
    print *,'open failed'
    stop
  endif

  buffer = -1
  print *,'calling fortran_fread'
  where = fortran_ftell(file)
  print *,'Where=',where
  status = fortran_fread(c_loc(buffer(1)) , elemsize , nelem , file)
  where = fortran_ftell(file)
  print *,'Where=',where
  print *,'status=',status
  
  print *,buffer(1:status)

  print *,'calling fortran_fseek'
  where = -200
  status = fortran_fseek(file,where,SEEK_CUR)
  print *,'status=',status
  where = fortran_ftell(file)
  print *,'Where=',where
  print *,'calling fortran_fread'
  status = fortran_fread(c_loc(buffer(1)) , elemsize , nelem , file)
  where = fortran_ftell(file)
  print *,'Where=',where
  print *,'status=',status
  
  print *,buffer(1:status)

  print *,'calling fortran_fclose'
  status = fortran_fclose(file)

end program
#endif