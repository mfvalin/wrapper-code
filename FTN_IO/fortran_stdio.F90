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

function fortran_fdopen(fd,mode) result(file)  ! not very useful
  use ISO_C_BINDING
  implicit none
  character(len=*), intent(IN) :: mode
  integer(C_INT), intent(IN) :: fd
  type(C_PTR) :: file

  interface
    function f__fdopen(fd,mode) result(file) bind(C,name='fdopen')
      import
      implicit none
      type(C_PTR) :: file
      integer(C_INT), intent(IN) :: fd
      character(C_CHAR), intent(IN), dimension(*) :: mode
    end function f__fdopen
  end interface

  character(len=1), dimension(:), allocatable :: lmode

  allocate( lmode(1+len_trim(mode)) )
  lmode = transfer(trim(mode)//achar(0),lmode)

  file = f__fdopen(fd , lmode )
  deallocate(lmode )

end function fortran_fdopen

#if defined(SELF_TEST)
! s.f90 -DSELF_TEST fortran_stdio.F90
! ./a.out
program test_fortran_stdio
  use ISO_C_BINDING
  include 'fortran_stdio.inc'
  type(C_PTR) :: file
  integer, dimension(100), target :: buffer
  integer(C_SIZE_T) :: status, nelem, elemsize, nbytes
  integer :: i
  integer(C_LONG) :: where
  character(len=1), dimension(11), target :: text
  integer(C_INT) :: fd
  character(len=3) :: shortcircuit

  text = transfer('0123456789'//achar(10),text)

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

!  print *,'calling fortran_setbuffer to set file to unbuffered mode'
!  call fortran_setbuffer(file,C_NULL_PTR,0)
  print *,'calling fortran_fwrite'
  where = fortran_ftell(file)
  print *,'Where=',where
  status = fortran_fwrite(c_loc(buffer(1)) , elemsize , nelem , file)  ! write 1-> 100
  print *,'status=',status
  where = fortran_ftell(file)
  print *,'Where=',where

  print *,'calling fortran_fileno'
  fd = fortran_fileno(file)
  print *,'calling fortran_write on fd=',fd,' (maybe short-circuiting buffer)'
  nbytes = 40
  status = fortran_write( fd, c_loc(buffer(1)) , nbytes)  ! write 1-> 10 (maybe short-circuiting buffer)
  print *,'status=',status

  print *,'calling fortran_fflush'
  status = fortran_fflush(file)
  print *,'status=',status
  print *,'calling fortran_fileno'
  fd = fortran_fileno(file)
  print *,'calling fortran_write on fd=',fd
  status = fortran_write( fd, c_loc(buffer(1)) , nbytes)  ! write 1-> 10 at end
  print *,'status=',status
  

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

  print *,'expecting : 1->10, 1 -> 90 or 1->100'
  print *,buffer(1:status)
  shortcircuit = '   '
  if(buffer(status) == 100) shortcircuit = 'NOT'
  print *,'buffer was '//trim(shortcircuit)//' short-circuited'

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
  
  print *,'expecting : 41 -> 100, 1->10 or 51->100, 1->10, 1->10'
  print *,buffer(1:status)
  shortcircuit = '   '
  if(buffer(1) == 51) shortcircuit = 'NOT'
  print *,'buffer was '//trim(shortcircuit)//' short-circuited'

  print *,'calling fortran_fseek'
  where = 100
  status = fortran_fseek(file,where,SEEK_SET)
  print *,'status=',status
  where = fortran_ftell(file)
  print *,'Where=',where
  print *,'calling fortran_fileno'
  fd = fortran_fileno(file)
  print *,'calling fortran_read on fd=',fd
  nbytes = 80
  status = fortran_read(fd, c_loc(buffer(1)) , nbytes)
  
  print *,'expecting : 16 -> 35 or 26 -> 45'
  print *,buffer(1:status/4)
  shortcircuit = '   '
  if(buffer(1) == 26) shortcircuit = 'NOT'
  print *,'buffer was '//trim(shortcircuit)//' short-circuited'

  print *,'calling fortran_fclose'
  status = fortran_fclose(file)

  print *,'calling fortran_write on stderr'
  nbytes = 11
  status = fortran_write(2, C_LOC(text(1)), nbytes)
  print *,'status=',status

  print *,'calling fortran_write on stdout'
  status = fortran_write(1, C_LOC(text(1)), nbytes)
  print *,'status=',status

end program
#endif