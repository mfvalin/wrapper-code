!
! fortran interface to libc's fopen/fclose/fseek/ftell/fread/fwrite/read/write
!
! libc function xxx is called fortran_xxx
! the calling sequence is mimicked from C calling sequence
! (man fopen/fclose/fseek/ftell/fread/fwrite/read/write provides the documentation)
!
! char * is character(len=*)
! int    is TYPE(C_INT)
! long   is TYPE(C_LONG)
! size_t is TYPE(C_SIZE_T)
! void * is type(C_PTR)
! FILE * is type(C_PTR)
!
! SEEK constants are defined
!
integer(C_INT), parameter :: SEEK_SET = 0
integer(C_INT), parameter :: SEEK_CUR = 1
integer(C_INT), parameter :: SEEK_END = 2

interface

! FILE *fopen(const char *path, const char *mode)
  function fortran_fopen(path,mode) result(file)
    import
    implicit none
    character(len=*), intent(IN) :: path,mode
    type(C_PTR) :: file
  end function fortran_fopen

! FILE *fdopen(int fd, const char *mode)
  function fortran_fdopen(fd,mode) result(file)
    import
    implicit none
    character(len=*), intent(IN) :: mode
    integer(C_INT), intent(IN) :: fd
    type(C_PTR) :: file
  end function fortran_fdopen

! int fclose(FILE *stream)
  function fortran_fclose(file) result(status) bind(C,name='fclose')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
    integer(C_INT) :: status
  end function fortran_fclose

! size_t fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream)
  function fortran_fwrite(what, itemsize, n, file) result(status) bind(C,name='fwrite')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: what
    integer(C_SIZE_T), value, intent(IN) :: itemsize, n
    type(C_PTR), value, intent(IN) :: file
    integer(C_SIZE_T) :: status
  end function fortran_fwrite

! size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
  function fortran_fread(what,itemsize,n,file) result(status) bind(C,name='fread')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: what
    integer(C_SIZE_T), value, intent(IN) :: itemsize,n
    type(C_PTR), value, intent(IN) :: file
    integer(C_SIZE_T) :: status
  end function fortran_fread

! int fflush(FILE *stream)
  function fortran_fflush(file) result(status) bind(C,name='fflush')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
    integer(C_INT) :: status
  end function fortran_fflush

! int fseek(FILE *stream, long offset, int whence)
  function fortran_fseek(file,offset,whence) result(status) bind(C,name='fseek')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
    integer(C_LONG), value, intent(IN) :: offset
    integer(C_INT), value, intent(IN) :: whence
    integer(C_INT) :: status
  end function fortran_fseek

! long ftell(FILE *stream, long offset, int whence)
  function fortran_ftell(file) result(where) bind(C,name='ftell')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
    integer(C_LONG) :: where
  end function fortran_ftell

! void setbuffer(FILE *stream, char *buf, size_t size)
  subroutine fortran_setbuffer(file,buf,size) bind(C,name='setbuffer')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
    type(C_PTR), value, intent(IN) :: buf
    integer(C_SIZE_T), value, intent(IN) :: size
  end subroutine fortran_setbuffer

! int fileno(FILE *stream)
  function fortran_fileno(file) result(fd) bind(C,name='fileno')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
    integer(C_INT) :: fd
  end function fortran_fileno
  
! int feof(FILE *stream)
  function fortran_feof(file) result(status) bind(C,name='feof')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
    integer(C_INT) :: status
  end function fortran_feof

! int ferror(FILE *stream)
  function fortran_ferror(file) result(status) bind(C,name='ferror')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
    integer(C_INT) :: status
  end function fortran_ferror

! void clearerr(FILE *stream)
  subroutine fortran_clearerr(file) bind(C,name='clearerr')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: file
  end subroutine fortran_clearerr

! ssize_t read(int fd, void *buf, size_t count)
  function fortran_read(fd,what,n) result(status) bind(C,name='read')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: what
    integer(C_SIZE_T), value, intent(IN) :: n
    integer(C_INT), value, intent(IN) :: fd
    integer(C_SIZE_T) :: status
  end function fortran_read

! ssize_t write(int fd, const void *buf, size_t count)
  function fortran_write(fd,what,n) result(status) bind(C,name='write')
    import
    implicit none
    type(C_PTR), value, intent(IN) :: what
    integer(C_SIZE_T), value, intent(IN) :: n
    integer(C_INT), value, intent(IN) :: fd
    integer(C_SIZE_T) :: status
  end function fortran_write

! int close(int fd)
  function fortran_close(fd) result(status) bind(C,name='close')
    import
    implicit none
    integer(C_INT), value, intent(IN) :: fd
    integer(C_INT) :: status
  end function fortran_close

end interface
