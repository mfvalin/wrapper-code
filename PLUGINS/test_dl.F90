program test_dl
  use ISO_C_BINDING
  implicit none
  interface
    function c_dlopen(name) result(handle) BIND(C,name='cdlopen')
      import C_INT, C_CHAR, C_PTR
      character(C_CHAR), dimension(*), intent(IN) :: name
      type(C_PTR) :: handle
    end function c_dlopen
    function c_dlsym(handle,name) result(fptr) BIND(C,name='dlsym')
      import C_CHAR, C_PTR, c_FUNPTR
      type(C_PTR), value, intent(IN) :: handle
      character(C_CHAR), dimension(*), intent(IN) :: name
      type(C_FUNPTR) :: fptr
    end function c_dlsym
  end interface
  abstract interface
    function procval(arg) result(val)
      import C_INT
      integer(C_INT), value, intent(IN) :: arg
      integer(C_INT) :: val
    end function procval
  end interface

  type(C_PTR) :: handle1, handle2
  type(C_FUNPTR) :: fptr1, fptr2
  procedure(integer(C_INT)), pointer :: fptradr => NULL()
  procedure(procval), pointer :: fptrval => NULL()
  integer(C_INT) :: temp

  handle1 = c_dlopen(C_CHAR_'libdyn_c.so'//C_NULL_CHAR)
  if(.not. c_associated(handle1)) then
    print *,'ERROR: handle1 is NULL'
    stop
  endif
  fptr1 = c_dlsym(handle1,C_CHAR_'by_adr'//C_NULL_CHAR)
  fptr2 = c_dlsym(handle1,C_CHAR_'by_val'//C_NULL_CHAR)
  if(.not. ( c_associated(fptr1) .and.  c_associated(fptr2) ) ) then
    print *,'ERROR: fptr1 or ftpr2 is NULL'
    stop
  endif
  call C_F_PROCPOINTER(fptr1,fptradr)
  temp = fptradr(123)
  call C_F_PROCPOINTER(fptr2,fptrval)
  temp = fptrval(123)

  handle2 = c_dlopen(C_CHAR_'libdyn_f90.so'//C_NULL_CHAR)
  if(.not. c_associated(handle2)) then
    print *,'ERROR: handle2 is NULL'
    stop
  endif
  fptr1 = c_dlsym(handle2,C_CHAR_'by_adr'//C_NULL_CHAR)
  fptr2 = c_dlsym(handle2,C_CHAR_'by_val'//C_NULL_CHAR)
  if(.not. ( c_associated(fptr1) .and.  c_associated(fptr2) ) ) then
    print *,'ERROR: fptr1 or ftpr2 is NULL'
    stop
  endif
  call C_F_PROCPOINTER(fptr1,fptradr)
  temp = fptradr(123)
  call C_F_PROCPOINTER(fptr2,fptrval)
  temp = fptrval(123)
end program