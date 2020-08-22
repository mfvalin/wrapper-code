program test_dl
  use ISO_C_BINDING
  implicit none
  
  integer, parameter :: RTLD_LAZY = 1
  integer, parameter :: RTLD_NOW = 2
  integer, parameter :: RTLD_GLOBAL = 256
  integer, parameter :: RTLD_LOCAL = 0
  interface
  function c_dlopen(path, flags) result(handle) bind(C,name='dlopen')
    import :: C_INT, C_PTR, C_CHAR
    integer(C_INT), intent(IN),    value :: flags
    character(C_CHAR), dimension(*), intent(IN) :: path
    type(C_PTR) :: handle
  end function c_dlopen
  function c_dlsym(handle, name) result(address) bind(C,name='dlsym')
    import :: C_PTR, C_CHAR, C_FUNPTR
    type(C_PTR), intent(IN),    value :: handle
    character(C_CHAR), dimension(*), intent(IN) :: name
    type(C_FUNPTR) :: address
  end function c_dlsym
  function c_dlclose(handle) result(ok) bind(C,name='dlclose')
    import :: C_INT, C_PTR
    type(C_PTR), intent(IN),    value :: handle
    integer(C_INT) :: ok
  end function c_dlclose
  function c_dlerror() result(cstring) bind(C,name='dlerror')
    import :: C_PTR
    type(C_PTR) :: cstring
  end function c_dlerror
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

  handle1 = c_dlopen(C_CHAR_'libdyn_c.so'//C_NULL_CHAR, RTLD_LAZY)
  if(.not. c_associated(handle1)) then
    print *,'ERROR: handle1 is NULL, libdyn_c.so not found'
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

  handle2 = c_dlopen(C_CHAR_'libdyn_f90.so'//C_NULL_CHAR, RTLD_LAZY)
  if(.not. c_associated(handle2)) then
    print *,'ERROR: handle2 is NULL, libdyn_f90.so not found'
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