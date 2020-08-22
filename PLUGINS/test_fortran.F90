program test_plugin
  use ISO_C_BINDING
  implicit none
  include 'plugins.inc'
  interface
    function c_strlen(string) result(length) BIND(C,name='strlen')
      import C_SIZE_T, C_PTR
      type(C_PTR), value :: string
      integer(C_SIZE_T) :: length
    end function c_strlen
  end interface
  type(C_FUNPTR) :: fptr1, fptr2                  ! C pointer to function

  procedure(procadr) , pointer :: fpl1 => NULL()  ! pointer to generic integer function with one integer argument

  procedure(procval) , pointer :: fpl2 => NULL()  ! pointer to integer function with one value integer argument

  abstract interface                              ! abstract interface needed fot pointer to function
    integer function procval(arg)
      integer, intent(IN), value :: arg
    end function procval
  end interface
  abstract interface                              ! abstract interface needed for pointer to function
    integer function procadr(arg)
      integer, intent(IN) :: arg
    end function procadr
  end interface

  type(C_PTR) :: handle1, handle2
  character(C_CHAR), dimension(:), pointer :: name
  integer :: ival, ilen
  type(C_PTR) :: entry_ptr
  character(len=128) :: longstr

  call set_plugin_diag(1)
! handle1 = load_plugin( transfer('libsharedf1.so'//C_NULL_CHAR,name) )  ! with transfer
  handle1 = load_plugin(C_CHAR_'libsharedf1.so'//C_NULL_CHAR)           ! load handle1 (with constant retyping)
  ival = unload_plugin(handle1);                                             ! unload handle1  (slot 0)

  handle2 = load_plugin( C_CHAR_'libshared2.so'//C_NULL_CHAR)           ! load handle2    (slot 0)
  handle1 = load_plugin(C_CHAR_'libsharedf1.so'//C_NULL_CHAR)           ! load handle1    (slot 1)
  ival = unload_plugin(handle2);                                             ! unload handle2  (slot 0)

  print *,'looking up name4f in libsharedf1.so'
  fptr1 = plugin_function(handle1, transfer('name4f'//C_NULL_CHAR,name) )  ! good name
  if(.not. c_associated(fptr1)) then
    print *,'name4f not found (NOT EXPECTED)'
  else
    print *,'number of entries in plugin libsharedf1.so =',plugin_n_functions(handle1)
    do ival = 1, plugin_n_functions(handle1)
      entry_ptr = plugin_function_name(handle1,ival)
      call c_f_pointer(entry_ptr,name,[128])           ! make Fortran pointer from C pointer
      ilen = int( c_strlen(entry_ptr) )
      longstr = ""
      longstr(1:ilen) = transfer(name(1:ilen),longstr)
     print 100,'entry ',ival,' = ',transfer(plugin_function(handle1,name),1_8),trim(longstr)
100   format(A,I3,A,Z16.16,2X,A)
    enddo
    call c_f_procpointer(fptr1,fpl1)      ! make Fortran function pointer from C function pointer
   ival = fpl1(123)
    print *, 'ival =',ival
  endif
  fptr1 = plugin_function(handle1, transfer('unadvertised'//C_NULL_CHAR,name) )  ! unadvertised name
  if( c_associated(fptr1) ) then
    call c_f_procpointer(fptr1,fpl2)
    ival = fpl2(789)
    print *, 'ival =',ival
  else
  endif

  print *,'looking up Name_2 in libsharedf1.so'
  fptr2 = plugin_function(handle1, C_CHAR_'Name_2'//C_NULL_CHAR )  ! bad name because wrong plugin
  if(.not. c_associated(fptr2)) then
    print *,'Name_2 not found (expected)'
  else
    print *, 'This should not print'
  endif

  handle2 = load_plugin( C_CHAR_'libshared2.so'//C_NULL_CHAR )           ! load handle2    (slot 0)
  
  print *,'looking up Name_2 in libshared2.so'
  fptr2 = plugin_function(handle2, C_CHAR_'Name_2'//C_NULL_CHAR )
  if(.not. c_associated(fptr2)) then
    print *,'Name_2 not found (NOT EXPECTED)'
  else
    print *,'number of entries in plugin libshared2.so =',plugin_n_functions(handle2)
    call c_f_procpointer(fptr2,fpl2)      ! make Fortran function pointer from C function pointer
    ival = fpl2(456)
    print *, 'ival =',ival
  endif
  
  print *,'looking up Name_3 in ANY plugin'
  fptr2 = plugin_function(C_NULL_PTR, C_CHAR_'Name_3'//C_NULL_CHAR )  ! blind call
  if(.not. c_associated(fptr2)) then
    print *,'Name_3 not found (NOT EXPECTED)'
  else
    call c_f_procpointer(fptr2,fpl2)
    ival = fpl2(789)
    print *, 'ival =',ival
  endif

  ival = unload_plugin(handle1);   ! unload handle1  (slot 1)
  ival = unload_plugin(handle2);   ! unload handle2  (slot 0)
end