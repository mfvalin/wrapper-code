! Copyright (C) 2021  Environnement et Changement climatique Canada
! 
! This is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
! 
! This software is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! Author:
!     M. Valin,   Environnement et Changement climatique Canada, 2020/2021
!
program test_plugin   ! test of fortran plugin module
  use ISO_C_BINDING
  use fortran_plugins
  implicit none

  type(plugin) :: sharedf1, sharedf2, sharedf3, shared2, anonymous
  logical :: status
  type(C_FUNPTR) :: fptr1, fptr2                  ! C pointer to function
  type(C_PTR) :: entry_ptr
  integer :: nsym, ilen, i, answer
  character(C_CHAR), dimension(:), pointer :: symbol_name
  character(len=128) :: longstr
  integer(C_INT64_T) :: ifptr

  interface
    subroutine c_exit(code) bind(C,name='exit')
      import :: C_INT
      integer(C_INT), intent(IN), value :: code
    end subroutine c_exit
  end interface

  procedure(procval) , pointer :: fpl2 => NULL()  ! pointer to integer function with one value integer argument
  abstract interface                              ! abstract interface needed for pointer to function
    integer function procval(arg) BIND(C)         ! with argument passed by value (BIND(C) is MANDATORY)
      integer, intent(IN), value :: arg
    end function procval
  end interface

  procedure(procadr) , pointer :: fpl1 => NULL()  ! pointer to generic integer function with one integer argument
  abstract interface                              ! abstract interface needed for pointer to function
    integer function procadr(arg) BIND(C)         ! with argument passed by address
      integer, intent(IN) :: arg
    end function procadr
  end interface

  call sharedf1 % diag(VERBOSE)                                         ! set verbose diagnostics
! call sharedf1 % diag(SILENT)                                          ! set non verbose diagnostics
  status = sharedf1 % load('libsharedf1.so')                            ! load sharedf1    (slot 0)
  status = sharedf1 % unload()                                          ! unload sharedf1  (slot 0)

  status = sharedf3 % load('libsharedf3.so')                            ! load sharedf2    (slot 0)
  print *,'load libsharedf3, status =',status
  status = sharedf2 % load('libsharedf2.so')                            ! load sharedf2    (slot 0)
  print *,'load libsharedf2, status =',status
  status = sharedf1 % load('libsharedf1.so')                            ! load sharedf1    (slot 1)
  print *,'load libsharedf1, status =',status
  status = sharedf3 % unload()                                          ! unload sharedf3
  print *,'unload libsharedf3, status =',status

  print *,"========================================"
  print *,'looking up name4f in libsharedf1.so'
  fptr1 = sharedf1 % fnptr('name4f')                                    ! known valid name in this plugin library
  if(.not. c_associated(fptr1)) then                                    ! OOPS !
    print *,'name4f not found (NOT EXPECTED)'
  else
    nsym = sharedf1 % symbols()                                         ! get number of advertized symbols
    print *,'number of advertized entries in plugin libsharedf1.so =', nsym
    do i = 1, nsym                                                      ! loop over number of symbols
      longstr = ""
      status = sharedf1 % fname(i, longstr)                             ! name of symbol at position i in list
      fptr2 = sharedf1 % fnptr(trim(longstr))
      ifptr = transfer(fptr2, ifptr)
     print 100,'entry ',i," = '",trim(longstr)//"'",ifptr
100   format(A,I3,A,A,2X,Z16.16)
    enddo
  endif
  if( c_associated(fptr1) ) then
    call c_f_procpointer(fptr1,fpl1)      ! make Fortran function pointer from C function pointer
    answer = fpl1(123)                    ! call by address (reference) type
    print *, 'fpl1(123) =', answer
    if(answer .ne. 123) print *,"ERROR: expecting 123, got",answer
  else
    print *, "ERROR: failed to find entry 'name4f'"
  endif
  fptr1 = sharedf1 % fnptr('unadvertised')  ! unadvertised name
  if( c_associated(fptr1) ) then
    call c_f_procpointer(fptr1,fpl2)      ! make Fortran function pointer from C function pointer
    answer = fpl2(789)                    ! call by value type
    print *, 'fpl2(789) =', answer
    if(answer .ne. 789) print *,"ERROR: expecting 789, got",answer
  else
    print *, "ERROR: failed to find entry 'unadvertised'"
  endif
  print *,"========================================"

  print *,'looking up name4f in libsharedf2.so'
  fptr1 = sharedf2 % fnptr('name4f')                                    ! known valid name in this plugin library
  if(.not. c_associated(fptr1)) then                                    ! OOPS !
    print *,'name4f not found (NOT EXPECTED)'
  else
    nsym = sharedf2 % symbols()                                         ! get number of advertized symbols
    print *,'number of advertized entries in plugin libsharedf2.so =', nsym
    do i = 1, nsym                                                      ! loop over number of symbols
      longstr = ""
      status = sharedf2 % fname(i, longstr)                             ! name of symbol at position i in list
      fptr2 = sharedf2 % fnptr(trim(longstr))
      ifptr = transfer(fptr2, ifptr)
     print 100,'entry ',i," = '",trim(longstr)//"'",ifptr
    enddo
  endif
  if( c_associated(fptr1) ) then
    call c_f_procpointer(fptr1,fpl1)      ! make Fortran function pointer from C function pointer
    answer = fpl1(345)                    ! call by address (reference) type
    print *, 'fpl1(345) =', answer
    if(answer .ne. 345) print *,"ERROR: expecting 345, got",answer
  else
    print *, "ERROR: failed to find entry 'name4f'"
  endif
  fptr1 = sharedf2 % fnptr('unadvertised')  ! unadvertised name
  if( c_associated(fptr1) ) then
    call c_f_procpointer(fptr1,fpl2)      ! make Fortran function pointer from C function pointer
    answer = fpl2(678)                    ! call by value type
    print *, 'fpl2(678) =', answer
    if(answer .ne. 678) print *,"ERROR: expecting 678, got",answer
  else
    print *, "ERROR: failed to find entry 'unadvertised'"
  endif
  print *,"========================================"

  print *,'looking up Name_2 in libsharedf1.so'
  fptr2 = sharedf1 % fnptr('Name_2')      ! bad name because wrong plugin library
  if(.not. c_associated(fptr2)) then
    print *,'Name_2 not found (AS EXPECTED)'
  else
    print *, 'ERROR: This should not print'
  endif
  print *,'looking up Name_2 in libsharedf2.so'
  fptr2 = sharedf1 % fnptr('Name_2')      ! bad name because wrong plugin library
  if(.not. c_associated(fptr2)) then
    print *,'Name_2 not found (AS EXPECTED)'
  else
    print *, 'ERROR: This should not print'
  endif
  print *,"========================================"

  status = shared2 % load('libshared2.so')           ! load shared2    (slot 0)
  
  print *,'looking up Name_2 in libshared2.so'
  fptr2 = shared2 % fnptr('Name_2')
  if(.not. c_associated(fptr2)) then
    print *,'Name_2 not found (NOT EXPECTED)'
  else
    nsym = shared2 % symbols()
    print *,'number of advertized entries in plugin libshared2.so =', nsym
    call c_f_procpointer(fptr2,fpl2)      ! make Fortran function pointer from C function pointer
    answer = fpl2(456)                    ! call by value type
    print *, 'fpl2(456) =',answer
    if(answer .ne. 456) print *,"ERROR: expecting 456, got",answer
  endif
  print *,"========================================"

  print *,'looking up Name_3 in ANY plugin'
  fptr2 = anonymous % fnptr('Name_3')     ! blind call, anonymous is not initialized
  if(.not. c_associated(fptr2)) then
    print *,'Name_3 not found (NOT EXPECTED)'
  else
    call c_f_procpointer(fptr2,fpl2)
    answer = fpl2(789)                    ! call by value type
    print *, 'fpl2(789) =',answer
    if(answer .ne. 789) print *,"ERROR: expecting 789, got",answer
  endif

  print *,'looking up unadvertised in ANY plugin'
  fptr2 = anonymous % fnptr('unadvertised')     ! blind call, anonymous is not initialized
  if(.not. c_associated(fptr2)) then
    print *,'unadvertised not found (NOT EXPECTED)'
  else
    call c_f_procpointer(fptr2,fpl2)
    answer = fpl2(567)                    ! call by value type
    print *, 'fpl2(567) =',answer
    if(answer .ne. 567) print *,"ERROR: expecting 567, got",answer
  endif
  print *,"========================================"

  status = sharedf1 % library(longstr)
  print *,"closing plugin library '"//trim(longstr)//"'"
  status = sharedf1 % unload()                ! unload sharedf1  (slot 1)
  print *,'unload '//trim(longstr)//', status =',status

  status = sharedf2 % library(longstr)
  print *,"closing plugin library '"//trim(longstr)//"'"
  status = sharedf2 % unload()                ! unload sharedf1  (slot 1)
  print *,'unload '//trim(longstr)//', status =',status

  status = shared2 % library(longstr)
  print *,"closing plugin library '"//trim(longstr)//"'"
  status = shared2 % unload()                 ! unload shared2  (slot 0)
  print *,'unload '//trim(longstr)//', status =',status

end program
