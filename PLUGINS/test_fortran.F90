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

  type(plugin) :: sharedf1, sharedf2, shared2, anonymous
  logical :: status
  type(C_FUNPTR) :: fptr1, fptr2                  ! C pointer to function
  type(C_PTR) :: entry_ptr
  integer :: nsym, ilen, i, answer
  character(C_CHAR), dimension(:), pointer :: symbol_name
  character(len=128) :: longstr

  procedure(procval) , pointer :: fpl2 => NULL()  ! pointer to integer function with one value integer argument
  abstract interface                              ! abstract interface needed for pointer to function
    integer function procval(arg)                 ! with argument passed by value
      integer, intent(IN), value :: arg
    end function procval
  end interface

  procedure(procadr) , pointer :: fpl1 => NULL()  ! pointer to generic integer function with one integer argument
  abstract interface                              ! abstract interface needed for pointer to function
    integer function procadr(arg)                 ! with argument passed by address
      integer, intent(IN) :: arg
    end function procadr
  end interface

  call sharedf1 % diag(VERBOSE)                                         ! set verbose diagnostics
  status = sharedf1 % load('libsharedf1.so')                            ! load sharedf1    (slot 0)
  status = sharedf1 % unload()                                          ! unload sharedf1  (slot 0)

  status = sharedf2 % load('libsharedf2.so')                            ! load sharedf2    (slot 0)
  status = sharedf1 % load('libsharedf1.so')                            ! load sharedf1    (slot 1)
  status = sharedf2 % unload()                                          ! unload sharedf2  (slot 0)

  print *,'looking up name4f in libsharedf1.so'
  fptr1 = sharedf1 % fnptr('name4f')                                    ! known valid name in this plugin library
  if(.not. c_associated(fptr1)) then                                    ! OOPS !
    print *,'name4f not found (NOT EXPECTED)'
  else
    nsym = sharedf1 % symbols()                                         ! get number of advertized symbols
    print *,'number of entries in plugin libsharedf1.so =', nsym
    do i = 1, nsym                                                      ! loop over number of symbols
      longstr = ""
      status = sharedf1 % fname(i, longstr)                             ! name of symbol at position i in list
     print 100,'entry ',i," = '",trim(longstr)//"'"
100   format(A,I3,A,A)
    enddo
  endif

  if( c_associated(fptr1) ) then
    call c_f_procpointer(fptr1,fpl1)      ! make Fortran function pointer from C function pointer
    answer = fpl1(123)                    ! call by address (reference) type
    print *, 'fpl1(123) =', answer
  else
    print *, "ERROR: failed to find entry 'name4f'"
  endif

  fptr1 = sharedf1 % fnptr('unadvertised')  ! unadvertised name
  if( c_associated(fptr1) ) then
    call c_f_procpointer(fptr1,fpl2)      ! make Fortran function pointer from C function pointer
    answer = fpl2(789)                    ! call by value type
    print *, 'fpl2(789) =', answer
  else
    print *, "ERROR: failed to find entry 'unadvertised'"
  endif

  print *,'looking up Name_2 in libsharedf1.so'
  fptr2 = sharedf1 % fnptr('Name_2')      ! bad name because wrong plugin library
  if(.not. c_associated(fptr2)) then
    print *,'Name_2 not found (AS EXPECTED)'
  else
    print *, 'ERROR: This should not print'
  endif

  status = shared2 % load('libshared2.so')           ! load shared2    (slot 0)
  
  print *,'looking up Name_2 in libshared2.so'
  fptr2 = shared2 % fnptr('Name_2')
  if(.not. c_associated(fptr2)) then
    print *,'Name_2 not found (NOT EXPECTED)'
  else
    nsym = shared2 % symbols()
    print *,'number of entries in plugin libshared2.so =', nsym
    call c_f_procpointer(fptr2,fpl2)      ! make Fortran function pointer from C function pointer
    answer = fpl2(456)                    ! call by value type
    print *, 'fpl2(456) =',answer
  endif
  
  print *,'looking up Name_3 in ANY plugin'
  fptr2 = anonymous % fnptr('Name_3')     ! blind call, anonymous is not initialized
  if(.not. c_associated(fptr2)) then
    print *,'Name_3 not found (NOT EXPECTED)'
  else
    call c_f_procpointer(fptr2,fpl2)
    answer = fpl2(789)                    ! call by value type
    print *, 'fpl2(789) =',answer
  endif

  status = sharedf1 % library(longstr)
  print *,"closing plugin library '"//trim(longstr)//"'"
  status = sharedf1 % unload()                ! unload sharedf1  (slot 1)

  status = shared2 % library(longstr)
  print *,"closing plugin library '"//trim(longstr)//"'"
  status = shared2 % unload()                 ! unload shared2  (slot 0)

end program
