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
!> \file
!> \brief interface  to C plugin library management routines (object oriented)
module fortran_plugins
  use ISO_C_BINDING
  implicit none
  include 'plugins.inc'

  !> \private
  interface
    function c_strlen_(string) result(length) BIND(C,name='strlen')
      import C_SIZE_T, C_PTR
      type(C_PTR), value :: string
      integer(C_SIZE_T) :: length
    end function c_strlen_
  end interface

  integer, parameter :: SILENT  = 0               !< constant for verbosity control
  integer, parameter :: VERBOSE = 1               !< constant for verbosity control

  !> \brief plugin library user defined type
  type :: plugin
    !> \private
    private
    type(C_PTR) :: h = C_NULL_PTR                               !< plugin handle
    character(len=:), allocatable :: libname                    !< plugin library name (deferred-length-string)
  contains
    !> \return true if successful, false otherwise
    procedure         :: load    => load_plugin_library         !< load plugin library
    !> \return true if successful, false otherwise
    procedure         :: unload  => unload_plugin_library       !< unload plugin library
    !> \return true if valid, false otherwise
    procedure         :: valid   => valid_plugin_library        !< check validity (open plugin library)
    !> \return true if valid library, false otherwise
    procedure         :: library => plugin_library_name         !< name of this open plugin library
    !> \return none
    procedure, NOPASS :: diag    => plugin_library_diag         !< verbose / silent mode
    !> \return function address, C_NULL_PTR if not found
    procedure         :: fnptr   => plugin_function_pointer     !< get function address
    !> \return true if valid library, false otherwise
    procedure         :: fname   => plugin_function_symbol      !< get function name
    !> \return number of symbols found in library, 0 if nothing found
    procedure         :: symbols => plugin_function_symbols     !< get number of advertized functions
  end type

  contains

  function load_plugin_library(this, libname) result(status)  !< load (open) plugin library called libname
    implicit none
    class(plugin), intent(INOUT) :: this                       !< plugin type
    character(len=*), intent(IN) :: libname                    !< file name of library to open
    logical :: status                                          !< true if successful, false otherwise
    integer :: lng
    status = .false.
    if(trim(libname) .ne. "") then
      this % h = load_plugin(trim(libname)//C_NULL_CHAR)
      this % libname = trim(libname)
    else
      this % h = C_NULL_PTR
    endif
    status = C_ASSOCIATED(this %h)
  end function load_plugin_library

  function plugin_library_name(this, libname) result(status)   !< get name of (open) plugin library
    implicit none
    class(plugin), intent(INOUT)  :: this                      !< plugin type
    character(len=*), intent(OUT) :: libname                   !< file name of library
    logical :: status                                          !< true if successful, false otherwise
    status = .false.
    if(C_ASSOCIATED(this %h) .and. trim(this % libname) .ne. "") then
      libname = trim(this % libname)
      status = .true.
    else
      libname = ""
    endif
  end function plugin_library_name

  function unload_plugin_library(this) result(status)          !< unload (close) this plugin library
    implicit none
    class(plugin), intent(INOUT) :: this                       !< plugin type
    logical :: status                                          !< true if successful, false otherwise
    integer :: ival
    ival = unload_plugin(this % h)
    this % h = C_NULL_PTR
    this % libname = ""
!     if(associated(this % lib)) deallocate(this % lib)
    status = ival == 0
  end function unload_plugin_library

  function valid_plugin_library(this) result(status)           !< is this associated to an open shared library
    implicit none
    class(plugin), intent(INOUT) :: this                       !< plugin type
    logical :: status                                          !< true if successful, false otherwise
    status = C_ASSOCIATED(this % h)
  end function valid_plugin_library

  subroutine plugin_library_diag(diag)                         !< set verbose / silent mode
    implicit none
    integer, intent(IN), value :: diag                         !< verbosity level SILENT / VERBOSE
    call set_plugin_diag(diag)
  end subroutine plugin_library_diag

  function plugin_function_symbols(this) result(nsym)          !< how many functions are advertized in this plub=gin library
    implicit none
    class(plugin), intent(IN) :: this                          !< plugin type
    integer :: nsym                                            !< number of advertized functions
    if(C_ASSOCIATED(this % h)) then
      nsym = plugin_n_functions(this % h)
    else
      nsym = 0
    endif
  end function plugin_function_symbols

  function plugin_function_pointer(this, function_name) result(fptr)   !< get pointer to function called function_name in this plugin
    class(plugin), intent(IN) :: this                          !< plugin type
    character(len=*), intent(IN) :: function_name              !< function name
    type(C_FUNPTR) :: fptr                                     !< pointer to function (C_NULL_PTR if not found)
    if(C_ASSOCIATED(this % h) .and. trim(function_name) .ne. "") then
      fptr = plugin_function(this % h, trim(function_name)//C_NULL_CHAR )   ! lookup in this plugin library
    else
      fptr = plugin_function(this % h, trim(function_name)//C_NULL_CHAR )   ! anonymous (all known plugins) lookup if NULL pointer
    endif
  end function plugin_function_pointer

  function plugin_function_symbol(this, ordinal, symbol_name) result(status)  !< get name of symbol #ordinal in this plugin
    implicit none
    class(plugin), intent(IN) :: this                          !< plugin type
    integer, intent(IN), value :: ordinal                      !< ordinal in advertized function table
    character(len=*), intent(INOUT) :: symbol_name             !< symbol (function) name
    logical :: status
    type(C_PTR) :: str
    integer :: lng
    character(C_CHAR), dimension(:), pointer :: temp
    symbol_name = ""
    if(C_ASSOCIATED(this % h) .and. ordinal >= 0) then
      str = plugin_function_name(this % h, ordinal)
      lng = int( c_strlen_(str) )
      lng = min(lng,len(symbol_name))
      call c_f_pointer(str,temp,[lng])
      symbol_name(1:lng) = transfer(temp(1:lng),symbol_name)
      status = .true.
    else
      status = .false.
    endif
  end function plugin_function_symbol

end module
