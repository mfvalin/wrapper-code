! Example og plugin written in Fortran
! this plugin supplies (and advertises) the following callable entry points 
! name1f, name2f, name3f, name4f
! BIND(C,name=...) mandatory to make sure names match entry point names
!
! beginning of routines in plugin
! integer(C_INT) function  fn1(arg) BIND(C,name='name1f')
! use ISO_C_BINDING
! implicit none
! integer(C_INT), intent(IN) :: arg
! print *,'in fortran function name1f from sharedf2 : arg =',arg
! fn1 = arg
! return
! end

integer(C_INT) function  fn2(arg) BIND(C,name='name2f')
use ISO_C_BINDING
implicit none
integer(C_INT), intent(IN) :: arg
print *,'in fortran function name2f from sharedf2 : arg =',arg
fn2 = arg
return
end

integer(C_INT) function  fn3(arg) BIND(C,name='name3f')
use ISO_C_BINDING
implicit none
integer(C_INT), intent(IN) :: arg
print *,'in fortran function name3f from sharedf2 : arg =',arg
fn3 = arg
return
end

integer(C_INT) function  fn4(arg) BIND(C,name='name4f')
use ISO_C_BINDING
implicit none
integer(C_INT), intent(IN) :: arg
print *,'in fortran function name4f from sharedf2 : arg =',arg
fn4 = arg
return
end

integer(C_INT) function  fn5(arg) BIND(C,name='unadvertised')  ! unadvertised entry, call by value
use ISO_C_BINDING
implicit none
integer(C_INT), intent(IN), value :: arg
print *,'in fortran function unadvertised from sharedf2 : arg =',arg
fn5 = arg
return
end
! end of routines in plugin
!
! what follows is boiler plate code
! to be adjusted by user : MAX_NAMES, MAX_NAME_LENGTH, 
!                          calls to insert_in_name_table in subroutine user_symbols
#define MAX_NAMES 4
#define MAX_NAME_LENGTH 10
#include <library_plugin_mod.hf>
subroutine fortran_constructor() bind(C,name='fortran_constructor')
  use library_plugin_mod
  implicit none
print *,'loading symbols in table for  sharedf2'
!   call insert_in_name_table('name1f')
  call insert_in_name_table('name2f')
  call insert_in_name_table('name3f')
  call insert_in_name_table('name4f')
! perform any tasks needed for library initialization here
  return
end
subroutine fortran_destructor() bind(C,name='fortran_destructor')
print *,'in subroutine fortran_destructor for sharedf1'
end 
