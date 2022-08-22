! Example og plugin written in Fortran
! this plugin supplies (and advertises) the following callable entry points 
! name1f, name2f, name3f, name4f
! BIND(C,name=...) mandatory to make sure names match entry point names
!
! beginning of routines in plugin
integer(C_INT) function  fn1(arg) BIND(C,name='name1f')
use ISO_C_BINDING
implicit none
integer(C_INT), intent(IN) :: arg
print *,'in fortran function name1f from sharedf2 : arg =',arg
fn1 = arg
return
end

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
! to be adjusted by user : MAX_NAMES, MAX_NAME_LENGTH, calls to insert_in_name_table in subroutine symbols
module interop
  use ISO_C_BINDING
  implicit none
  save
! start of user adjusted code
#define MAX_NAMES 4
#define MAX_NAME_LENGTH 10
! end of user adjusted code
  type(C_PTR), dimension(MAX_NAMES+1), BIND(C,name='EntryList_') :: name_table ! mandatory external name
  character(len=1), dimension(MAX_NAME_LENGTH+1,MAX_NAMES), target :: names
  integer :: nargs = 0
  interface
    subroutine constructor() bind(C,name='FortranConstructor')
    end subroutine constructor
  end interface
contains
  subroutine insert_in_name_table(name)  ! add name to name table and name pointers
    use ISO_C_BINDING
    implicit none
    character(len=*) :: name
    integer :: l
    l = len(trim(name)) + 1
    nargs = min(nargs+1, MAX_NAMES)   ! make sure that list can be null terminated
    names(1:l,nargs) = transfer(trim(name)//achar(0) , names, l)
    name_table(nargs) = C_LOC(names(1,nargs))
    name_table(nargs+1) = C_NULL_PTR   ! make sure list is null terminated
    return
  end subroutine insert_in_name_table
  function symbols() bind(C,name='get_symbol_number') result(number)  ! this function is mandatory in a valid Fortran plugin
    use ISO_C_BINDING
    implicit none
    integer(C_INT) :: number
    if(nargs == 0) then
      print *,"calling insert_symbols from symbols(sharedf2)"
      call insert_symbols()
    endif
    number = nargs   ! return number of arguments
    return
  end function symbols
  subroutine insert_symbols() bind(C,name='fortran_constructor')
    use ISO_C_BINDING
    implicit none
    nargs = 0
! start of user adjusted code
print *,'automatic insertion of symbols for sharedf2'
    call insert_in_name_table('name1f')
    call insert_in_name_table('name2f')
    call insert_in_name_table('name3f')
    call insert_in_name_table('name4f')
! end of user adjusted code
    if(nargs == -1) call constructor  ! never called, but forces FortranConstructor to be loaded from library
    return
  end subroutine insert_symbols
end module interop
