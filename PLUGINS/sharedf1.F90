! Example og plugin written in Fortran
! this plugin supplies (and advertises) the following callable entry points 
! name1f, name2f, name3f, name4f
! BIND(C,name=...) mandatory to make sure names match etry pooint names
!
! beginning of routines in plugin
integer function  fn1(arg) BIND(C,name='name1f')
integer, intent(IN) :: arg
print *,'fortran name1 =',arg
fn1 = arg
return
end

integer function  fn2(arg) BIND(C,name='name2f')
integer, intent(IN) :: arg
print *,'fortran name2 =',arg
fn2 = arg
return
end

integer function  fn3(arg) BIND(C,name='name3f')
integer, intent(IN) :: arg
print *,'fortran name3f =',arg
fn3 = arg
return
end

integer function  fn4(arg) BIND(C,name='name4f')
integer, intent(IN) :: arg
print *,'fortran name4f =',arg
fn4 = arg
return
end

integer function  fn5(arg) BIND(C,name='unadvertised')  ! unadvertised entry, call by value
integer, intent(IN), value :: arg
print *,'fortran unadvertised =',arg
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
! start of user adjusted code
#define MAX_NAMES 4
#define MAX_NAME_LENGTH 10
! end of user adjusted code
  type(C_PTR), dimension(MAX_NAMES+1), save, target, BIND(C,name='entry_list') :: name_table
  character(len=1), dimension(MAX_NAME_LENGTH+1,MAX_NAMES), save, target :: names
  integer, save :: nargs
  contains
  subroutine insert_in_name_table(name)  ! add name to name table and neme pointers
    use ISO_C_BINDING
    implicit none
    character(len=*) :: name
    integer :: l
    l = len(trim(name)) + 1
    nargs = nargs + 1
    names(1:l,nargs) = transfer(trim(name)//achar(0) , names, l)
    name_table(nargs) = C_LOC(names(1,nargs))
    return
  end subroutine insert_in_name_table
  function symbols() bind(C,name='get_symbol_number') result(number)
    use ISO_C_BINDING
    implicit none
    integer(C_INT) :: number
    nargs = 0
! start of user adjusted code
    call insert_in_name_table('name1f')
    call insert_in_name_table('name2f')
    call insert_in_name_table('name3f')
    call insert_in_name_table('name4f')
! end of user adjusted code
    number = nargs   ! return number of arguments
    return
  end function symbols
end module interop
