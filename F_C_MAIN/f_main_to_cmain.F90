program fcmain
  use ISO_C_BINDING
  implicit none
  integer(C_INT) :: nargs
  integer :: i, length, status
  character(len=4096) :: argument
  character(len=1), dimension(:), pointer :: arg1
  type(C_PTR), dimension(:), pointer :: argv
  type(C_PTR) :: argtab
  interface
    function c_main(nargs,argv) result(status) BIND(C,name='MY_C_MAIN')
    import
    implicit none
    integer, intent(IN), value :: nargs
    type(C_PTR), intent(IN), value :: argv
    integer :: status
    end function c_main
  end interface

  nargs = command_argument_count()
  allocate(argv(0:nargs+1))
  argv = C_NULL_PTR
  do i=0,nargs
    call get_command_argument(i,argument,length,status)
    allocate(arg1(length+1))
    arg1 = transfer(trim(argument)//achar(0),arg1,length+1)
    argv(i) = C_LOC(arg1(1))
  enddo
  argtab = C_LOC(argv(0))
  status = c_main(nargs,argtab)
  stop
end
