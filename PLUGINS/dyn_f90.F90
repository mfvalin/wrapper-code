function by_value(arg) result(val) BIND(C,name='by_val')
  use ISO_C_BINDING
  implicit none
  integer(C_INT), value, intent(IN) :: arg
  integer(C_INT) :: val
  print *,'by_val F90 arg=',arg
  val = arg
end function by_value

function by_address(arg) result(val) BIND(C,name='by_adr')
  use ISO_C_BINDING
  implicit none
  integer(C_INT),intent(IN) :: arg
  integer(C_INT) :: val
  print *,'by_adr F90 arg=',arg
  val = arg
end function by_address
