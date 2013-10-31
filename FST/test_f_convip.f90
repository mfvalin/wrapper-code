program test_f_convip
! s.cc -c test_c_convip.c
! rm -f a.out 
! s.f90 test_c_convip.o  convert_ip123.f90  test_f_convip.f90
! ./a.out
call convip_all_tests
stop
end
!===============================================================================================
subroutine convip_all_tests
  use convert_ip123, only : test_value_to_string, test_convip_plus, value_to_string
  implicit none
  character(len=15) :: string
  integer :: fmtcode
  print *,'========== convip_unit_tests =========='
  call convip_unit_tests
  print *,'========== test_value_to_string =========='
  call test_value_to_string
  fmtcode = value_to_string(1.0,string,len(string))
  print *,'format code for 1.0 is',fmtcode,' (expected 2)'
  print *,'========== test_convip_plus =========='
  call test_convip_plus
end subroutine convip_all_tests
!===============================================================================================
subroutine convip_unit_tests

interface main1
 subroutine main1() BIND(C,name='Cmain1')
 end subroutine main1
end interface
interface main2
 subroutine main2() BIND(C,name='Cmain2')
 end subroutine main2
end interface

call main1
call main2
return
end subroutine convip_unit_tests
