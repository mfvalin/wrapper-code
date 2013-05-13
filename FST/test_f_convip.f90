program test_ip_conversion
! s.cc -c test_c_convip.c
! rm -f a.out 
! s.f90 test_f_convip.f90 conv_kind_15.f90  test_c_convip.o convert_ip.o convert_ip123.f90  -lrmn_013 
! ./a.out
interface cmain
 subroutine cmain() BIND(C,name='Cmain')
 end subroutine cmain
end interface

call cmain
stop
end
subroutine qqexit  ! from rmnlib
return
end
subroutine cigaxg  ! from rmnlib
stop
end
