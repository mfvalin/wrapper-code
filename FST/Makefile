
default:
	@echo dummy target default

test:
	rm -f a.out *.o *.mod
#	s.cc -c -I. *.c 
#	s.f77 -c convert_ip.f
#	s.f90 -o a.out test_f_convip.f90  test_c_convip.o convert_ip123.f90
	s.compile -src test_c_convip.c convert_ip123.f90 -bidon pure_fortran -main convip_all_tests -o a.out
	./a.out
	rm -f a.out *.o *.mod

clean:
	rm -f *.o *.out *.mod
