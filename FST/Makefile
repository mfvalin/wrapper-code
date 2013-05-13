
default:
	@echo dummy target default

test:
	rm -f a.out *.o *.mod
	s.cc -c -I. *.c 
	s.f77 -c convert_ip.f
	s.f90 -o a.out test_f_convip.f90 conv_kind_15.f90  test_c_convip.o convert_ip.o convert_ip123.f90
	./a.out

clean:
	rm -f *.o *.out *.mod
