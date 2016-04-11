
FC = gfortran
CC = cc
FFLAGS = -O
CFLAGS = -O
#FFLAGS = -Wall -O
#CFLAGS = -Wall -O

all:	obj libs abs

libs:	libshared0.so libshared1.so libshared2.so libsharedf1.so libsharednone.so

obj: plugin_code.o

abs: test_c test_fortran

plugin_code.o:	plugin_code.c plugins.h
	$(CC) $(CFLAGS) -c -I. $<

libshared0.so:	shared0.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libshared1.so:	shared1.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libshared2.so:	shared2.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libsharednone.so:	sharednone.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libsharedf1.so:	sharedf1.F90
	$(FC) $(FFLAGS) -I. -shared -o $@ -fpic $<

test_c:	test_c.c plugin_code.o
	$(CC) $(CFLAGS) -I. -o $@ $< plugin_code.o -ldl -Wl,-rpath,.

test_fortran:	test_fortran.F90 plugin_code.o
	$(FC) $(FFLAGS) -I. -o $@ $< plugin_code.o -ldl -Wl,-rpath,.

clean:
	rm -f *.so *.o *.mod test_fortran test_c