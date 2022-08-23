
FC = gfortran
CC = cc
FFLAGS = -O -fpic
CFLAGS = -O -fpic
#FFLAGS = -Wall -O
#CFLAGS = -Wall -O

all:	inc obj libs abs

libs:	libplugin.a libbuildplugin.a libshared0.so libshared1.so libshared2.so \
	libsharedf1.so libsharedf2.so libsharednone.so libdyn_c.so libdyn_f90.so

obj: plugin_code.o

inc:
	grep InTf *.c | sed -e 's:^.*//::' -e 's/[ ]*!InTf!.*//' >plugins.inc
	grep InTc plugin_code.c | sed 's:[ ]*//InTc.*:;:'        >plugins.h

abs: test_c test_dl_f test_fortran

plugin_code.o:	plugin_code.c plugins.h
	$(CC) $(CFLAGS) -c -I. $<

libplugin.a : plugin_code.o
	ar rcv $@ $^

libbuildplugin.a : FortranConstructor.o
	ar rcv $@ FortranConstructor.o

libshared0.so:	shared0.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libdyn_c.so:	dyn_c.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libdyn_f90.so: dyn_f90.F90
	$(FC) $(FFLAGS)  -I. -shared -o $@ -fpic $<

libshared1.so:	shared1.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libshared2.so:	shared2.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libsharednone.so:	sharednone.c
	$(CC) $(CFLAGS) -I. -shared -o $@ -fpic $<

libsharedf1.so:	sharedf1.F90 libplugin.a
	$(FC) $(FFLAGS) -I. -shared -o $@ -fpic $< # -L. -lbuildplugin

libsharedf2.so:	sharedf2.F90 libplugin.a
	$(FC) $(FFLAGS) -I. -shared -o $@ -fpic $< # -L. -lbuildplugin

test_c:	test_c.c plugin_code.o libs
	$(CC) $(CFLAGS) -I. -o $@ $< plugin_code.o -ldl -Wl,-rpath,.

test_dl_f: test_dl.F90 libs
	$(FC) $(FFLAGS) '-Wl,-rpath,.' -o $@ $< -ldl

test_fortran: fortran_plugins.F90 test_fortran.F90 plugin_code.o libs
	$(FC) $(FFLAGS) -I. -o $@ fortran_plugins.F90 test_fortran.F90 plugin_code.o -ldl '-Wl,-rpath,$$ORIGIN'

clean:
	rm -f *.so *.o *.mod test_fortran test_c test_dl_c test_dl_f a.out *.a
