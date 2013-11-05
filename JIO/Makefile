
linux64: ./s.jio-prof-$(BASE_ARCH).so

linux32: ./s.jio-prof-Linux.so

./s.jio-prof-$(BASE_ARCH).so: jio-uni-2.c
	s.cc -shared -o $@ -openmp jio-uni-2.c -ldl

./s.jio-prof-Linux.so: jio-uni-2.c
	cc -m32 -fpic -shared -o ./s.jio-prof-Linux.so jio-uni-2.c -ldl

linux: linux32 linux64

clean:
	rm -f *.so *.o

install: 
	echo install

