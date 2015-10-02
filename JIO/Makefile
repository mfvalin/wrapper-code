
linux64: ./s.jio-prof-$(BASE_ARCH).so

linux32: ./s.jio-prof-Linux.so

./s.jio-prof-$(BASE_ARCH).so: jio-uni-2.c
#	gcc      -fpic -shared          -o ./s.jio-prof-$(BASE_ARCH).so jio-uni-2.c -ldl
	gcc      -fpic -shared -pthread -o ./s.jio-prof-$(BASE_ARCH).so jio-uni-2.c -ldl

./s.jio-prof-Linux.so: jio-uni-2.c
#	gcc -m32 -fpic -shared          -o ./s.jio-prof-Linux.so jio-uni-2.c -ldl
	gcc -m32 -fpic -shared -pthread -o ./s.jio-prof-Linux.so jio-uni-2.c -ldl

AIX: ./s.jio-prof-libjiop7.a

./jio-uni-ibmp7.so: jio-uni-2.c
	s.cc -openmp -qcpluscmt -shared jio-uni-2.c -o jio-uni-ibmp7.so -lm

s.jio-prof-libjiop7.a: jio-uni-ibmp7.so
	ar rcv  s.jio-prof-libjiop7.a jio-uni-ibmp7.so

linux: linux32 linux64

clean:
	rm -f *.so *.o

install: 
	echo install

