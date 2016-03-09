#!/bin/bash
set -x
rm -f *.o a.out 
s.cc -c -Dtest_main=MY_C_MAIN c_main.c 
s.f90 f_main_to_cmain.F90 c_main.o 
./a.out a b c d
rm -f *.o a.out 
