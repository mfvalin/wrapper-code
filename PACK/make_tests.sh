#!/bin/bash
set -x
for mod in aocc intel gnu ; do
  make OPT='-O3 -mavx2 -DWITH_SIMD' test_floats MOD=$mod
  mv test_floats.$mod.exe SIMD/.
done
#
for mod in aocc intel gnu ; do
  make OPT='-O3 -mavx2' test_floats MOD=$mod
  mv test_floats.$mod.exe NO_SIMD/.
done
