// Hopefully useful code for C
// Copyright (C) 2022  Recherche en Prevision Numerique
//
// This code is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation,
// version 2.1 of the License.
//
// This code is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
#include <misc_timers.h>

uint64_t ElapsedUs(void){
  return elapsed_us() ;
}

uint64_t ElapsedCycles(void){
  return elapsed_cycles() ;
}

uint64_t CyclesCounterFreq(void){
  return cycles_counter_freq() ;
}
