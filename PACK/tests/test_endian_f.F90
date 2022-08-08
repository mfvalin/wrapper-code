program test_endian_f
  use ISO_C_BINDING
  implicit none
#include <misc_endian.hf>
#include <misc_timers.hf>
  integer, parameter :: NTIMES = 10
  integer(C_INT64_T) :: t0, t1, tmin, t(NTIMES*100), tmax, freq
  integer(C_INT64_T), dimension(4096*4096*2) :: s, d
  integer :: i, isize
  real(kind=8) :: nano

  freq = cycles_counter_freq()
  nano = 1.0E+9 / freq
  s = 0
  d = 0
  isize = size(s)
  call Swap_8_in_32(s, d, isize)
  call Swap_8_in_32(s, d, isize)
  tmin = 999999999
  tmax = 0
  do i = 1, NTIMES
    t0 = elapsed_cycles()
    call Swap_8_in_32(s, d, isize)
    t1 = elapsed_cycles()
    t(i) = t1 - t0
    tmin = min(tmin, t(i))
    tmax = max(tmax, t(i))
  enddo
  print *,'in memory'
  print *,'tmin, tmax =',tmin,tmax
  print 1,'tmin/word, tmax/word ', tmin*nano/isize,tmax*nano/isize

  print *,'in cache'
  tmin = 999999999
  tmax = 0
  isize = 4096
  do i = 1, NTIMES*100
    t0 = elapsed_cycles()
    call Swap_8_in_32(s, d, isize)
    t1 = elapsed_cycles()
    t(i) = t1 - t0
    tmin = min(tmin, t(i))
    tmax = max(tmax, t(i))
  enddo
  print *,'tmin, tmax =',tmin,tmax
  print 1,'tmin/word, tmax/word ', tmin*nano/isize,tmax*nano/isize
1 format(A,2F8.2)
end
