program lorenzo_test
  use lorenzo_mod
  implicit none
#include <misc_timers.hf>
  integer, parameter :: ni = 6
  integer, parameter :: lni = 7
  integer, parameter :: nj = 5
  integer :: i, j, errors
  integer, dimension(lni, nj) :: f, fr
  integer(C_INT64_T) :: t0, t1, freq
  real(kind=4) :: ns

  freq = cycles_counter_freq()
  ns = freq
  ns = 1.0E9/freq
  f = 999
  do j = 1, nj
  do i = 1, ni
    f(i,j) = i*i + j*j
    fr(i,j) = i*i + j*j
  enddo
  enddo
  print *,"============ original ============"
  do j = nj, 1, -1
    print 1, j, f(1:lni,j)
  enddo

  t0 = elapsed_cycles()
  call lorenzo_predict_2d_inplace(f, ni, lni, nj)
  t1 = elapsed_cycles()
  print *,"============ predicted ============"," time =",(t1-t0)*ns,' ns'
  do j = nj, 1, -1
    print 1, j, f(1:lni,j)
  enddo

  t0 = elapsed_cycles()
  call lorenzo_unpredict_2d_inplace(f, ni, lni, nj)
  t1 = elapsed_cycles()
  print *,"============ restored ============"," time =",(t1-t0)*ns,' ns'
  do j = nj, 1, -1
    print 1, j, f(1:lni,j)
  enddo
  errors = 0
  do j = 1, nj
  do i = 1, ni
    if( f(i,j) .ne. fr(i,j) ) errors = errors + 1
  enddo
  enddo
  print *,'errors =',errors
1 format(I6,' |',20I6)
end program
