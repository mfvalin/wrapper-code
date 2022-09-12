program lorenzo_test
  use lorenzo_mod
  use ISO_C_BINDING
  implicit none
#include <misc_timers.hf>

#if ! defined(NI)
#define NI 128
#endif
#if ! defined(NJ)
#define NJ NI
#endif
#define NTIMES 10000

  integer, parameter :: ni = NI
  integer, parameter :: lnio = ni+2
  integer, parameter :: lnid = ni+1
  integer, parameter :: nj = NJ
  integer, parameter :: np = NI*NJ
  integer :: i, j, k, errors
  integer, dimension(lnio, nj) :: f, fr, fp
  integer, dimension(lnid, nj) :: pred
  integer(C_INT64_T) :: t0, t1, freq
  integer(C_INT64_T), dimension(NTIMES) :: t
  real(kind=4) :: ns, avg

  freq = cycles_counter_freq()
  ns = freq
  ns = 1.0E9/freq
  f  = 888
  fr = 999
  pred = 999
  do j = 1, nj
  do i = 1, ni
!     f(i,j) = i*i + j*j
    fr(i,j) = i*i + j*j
  enddo
  enddo
  f = fr
  if(ni < 18 .and. nj < 10) then
    print 2,"=== original ==="
    do j = nj, 1, -1
      print 1, j, f(1:lnio,j)
    enddo
  endif

  do k = 1, NTIMES
    f = fr
    t(k) = elapsed_cycles()
    call lorenzopredict_inplace(f, ni, lnio, nj)
    t(k) = elapsed_cycles() - t(k)
  enddo
  avg = SUM(t*ns)/NTIMES
  print 2,"=== pred in place ==="," time =",avg/np,' ns/pt (min =',MINVAL(t*ns)/np,')'
  if(ni < 18 .and. nj < 10) then
    do j = nj, 1, -1
      print 1, j, f(1:lnio,j)
    enddo
  endif

  f = fr
  do k = 1, NTIMES
    t(k) = elapsed_cycles()
    call lorenzopredict(f, pred, ni, lnio, lnid, nj)
    t(k) = elapsed_cycles() - t(k)
  enddo
  avg = SUM(t*ns)/NTIMES
  print 2,"=== generic       ==="," time =",avg/np,' ns/pt (min =',MINVAL(t*ns)/np,')'
  if(ni < 18 .and. nj < 10) then
    do j = nj, 1, -1
      print 1, j, pred(1:lnid,j)
    enddo
  endif

  f = fr
  do k = 1, NTIMES
    t(k) = elapsed_cycles()
    call lorenzopredict_f(f, pred, ni, lnio, lnid, nj)
    t(k) = elapsed_cycles() - t(k)
  enddo
  avg = SUM(t*ns)/NTIMES
  print 2,"=== predicted F90 ==="," time =",avg/np,' ns/pt (min =',MINVAL(t*ns)/np,')'
  if(ni < 18 .and. nj < 10) then
    do j = nj, 1, -1
      print 1, j, pred(1:lnid,j)
    enddo
  endif

  f = fr
  do k = 1, NTIMES
    t(k) = elapsed_cycles()
    call lorenzopredict_c(f, pred, ni, lnio, lnid, nj)
    t(k) = elapsed_cycles() - t(k)
  enddo
  avg = SUM(t*ns)/NTIMES
  print 2,"=== predicted C   ==="," time =",avg/np,' ns/pt (min =',MINVAL(t*ns)/np,')'
  if(ni < 18 .and. nj < 10) then
    do j = nj, 1, -1
      print 1, j, pred(1:lnid,j)
    enddo
  endif

  do k = 1, NTIMES
    fp = fr
    call lorenzopredict_inplace(fp, ni, lnio, nj)
    t(k) = elapsed_cycles()
    call lorenzounpredict_inplace(fp, ni, lnio, nj)
    t(k) = elapsed_cycles() - t(k)
  enddo
  avg = SUM(t*ns)/NTIMES
  print 2,"=== in place F90  ==="," time =",avg/np,' ns/pt (min =',MINVAL(t*ns)/np,')'
  if(ni < 18 .and. nj < 10) then
    do j = nj, 1, -1
      print 1, j, fp(1:lnio,j)
    enddo
  endif
  errors = 0
  do j = 1, nj
  do i = 1, ni
    if( fp(i,j) .ne. fr(i,j) ) errors = errors + 1
  enddo
  enddo
  print *,'errors =',errors

  f = 999
  do k = 1, NTIMES
    t(k) = elapsed_cycles()
    call lorenzounpredict(f, pred, ni, lnio, lnid, nj)
    t(k) = elapsed_cycles() - t(k)
  enddo
  avg = SUM(t*ns)/NTIMES
  print 2,"=== restored F90  ==="," time =",avg/np,' ns/pt (min =',MINVAL(t*ns)/np,')'
  if(ni < 18 .and. nj < 10) then
    do j = nj, 1, -1
      print 1, j, f(1:lnio,j)
    enddo
  endif
  errors = 0
  do j = 1, nj
  do i = 1, ni
    if( f(i,j) .ne. fr(i,j) ) errors = errors + 1
  enddo
  enddo
  print *,'errors =',errors

  f = 999
  do k = 1, NTIMES
    t(k) = elapsed_cycles()
    call lorenzounpredict_c(f, pred, ni, lnio, lnid, nj)
    t(k) = elapsed_cycles() - t(k)
  enddo
  avg = SUM(t*ns)/NTIMES
  print 2,"=== restored C    ==="," time =",avg/np,' ns/pt (min =',MINVAL(t*ns)/np,')'
  if(ni < 18 .and. nj < 10) then
    do j = nj, 1, -1
      print 1, j, f(1:lnio,j)
    enddo
  endif
  errors = 0
  do j = 1, nj
  do i = 1, ni
    if( f(i,j) .ne. fr(i,j) ) errors = errors + 1
  enddo
  enddo
  print *,'errors =',errors
1 format(I6,' |',20I6)
2 format(A,A,F6.2,A,F6.2,A)
end program
