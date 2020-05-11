program test_cdf
  use ISO_C_BINDING
  implicit none
#define FORTRAN_SOURCE
#include <cdf97.h>
#include <cdf53.h>
#if ! defined(NPTS)
#define NPTS 15
#endif
  real, dimension(NPTS,NPTS) :: x
  integer, dimension(NPTS,NPTS) :: xi
  integer :: i, j
  real :: quantum = .05

  do j = 0, NPTS-1
  do i = 0, NPTS-1
    x(i+1,j+1) = (3+i+0.4*i*i-0.02*i*i*i) * (3+j+0.4*j*j-0.02*j*j*j)
  enddo
  enddo
  print *,'Original'
  do j = NPTS, 1, -1
    print 1,x(:,j)
  enddo
  call F_CDF97_2D_split_inplace_n(x, NPTS,  NPTS, NPTS, 3)
  print *,'After transform and quantification'
  xi = (x / quantum) + .5
  x = xi * quantum
  do j = NPTS, 1, -1
    print 1,x(:,j)
  enddo
  call I_CDF97_2D_split_inplace_n(x, NPTS,  NPTS, NPTS, 3)
  print 2,'Error after quantification by',quantum
  do j = NPTS-1, 0, -1
    do i = 0, NPTS-1
      x(i+1,j+1) = x(i+1,j+1) - (3+i+0.4*i*i-0.02*i*i*i) * (3+j+0.4*j*j-0.02*j*j*j)
    enddo
    print 1,x(:,j+1)
    enddo
1 format(20F9.2)
2 format(A,F9.3)
end
