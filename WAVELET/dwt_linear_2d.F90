#define ARGUMENT_TYPE real*4
subroutine inv_linear53_dwt2d_r4(z,ni,nj)   ! 2D inverse transform , along j first, then along i
! in place INVERSE lifting transform using linear prediction wavelet for ARGUMENT_TYPE numbers
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  ARGUMENT_TYPE, intent(INOUT), dimension(0:ni-1,0:nj-1) :: z

  ARGUMENT_TYPE, dimension(-1:ni) :: even, odd
  integer :: i, j, j00, jm1, jm2, jp1, jj
  integer :: nodd, neven

  nodd  = ishft(ni,-1)     ! number of odd terms
  neven = ishft(ni+1,-1)   ! number of even terms (nodd + 1 if ni is odd)
  j00 = 0

  do j = 0 , nj-1 , 2
    jm2 = j - 2           ! not used for j == 0
    jm1 = abs(j - 1)      ! mirror condition at lower boundary
    jp1 = j + 1
    if(jp1 == nj) jp1 = nj - 2                     ! upper mirror boundary condition
    z(:,j) = z(:,j) - .25 * (z(:,jm1) + z(:,jp1))  ! unupdate even row j
    if(j > 0) then
      z(:,jm1) = z(:,jm1) + .5 * (z(:,jm2) + z(:,j)) ! unpredict odd row jm1
      call inv_linear53_dwt1d_r4(z(:,j-2))
      call inv_linear53_dwt1d_r4(z(:,j-1))
      j00 = j
    endif
  enddo
  if(iand(nj,1)==0) z(:,nj-1) = z(:,nj-1) + z(:,nj-2) ! last row is odd, unpredict it
  do j = j00 , nj-1
    call inv_linear53_dwt1d_r4(z(:,j))
  enddo
contains
  subroutine inv_linear53_dwt1d_r4(f)  ! 1D along i transform
    use ISO_C_BINDING
    implicit none
    ARGUMENT_TYPE, intent(INOUT), dimension(0:ni-1) :: f

    integer :: i

    do i=0,nodd-1           ! split into even / odd
      even(i) = f(i)
      odd(i)  = f(neven+i)
    enddo
    if(iand(ni,1) .ne. 0) even(nodd) = f(nodd)   ! one more even values than odd values if ni is odd
    odd(-1)   = odd(0)         ! mirror condition at lower boundary
    odd(nodd) = odd(nodd-1)    ! mirror condition at upper boundary, used only of ni is odd

    even(0) = even(0) - odd(0)
    do i = 1,neven-1  ! unupdate even points
      even(i) = even(i) - .25 * (odd(i-1) + odd(i))
    enddo

    if(iand(ni,1) .eq. 0) even(nodd) = even(nodd-1) ! mirror condition at upper boundary if ni is even
    do i = 0,nodd-1  ! unpredict odd points
      odd(i) = odd(i) + .5 * (even(i) + even(i+1))
    enddo
    do i = 0,nodd-1
      f(2*i)   = even(i)
      f(2*i+1) = odd(i)
    enddo
    if(iand(ni,1) .ne. 0) f(ni-1) = even(neven-1)  ! odd number of points, one extra even element
  end subroutine inv_linear53_dwt1d_r4
end subroutine inv_linear53_dwt2d_r4

subroutine fwd_linear53_dwt2d_r4(z,ni,nj)   ! 2D forward transform, along i first, then along j
! in place FORWARD lifting transform using linear prediction wavelet for ARGUMENT_TYPE numbers
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  ARGUMENT_TYPE, intent(INOUT), dimension(0:ni-1,0:nj-1) :: z

  ARGUMENT_TYPE, dimension(-1:ni) :: even, odd
  integer :: i, j, j00, jm1, jm2, jp1, jj
  integer :: nodd, neven

  nodd  = ishft(ni,-1)     ! number of odd terms
  neven = ishft(ni+1,-1)   ! number of even terms (nodd + 1 if ni is odd)
  j00 = 0
  do j = 1,nj-1,2
    jm2 = abs(j - 2)
    jm1 = j - 1
    jp1 = j + 1
    if(jp1 == nj) jp1 = nj - 2                              ! upper mirror boundary condition
    do jj = j00,max(j,jp1)
      call fwd_linear53_dwt1d_r4(z(:,jj))
    enddo
    j00 = jp1 + 1
    z(:,  j) = z(:,  j) -  .5 * (z(:,jm1) + z(:,jp1))       ! predict odd rows
    z(:,jm1) = z(:,jm1) + .25 * (z(:,jm2) + z(:,  j))       ! update even rows (below odd row)
  enddo
  if(mod(nj,2)==1) z(:,nj-1) = z(:,nj-1) + .5 * z(:,nj-2)   ! odd number of rows, last row is en even row
  return
contains
  subroutine fwd_linear53_dwt1d_r4(f)   ! 1D along i transform
    implicit none
    ARGUMENT_TYPE, intent(INOUT), dimension(0:ni-1) :: f

    integer :: i

    do i=0,nodd-1           ! split into even / odd
      even(i) = f(2*i)
      odd(i)  = f(2*i+1)
    enddo
    if(iand(ni,1) .ne. 0) then  ! is ni odd ?
      even(nodd) = f(ni-1)      ! one more even values than odd values if ni is odd
    else
      even(nodd) = even(nodd-1) ! mirror condition at upper boundary if ni is even
    endif
    do i = 0, nodd-1           ! predict odd values (and copy updated values into f)
      odd(i) = odd(i) - .5 * (even(i) + even(i+1))
      f(neven+i) = odd(i)
    enddo
    odd(-1)   = odd(0)         ! mirror condition at lower boundary
    odd(nodd) = odd(nodd-1)    ! mirror condition at upper boundary, used only of ni is odd
    do  i=0,neven-1            ! update even values (and copy updated values into f)
      f(i) = even(i) + .25 * (odd(i) + odd(i-1))
    enddo
  end subroutine fwd_linear53_dwt1d_r4
end subroutine fwd_linear53_dwt2d_r4
#if defined(SELF_TEST)
#define NI 12001
#define NJ 9001
program test
  ARGUMENT_TYPE, dimension(0:NI-1,0:NJ-1) :: z, z0
  integer :: i, j
  real*8 :: T0,T1,T2,T3
  real*8, external :: MPI_WTIME
  do j = 0,NJ-1
  do i = 0,NI-1
    z(i,j) = 1 + i*1.2 + j*1.7
    z0(i,j) = z(i,j)
  enddo
  enddo
  print *,'min,max z0',minval(z0),maxval(z0)
1 format(12F6.2)
  if(NI <=10 .and. NJ <=10) then
    print *,'=========================================================='
    do j = NJ-1,0,-1
      print 1,z(:,j)
    enddo
  endif
  do irep = 1, 100
  T0 = MPI_WTIME()
  call fwd_linear53_dwt2d_r4(z,NI,NJ)
  T1 = MPI_WTIME()
  if(NI <=10 .and. NJ <=10) then
    print *,'=========================================================='
    do j = NJ-1,0,-1
      print 1,z(:,j)
    enddo
  endif
  T2 = MPI_WTIME()
  call inv_linear53_dwt2d_r4(z,NI,NJ)
  T3 = MPI_WTIME()
  enddo
  if(NI <=10 .and. NJ <=10) then
    print *,'=========================================================='
    do j = NJ-1,0,-1
      print 1,z(:,j)
    enddo
  endif
  print *,'transform time:',t1-t0,t3-t2,(t1-t0)/NI/NJ,(t3-t2)/NI/NJ
  print *,'min,max z0',minval(z0),maxval(z0)
  print *,'min,max z',minval(z),maxval(z)
  z0 = (z0-z)/z0
  print *,'total error:',sum(z0)/NI/NJ
  stop
end program test
#endif

