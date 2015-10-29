
subroutine dwt_fwd_lift_haar_r(z,ni,nj,alongx,alongy)
! in place FORWARD Haar lifting transform
! output is stored interleaved (shuffled) into array z
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
  logical, intent(IN) :: alongx, alongy

  integer :: i, j, ii, jj, jp2, jm1
!
! 1D forward lifting transform along i  (mirror condition at both ends)
! even points : "low" frequency, odd points : "high" frequency"
!
  if(alongx) then
    do j = 0 , nj-1
      z(1,j) = z(1,j) - .5*(z(0,j) + z(2,j))              ! predict odd points using even points
      z(0,j) = z(0,j) - z(1,j)                            ! update even points using residuals at odd points
      do i = 2 , ni-4 , 2
        z(i+1,j) = z(i+1,j) - .5*(z(i+2,j) + z(i  ,j))    ! predict odd points using even points
        z(i  ,j) = z(i  ,j) - .5*(z(i+1,j) + z(i-1,j))    ! update even points using residuals at odd points
      enddo
      z(ni-1,j) = z(ni-1,j) - z(ni-2,j)                   ! predict odd points using even points
      z(ni-2,j) = z(ni-2,j) - .5*(z(ni-1,j) + z(ni-3,j))  ! update even points using residuals at odd points
    enddo
  endif
!
! 1D forward lifting transform along j  (mirror condition at both ends)
! even points : "low" frequency, odd points : "high" frequency"
!
  if(alongy) then
    do j = 0 , nj - 2 , 2
      jp2 = j + 2
      if(jp2 > nj-1) jp2 = nj - 2
      jm1 = j - 1
      if(jm1 < 0) jm1 = 1
      z(:,j+1) = z(:,j+1) - .5*(z(:,jp2) + z(:,j  ))
      z(:,j  ) = z(:,j  ) - .5*(z(:,jm1) + z(:,j+1))
    enddo
  endif
  return
end subroutine dwt_fwd_lift_haar_r

subroutine dwt_inv_lift_haar_r(z,ni,nj,alongx,alongy)
! in place FORWARD inverse Haar lifting transform
! input is interleaved (shuffled) in array z
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
  logical, intent(IN) :: alongx, alongy

  integer :: i, j, ii, jj, jp2, jm1
!
! 1D inverse lifting transform along j  (mirror condition at both ends)
! even points : "low" frequency, odd points : "high" frequency"
!
  if(alongy) then
    do j = nj - 2 , 0 , -2
      jp2 = j + 2
      if(jp2 > nj-1) jp2 = nj - 2
      jm1 = j - 1
      if(jm1 < 0) jm1 = 1
      z(:,j  ) = z(:,j  ) + .5*(z(:,jm1) + z(:,j+1))      ! update even points using residuals at odd points
      z(:,j+1) = z(:,j+1) + .5*(z(:,jp2) + z(:,j  ))      ! predict odd points using even points
    enddo
  endif
!
! 1D inverse lifting transform along i  (mirror condition at both ends)
! even points : "low" frequency, odd points : "high" frequency"
!
  if(alongx) then
    do j = 0 , nj-1
      z(ni-2,j) = z(ni-2,j) + .5*(z(ni-1,j) + z(ni-3,j))  ! update even points using residuals at odd points
      z(ni-1,j) = z(ni-1,j) + z(ni-2,j)                   ! predict odd points using even points
      do i = ni-4 , 2 , -2
        z(i  ,j) = z(i  ,j) + .5*(z(i+1,j) + z(i-1,j))    ! update even points using residuals at odd points
        z(i+1,j) = z(i+1,j) + .5*(z(i+2,j) + z(i,j))      ! predict odd points using even points
      enddo
      z(0,j) = z(0,j) + z(1,j)                            ! update even points using residuals at odd points
      z(1,j) = z(1,j) + .5*(z(0,j) + z(2,j))              ! predict odd points using even points
    enddo
  endif
  return
end subroutine dwt_inv_lift_haar_r

subroutine dwt_fwd_lift_haar_i(z,ni,nj,alongx,alongy)
  implicit none
  integer, intent(IN) :: ni, nj
  integer, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
  logical, intent(IN) :: alongx, alongy

  integer :: i, j, ii, jj, jp2, jm1
!
! 1D forward lifting transform along i  (mirror condition at both ends)
!
  if(alongx) then
    do j = 0 , nj-1
      z(1,j) = z(1,j) - nint(.5*(z(0,j) + z(2,j)))              ! predict odd points using even points
      z(0,j) = z(0,j) - z(1,j)                                  ! update even points using residuals at odd points
      do i = 2 , ni-4 , 2
        z(i+1,j) = z(i+1,j) - nint(.5*(z(i+2,j) + z(i  ,j)))    ! predict odd points using even points
        z(i  ,j) = z(i  ,j) - nint(.5*(z(i+1,j) + z(i-1,j)))    ! update even points using residuals at odd points
      enddo
      z(ni-1,j) = z(ni-1,j) - z(ni-2,j)                         ! predict odd points using even points
      z(ni-2,j) = z(ni-2,j) - nint(.5*(z(ni-1,j) + z(ni-3,j)))  ! update even points using residuals at odd points
    enddo
  endif
!
! 1D forward lifting transform along j  (mirror condition at both ends)
!
  if(alongy) then
    do j = 0 , nj - 2 , 2
      jp2 = j + 2
      if(jp2 > nj-1) jp2 = nj - 2
      jm1 = j - 1
      if(jm1 < 0) jm1 = 1
      z(:,j+1) = z(:,j+1) - nint(.5*(z(:,jp2) + z(:,j  )))
      z(:,j  ) = z(:,j  ) - nint(.5*(z(:,jm1) + z(:,j+1)))
    enddo
  endif
  return
end subroutine dwt_fwd_lift_haar_i

subroutine dwt_inv_lift_haar_i(z,ni,nj,alongx,alongy)
  implicit none
  integer, intent(IN) :: ni, nj
  integer, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
  logical, intent(IN) :: alongx, alongy

  integer :: i, j, ii, jj, jp2, jm1
!
! 1D inverse lifting transform along j  (mirror condition at both ends)
!
  if(alongy) then
    do j = nj - 2 , 0 , -2
      jp2 = j + 2
      if(jp2 > nj-1) jp2 = nj - 2
      jm1 = j - 1
      if(jm1 < 0) jm1 = 1
      z(:,j  ) = z(:,j  ) + nint(.5*(z(:,jm1) + z(:,j+1)))      ! update even points using residuals at odd points
      z(:,j+1) = z(:,j+1) + nint(.5*(z(:,jp2) + z(:,j  )))      ! predict odd points using even points
    enddo
  endif
!
! 1D inverse lifting transform along i  (mirror condition at both ends)
!
  if(alongx) then
    do j = 0 , nj-1
      z(ni-2,j) = z(ni-2,j) + nint(.5*(z(ni-1,j) + z(ni-3,j)))  ! update even points using residuals at odd points
      z(ni-1,j) = z(ni-1,j) + z(ni-2,j)                         ! predict odd points using even points
      do i = ni-4 , 2 , -2
        z(i  ,j) = z(i  ,j) + nint(.5*(z(i+1,j) + z(i-1,j)))    ! update even points using residuals at odd points
        z(i+1,j) = z(i+1,j) + nint(.5*(z(i+2,j) + z(i,j)))      ! predict odd points using even points
      enddo
      z(0,j) = z(0,j) + z(1,j)                                  ! update even points using residuals at odd points
      z(1,j) = z(1,j) + nint(.5*(z(0,j) + z(2,j)))              ! predict odd points using even points
    enddo
  endif
  return
end subroutine dwt_inv_lift_haar_i

subroutine dwt_unshuffle(zs,zu,ni,nj)  ! 2D unshuffling (de-interleaving) of zs into zu
! input array zs is interlaved (shuffled)
! output array zu is in "quadrant" form
  implicit none
  integer, intent(IN) :: ni, nj
  integer, dimension(0:ni-1, 0:nj-1), intent(IN)  :: zs
  integer, dimension(0:ni-1, 0:nj-1), intent(OUT) :: zu

  integer :: ni2, nj2, i, j, ii, jj

  ni2 = ni/2
  nj2 = nj/2
  do j=0,nj2-1
    jj = j + j
    ii = 0
    do i=0,ni2-1
      zu(i    ,j    ) = zs(ii  ,jj  )  ! LL quadrant
      zu(ni2+i,j    ) = zs(ii+1,jj  )  ! HL quadrant
      zu(i    ,nj2+j) = zs(ii  ,jj+1)  ! LH quadrant
      zu(ni2+i,nj2+j) = zs(ii+1,jj+1)  ! HH quadrant
      ii = ii + 2
    enddo
  enddo
  return
end subroutine dwt_unshuffle

subroutine dwt_shuffle(zs,zu,ni,nj)  ! 2D shuffling (interleaving) of zu into zs
! input array zu is in "quadrant" form
! output array zs is interlaved (shuffled)
  implicit none
  integer, intent(IN) :: ni, nj
  integer, dimension(0:ni-1, 0:nj-1), intent(IN)  :: zu
  integer, dimension(0:ni-1, 0:nj-1), intent(OUT) :: zs

  integer :: ni2, nj2, i, j, ii, jj

  ni2 = ni/2
  nj2 = nj/2
  do j=0,nj2-1
    jj = j + j
    ii = 0
    do i=0,ni2-1
      zs(ii  ,jj  ) = zu(i    ,j    )   ! LL quadrant
      zs(ii+1,jj  ) = zu(ni2+i,j    )   ! HL quadrant
      zs(ii  ,jj+1) = zu(i    ,nj2+j)   ! LH quadrant
      zs(ii+1,jj+1) = zu(ni2+i,nj2+j)   ! HH quadrant
      ii = ii + 2
    enddo
  enddo
  return
end subroutine dwt_shuffle

subroutine dwt_qsplit(zs,ni,nj,ll,hl,lh,hh,nni,nnj)  ! split interleaved zs into quadrants
! split interleaved array zs into 4 quadrants (ll, lh, hl, hh)
! it is expected that ni >= 2*nni and nj >= 2*nnj
  implicit none
  integer, intent(IN) :: ni, nj, nni, nnj
  integer, dimension(0:ni-1, 0:nj-1), intent(IN) :: zs
  integer, dimension(0:nni-1, 0:nnj-1), intent(OUT) :: ll, lh, hl, hh

  integer :: i, j, ii, jj

  do j=0,nnj-1
    jj = j + j
    ii = 0
    do i=0,nni-1
      ll(i,j) = zs(ii  ,jj  )  ! LL quadrant
      hl(i,j) = zs(ii+1,jj  )  ! HL quadrant
      lh(i,j) = zs(ii  ,jj+1)  ! LH quadrant
      hh(i,j) = zs(ii+1,jj+1)  ! HH quadrant
      ii = ii + 2
    enddo
  enddo
  return
end subroutine dwt_qsplit

subroutine dwt_qmerge(zs,ni,nj,ll,hl,lh,hh,nni,nnj)  ! interleave  quadrants into zs
! merge and interleave 4 quadrants (ll, lh, hl, hh) into zs
! it is expected that ni >= 2*nni and nj >= 2*nnj
  implicit none
  integer, intent(IN) :: ni, nj, nni, nnj
  integer, dimension(0:nni-1, 0:nnj-1), intent(IN) :: ll, lh, hl, hh
  integer, dimension(0:ni-1, 0:nj-1), intent(OUT) :: zs

  integer :: i, j, ii, jj

  do j=0,nnj-1
    jj = j + j
    ii = 0
    do i=0,nni-1
      zs(ii  ,jj  ) = ll(i,j)   ! LL quadrant
      zs(ii+1,jj  ) = hl(i,j)   ! HL quadrant
      zs(ii  ,jj+1) = lh(i,j)   ! LH quadrant
      zs(ii+1,jj+1) = hh(i,j)   ! HH quadrant
      ii = ii + 2
    enddo
  enddo
  return
end subroutine dwt_qmerge

subroutine dwt_normalize(z,n,bigval,auto)
! simulate normalization to largest exponent of field z by
! dropping some least significant bits if exponent of number
! is smaller than bigval's (or largest value of field)
! if exponent is n less than normalization exponent, drop n LSBs
! if n > 23, set number to zero
  implicit none
  integer, intent(IN) :: n
  real, intent(IN) :: bigval
  real, dimension(n), intent(INOUT) :: z
  logical :: auto

  integer :: i, j, maxexp, myexp, ival
  real :: temp
! transfer real into integer
  ival = transfer(bigval,ival)
  if(auto) then
    temp = max( abs(maxval(z)) , abs(minval(z)) )
    ival = transfer(temp,ival)
  endif
! extract exponent (lower 8 of upper 9 bits for IEE754-32)
  maxexp = iand(255 , ishft(ival , -23) )   ! normalization exponent
  do i = 1 , n
    ival = transfer(z(i),ival)
    myexp = iand(255 , ishft(ival , -23) )
    if(maxexp-myexp > 0) then      ! nothing to drop if myexp >= maxexp
      if (maxexp-myexp > 23) then  ! drop maxexp-myexp least significant bits
        ival = 0                                    ! nothing left
      else
        ival = iand(ival , ishft(-1,maxexp-myexp))  ! drop LSBs
      endif
    endif
    z(i) = transfer(ival,z(i))
  enddo
end subroutine dwt_normalize

subroutine dwt_quantize(z,iz,n,nbts,error)
  implicit none
  integer, intent(IN) :: n
  integer, intent(IN) :: nbts
  real, intent(IN) :: error
  real, dimension(n), intent(IN) :: z
  integer, dimension(n), intent(OUT) :: iz
  integer :: nbits

  integer :: i, j, maxexp, myexp, ival, hidden, mask, ratio
  real :: temp, tempmin

! transfer real into integer
  tempmin = minval(z)
!  if(auto) then
    temp = maxval(z) - tempmin
    ival = transfer(temp,ival)
!  endif
! extract exponent (lower 8 of upper 9 bits for IEE754-32)
  maxexp = iand(255 , ishft(ival , -23) )   ! exponent used for forced normalization

  nbits = nbts
  if(nbits == -1) then
    nbits = 0
    ratio = nint(temp/error)
    do while(ratio > 0)
      nbits = nbits + 1
      ratio = ishft(ratio,-1)
    enddo
    print *,'AUTO: max-min, error, nbits=',temp,temp/error,nbits,ratio
  endif
  hidden = ishft(1,23)                       ! "hidden" 1 of IEEE754-32
  mask = not(ishft(-1,23))                   ! mask for mantissa (lower 23 bits)
  do i = 0 , n
    temp = z(i) - tempmin
    ival = transfer(temp,ival)               ! transfer real into integer
    myexp = iand(255 , ishft(ival , -23) )   ! extract exponent
    if (maxexp-myexp > 23) then
      ival = 0
    else
      ival = ior(hidden,iand(mask,ival))     ! get mantissa, add "hidden" 1
      ival = ishft(ival ,myexp-maxexp)       ! quantize absolute value
    endif
    ival = ishft(ival, nbits - 24)             ! reduce to nbits
!    if(z(i) < 0) ival = -ival                ! take care of sign
    iz(i) = ival
  enddo
end subroutine dwt_quantize
