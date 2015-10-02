
subroutine dwt_fwd_lift_haar_r(z,ni,nj,alongx,alongy)
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
  logical, intent(IN) :: alongx, alongy

  integer :: i, j, ii, jj, jp2, jm1
!
! 1D forward lifting transform along i  (mirror condition at both ends)
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
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
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
      z(:,j  ) = z(:,j  ) + .5*(z(:,jm1) + z(:,j+1))      ! update even points using residuals at odd points
      z(:,j+1) = z(:,j+1) + .5*(z(:,jp2) + z(:,j  ))      ! predict odd points using even points
    enddo
  endif
!
! 1D inverse lifting transform along i  (mirror condition at both ends)
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

subroutine dwt_unshuffle(zs,zu,ni,nj)  ! unshuffle zs into zu
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

subroutine dwt_shuffle(zs,zu,ni,nj)  ! shuffle zu into zs
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

subroutine dwt_qsplit(zs,ni,nj,ll,hl,lh,hh,nni,nnj)  ! split zs into quadrants
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

subroutine dwt_qmerge(zs,ni,nj,ll,hl,lh,hh,nni,nnj)  ! inerleave  quadrants into zs
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

subroutine dwt_normalize(z,ni,nj,bigval,auto)
  implicit none
  integer, intent(IN) :: ni, nj
  real, intent(IN) :: bigval
  real, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
  logical :: auto

  integer :: i, j, maxexp, myexp, ival
  real :: temp

  ival = transfer(bigval,ival)
  if(auto) then
    temp = max( abs(maxval(z)) , abs(minval(z)) )
    ival = transfer(temp,ival)
  endif
  maxexp = iand(255 , ishft(ival , -23) )
  do j = 0 , nj-1
  do i = 0 , ni-1
    ival = transfer(z(i,j),ival)
    myexp = iand(255 , ishft(ival , -23) )
    if (maxexp-myexp > 23) then
      ival = 0
    else
      ival = iand(ival , ishft(-1,maxexp-myexp))
    endif
    z(i,j) = transfer(ival,z(i,j))
  enddo
  enddo
end subroutine dwt_normalize

subroutine dwt_quantize(z,iz,ni,nj,bigval,auto)
  implicit none
  integer, intent(IN) :: ni, nj
  real, intent(IN) :: bigval
  real, dimension(0:ni-1, 0:nj-1), intent(IN) :: z
  integer, dimension(0:ni-1, 0:nj-1), intent(OUT) :: iz
  logical :: auto

  integer :: i, j, maxexp, myexp, ival, hidden, mask
  real :: temp

  hidden = ishft(1,23)
  mask = ishft(-1,23)
  ival = transfer(bigval,ival)
  if(auto) then
    temp = maxval( abs(z) )
    ival = transfer(temp,ival)
  endif
  maxexp = iand(255 , ishft(ival , -23) )   ! exponent used for forced normalization

  do j = 0 , nj-1
  do i = 0 , ni-1
    ival = transfer(z(i,j),ival)
    myexp = iand(255 , ishft(ival , -23) )
    if (maxexp-myexp > 23) then
      ival = 0
    else
      ival = iand(mask,ival)
      ival = ishft( ior(hidden,ival) ,myexp-maxexp)
    endif
    if(z(i,j) < 0) ival = -ival
    iz(i,j) = ival
  enddo
  enddo
end subroutine dwt_quantize
