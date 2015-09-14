program fst_wav
call fst_wav_1(256,256)
stop
end
subroutine fst_wav_1(ni,nj)
  implicit none
  integer, intent(IN) :: ni,nj

  real, dimension(0:ni-1,0:nj-1) :: z, zw, z0, zd
  real, dimension(0:ni/2-1,0:nj/2-1) :: zh, zh1, zh2
  real, dimension(0:ni/2-1,0:nj/2-1) :: ll, hl, lh, hh
  real, dimension(0:ni/4-1,0:nj/4-1) :: ll2, hl2, lh2, hh2
  real :: mult
  integer :: i, j, iblk, jblk, iun, status, nerr
  integer :: blk=64
  integer, external :: fnom, fstouv
  logical zer
  real, dimension(8) :: reals, reals2

  reals = [-1.0, 1.1, 1.01, -1.01, 1.0001, 1.00001, 1.000001, 1.0000001]

  goto 2
  reals2 = reals
  call dwt_normalize(reals2,8,1,1.1,.true.)
  print 1,1.1, reals,reals2
  reals(1) =-2.0; reals2 = reals ; 
  call dwt_normalize(reals2,8,1,2.0,.true.)
  print 1,2.0, reals,reals2
  reals(1) =-8.0; reals2 = reals ; 
  call dwt_normalize(reals2,8,1,8.0,.true.)
  print 1,8.0, reals,reals2
  reals(1) =-64.0; reals2 = reals ; 
  call dwt_normalize(reals2,8,1,64.0,.true.)
  print 1,64.0, reals,reals2
  reals(1) =-1024.0; reals2 = reals ; 
  call dwt_normalize(reals2,8,1,1024.0,.true.)
  print 1,1024.0, reals,reals2
1 format(/9G15.8/,15x,8G15.8)

  stop
2 continue

  blk=64
  zer = .false.
  z = 0.0
  mult = 4.0
  do j = 0, nj-1
    jblk = j/blk
    do i = 0, ni-1
      iblk = i/blk
      if( (mod(iblk,2) + mod(jblk,2)) == 1 ) then
        z(i,j) = 1.0
        if( mod(i,blk) >16 .and. mod(i,blk)<48 .and. mod(j,blk) >16 .and. mod(j,blk)<48 ) z(i,j) = mult
        if( mod(i,blk) >24 .and. mod(i,blk)<40 .and. mod(j,blk) >24 .and. mod(j,blk)<40 ) z(i,j) = mult*mult
      endif
      z(i,j) = sqrt( (i-ni*.5+.1)**2 + (j-nj*.5+.3)**2 )
      z(i,j) = sin(z(i,j)*.1) + sin(z(i,j)*.3) ! + sin(z(i,j)*1.7) + z(i,j)*.1  + mod(i+j,2)*.75
    enddo
  enddo
  z0 = z
  iun = 0
  status = fnom(iun,'data.fst','STD+RND',0)
  print *,'iun=',iun,' status=',status
  status = fstouv(iun,'RND')
  print *,'iun=',iun,' status=',status

  call fstecr(z,zw,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','ORIG','WAVELET','X',0,0,0,0,5,.false.)

  zw = z
  call dwt_diag(zw,ni,nj,16,'zw( pre)')
  call dwt_fwd_lift_haar_r(zw,ni,nj,.true.,.true.)              ! forward 2D transform of zw
!  call dwt_fwd_lift_haar_r(zw,ni,nj,.false.,.true.)
  call fstecr(zw,zw,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','FWH0','WAVELET','X',0,0,0,0,5,.false.)
  call dwt_qsplit(zw,ni,nj,ll,lh,hl,hh,ni/2,nj/2)

  call dwt_diag(zw,ni,nj,16,'zw( dwt)')
  call dwt_diag(ll,ni/2,nj/2,16,'ll   -')

  call dwt_fwd_lift_haar_r(ll,ni/2,nj/2,.true.,.true.)          ! forward 2D transform of ll
  call dwt_qsplit(ll,ni/2,nj/2,ll2,lh2,hl2,hh2,ni/4,nj/4)
  call dwt_diag(ll,ni/2,nj/2,16,'ll( dwt)-')

  call dwt_diag(ll2,ni/4,nj/4,16,'ll2  -')
  call dwt_diag(hl2,ni/4,nj/4,16,'hl2  -')
  call dwt_diag(lh2,ni/4,nj/4,16,'lh2  -')
  call dwt_diag(hh2,ni/4,nj/4,12,'hh2  -')

  call dwt_quant(ll2,ni/4,nj/4,18)
  call dwt_quant(hl2,ni/4,nj/4,16)
  call dwt_quant(lh2,ni/4,nj/4,16)
  call dwt_quant(hh2,ni/4,nj/4,8)

  call dwt_qmerge(ll,ni/2,nj/2,ll2,lh2,hl2,hh2,ni/4,nj/4)
  call dwt_inv_lift_haar_r(ll,ni/2,nj/2,.true.,.true.)          ! inverse 2D transform of ll

  call dwt_diag(hl,ni/2,nj/2,12,'hl   -')
  call dwt_diag(lh,ni/2,nj/2,12,'lh   -')
  call dwt_diag(hh,ni/2,nj/2,10,'hh   -')

!   call dwt_quant(ll,ni/2,nj/2,16)
  call dwt_quant(lh,ni/2,nj/2,12)
  call dwt_quant(hl,ni/2,nj/2,12)
  call dwt_quant(hh,ni/2,nj/2,10)

  call dwt_qmerge(zw,ni,nj,ll,lh,hl,hh,ni/2,nj/2)
!   zd = 99.0
  call dwt_unshuffle(zw,zd,ni,nj)     ! unshuffle zw into zd
!   do j = nj/2,nj-1
!   do i = ni/2,ni-1
!     zd(i,j) = 0.0
!   enddo
!   enddo
  call dwt_shuffle(zw,zd,ni,nj)        ! shuffle zd into zw
!  z = 88.0
!  call dwt_shuffle(z,zd,ni,nj)        ! shuffle zd into z
!  call dwt_verif(zw,z,ni,nj,'shuf')   ! compare z and zw

  call fstecr(zd,zd,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','FWH4','WAVELET','X',0,0,0,0,5,.false.)
  call fstecr(z,z,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','FWHU','WAVELET','X',0,0,0,0,5,.false.)
!  call dwt_inv_lift_haar_r(zw,ni,nj,.false.,.true.)

  call dwt_inv_lift_haar_r(zw,ni,nj,.true.,.true.)              ! inverse 2D transform of zw
  call dwt_diag(zw,ni,nj,16,'zw(post)-')
  call fstecr(zw,zw,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','INVH','WAVELET','X',0,0,0,0,5,.false.)

  call dwt_verif(zw,z0,ni,nj,'haar')
!  zd = zw - z0
!  print *,'haar min/max hor error',minval(zd),maxval(zd)
!  print *,'haar avg hor error',sum(zd)/ni/nj
!  print *,'haar avg abs hor error',sum(abs(zd))/ni/nj
!  z = abs(zd/z0)
!  print *,'haar min/max rel hor error',minval(z),maxval(z)
!  print *,'haar avg rel hor error',sum(z)/ni/nj
!  print *,'haar avg rel abs hor error',sum(abs(z))/ni/nj  

  goto 777
  zw = z
  call dwt_d4(z0,z,ni,nj,.true.,.false.)         ! horizontal forward transform
  call dwt_d4(z,zw,ni,nj,.false.,.true.)         ! forward vertical transform

  call idwt_d4(z,zd,ni,nj,.true.,.false.)        ! horizontal inverse transform

  call fstecr(z,z,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','HORI','WAVELET','X',0,0,0,0,5,.false.)
  call fstecr(zd,zd,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','BAKH','WAVELET','X',0,0,0,0,5,.false.)
  call fstecr(zd-z0,zd,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','DELH','WAVELET','X',0,0,0,0,5,.false.)

  zd = zd - z0
  print *,'min/max hor error',minval(zd),maxval(zd)
  print *,'avg hor error',sum(zd)/ni/nj
  print *,'avg abs hor error',sum(abs(zd))/ni/nj
  z = abs(zd/z0)
  print *,'min/max rel hor error',minval(z),maxval(z)
  print *,'avg rel hor error',sum(z)/ni/nj
  print *,'avg rel abs hor error',sum(abs(z))/ni/nj  

  call fstecr(zw,zw,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','FULL','WAVELET','X',0,0,0,0,5,.false.)

  zh2 = zw(0:ni/2-1,0:nj/2-1)
  call dwt_d4(zh2,zh1,ni/2,nj/2,.true.,.false.)         ! horizontal forward transform
  call dwt_d4(zh1,zh,ni/2,nj/2,.false.,.true.)         ! forward vertical transform
  zw(0:ni/2-1,0:nj/2-1) = zh
  call fstecr(zw,zw,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','FUL2','WAVELET','X',0,0,0,0,5,.false.)
  zw(0:ni/2-1,0:nj/2-1) = zh2

!  call dwt_trim(zh(0:ni/4-1,0:nj/4-1),ni/4,nj/4,.0001,.true.,zer)   ! LL trim values < .0001
!  call dwt_trim(zh(0:ni/4-1,nj/4:nj-1),ni/4,nj/4,.0002,.true.,zer)   ! LH trim values < .001
!  call dwt_trim(zh(ni/4:ni-1,0:nj/4-1),ni/4,nj/4,.0002,.true.,zer)   ! HL trim values < .001
  call dwt_trim(zh(ni/4:ni-1,nj/4:nj-1),ni/4,nj/4,.0002,.false.,zer)  ! HH trim values < .001
  call idwt_d4(zh,zh1,ni/2,nj/2,.false.,.true.)         ! inverse vertical transform
  call idwt_d4(zh1,zh,ni/2,nj/2,.true.,.false.)         ! horizontal inverse transform
  
  zw(0:ni/2-1,0:nj/2-1) = zh

  zh = zw(0:ni/2-1,0:nj/2-1)
  call fstecr(zh,zh,-32,iun,0,0,0,ni/2,nj/2,1,0,0,0, &
              'XX','LL00','WAVELET','X',0,0,0,0,5,.false.)

  zh = zw(0:ni/2-1,nj/2:nj-1)
  call fstecr(zh,zh,-32,iun,0,0,0,ni/2,nj/2,1,0,0,0, &
              'XX','LH00','WAVELET','X',0,0,0,0,5,.false.)

  zh = zw(ni/2:ni-1,0:nj/2-1)
  call fstecr(zh,zh,-32,iun,0,0,0,ni/2,nj/2,1,0,0,0, &
              'XX','HL00','WAVELET','X',0,0,0,0,5,.false.)

  zh = zw(ni/2:ni-1,nj/2:nj-1)
  call fstecr(zh,zh,-32,iun,0,0,0,ni/2,nj/2,1,0,0,0, &
              'XX','HH00','WAVELET','X',0,0,0,0,5,.false.)

  
  z = 0
  call dwt_trim(zw(0:ni/2-1,0:nj/2-1),ni/2,nj/2,.0001,.false.,zer) ! LL trim values < .0001
  call dwt_trim(zw(0:ni/2-1,nj/2:nj-1),ni/2,nj/2,.0002,.false.,zer)  ! LH trim values < .001
  call dwt_trim(zw(ni/2:ni-1,0:nj/2-1),ni/2,nj/2,.0002,.false.,zer)  ! HL trim values < .001
  call dwt_trim(zw(ni/2:ni-1,nj/2:nj-1),ni/2,nj/2,.0005,.true.,zer) ! HH trim values < .001
  call idwt_d4(zw,zd,ni,nj,.false.,.true.)          ! inverse vertical transform
  call idwt_d4(zd,z,ni,nj,.true.,.false.)           ! inverse horizontal transform

  call fstecr(z,z,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','BACK','WAVELET','X',0,0,0,0,5,.false.)

  zw = z0 - z
  call fstecr(zw,zw,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','DELT','WAVELET','X',0,0,0,0,5,.false.)

  print *,'min/max error',minval(zw),maxval(zw)
  print *,'avg error',sum(zw)/ni/nj
  print *,'avg abs error',sum(abs(zw))/ni/nj
  z = abs(zw/z0)
  call fstecr(z,z,-32,iun,0,0,0,ni,nj,1,0,0,0, &
              'XX','RERR','WAVELET','X',0,0,0,0,5,.false.)
  print *,'min/max rel error',minval(z),maxval(z)
  print *,'avg rel error',sum(z)/ni/nj
  print *,'avg rel abs error',sum(abs(z))/ni/nj  
  nerr = 0
  do j = 0 , nj-1
  do i = 1 , ni-1
    if(abs(zw(i,j)) > 1.0) then
      nerr = nerr + 1
    endif
  enddo
  enddo
  print *, 'number of large errors =',nerr
777 continue
  call fstfrm(iun)
  call fclos(iun)
  return
end subroutine fst_wav_1
subroutine dwt_d4(zi,zo,ni,nj,alongx,alongy)  ! daubechies DB4 forward transform
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(0:ni-1, 0:nj-1), intent(IN)  :: zi
  real, dimension(0:ni-1, 0:nj-1), intent(OUT) :: zo
  logical, intent(IN) :: alongx, alongy

  real, dimension(0:ni-1) :: ti
  integer :: i, j, half, ii, jj
  real, parameter :: h0=(1.0+sqrt(3.0))/(4.0*sqrt(2.0))
  real, parameter :: h1=(3.0+sqrt(3.0))/(4.0*sqrt(2.0))
  real, parameter :: h2=(3.0-sqrt(3.0))/(4.0*sqrt(2.0))
  real, parameter :: h3=(1.0-sqrt(3.0))/(4.0*sqrt(2.0))
  real, parameter :: g0=h3
  real, parameter :: g1=-h2
  real, parameter :: g2=h1
  real, parameter :: g3=-h0

  if(alongx) then
    half = ni/2
    do j = 0, nj-1
      ii = 0
      do i = 0, ni-4, 2
        ti(ii)      = zi(i,j)*h0 + zi(i+1,j)*h1 + zi(i+2,j)*h2 + zi(i+3,j)*h3
        ti(ii+half) = zi(i,j)*g0 + zi(i+1,j)*g1 + zi(i+2,j)*g2 + zi(i+3,j)*g3
        ii = ii + 1
      enddo
      ti(half-1) = zi(ni-2,j)*h0 + zi(ni-1,j)*h1 + zi(0,j)*h2 + zi(1,j)*h3
      ti(ni-1)   = zi(ni-2,j)*g0 + zi(ni-1,j)*g1 + zi(0,j)*g2 + zi(1,j)*g3
      zo(0:ni-1,j) = ti(0:ni-1)
    enddo
  endif

  if(alongy) then
    half = nj/2
      jj = 0
      do j = 0, nj-4, 2
        zo(:,jj)      = zi(:,j)*h0 + zi(:,j+1)*h1 + zi(:,j+2)*h2 + zi(:,j+3)*h3
        zo(:,jj+half) = zi(:,j)*g0 + zi(:,j+1)*g1 + zi(:,j+2)*g2 + zi(:,j+3)*g3
        jj = jj + 1
      enddo   ! j
      zo(:,half-1) = zi(:,nj-2)*h0 + zi(:,nj-1)*h1 + zi(:,0)*h2 + zi(:,1)*h3
      zo(:,nj-1) = zi(:,nj-2)*g0 + zi(:,nj-1)*g1 + zi(:,0)*g2 + zi(:,1)*g3
  endif
  return
end subroutine dwt_d4

subroutine idwt_d4(zi,zo,ni,nj,alongx,alongy)  ! daubechies DB4 inverse transform
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(0:ni-1, 0:nj-1), intent(IN)  :: zi
  real, dimension(0:ni-1, 0:nj-1), intent(OUT) :: zo
  logical, intent(IN) :: alongx, alongy

  real, dimension(0:ni-1) :: ti
!  real, dimension(0:nj-1) :: tj
  integer, parameter :: strip = 16
  real, dimension(0:ni-1,0:nj-1) :: tj
  integer :: i, j, half, ii, jj, i0, in, sl

  real, parameter :: h0=(1.0+sqrt(3.0))/(4.0*sqrt(2.0))
  real, parameter :: h1=(3.0+sqrt(3.0))/(4.0*sqrt(2.0))
  real, parameter :: h2=(3.0-sqrt(3.0))/(4.0*sqrt(2.0))
  real, parameter :: h3=(1.0-sqrt(3.0))/(4.0*sqrt(2.0))
  real, parameter :: g0=h3
  real, parameter :: g1=-h2
  real, parameter :: g2=h1
  real, parameter :: g3=-h0
  real, parameter :: Ih0=h2
  real, parameter :: Ih1=g2
  real, parameter :: Ih2=h0
  real, parameter :: Ih3=g0
  real, parameter :: Ig0=h3
  real, parameter :: Ig1=g3
  real, parameter :: Ig2=h1
  real, parameter :: Ig3=g1

  if(alongx) then
    half = ni/2
    do j = 0, nj-1
      ti(0) = zi(half-1,j)*Ih0 + zi(ni-1,j)*Ih1 + zi(0,j)*Ih2 + zi(half,j)*Ih3
      ti(1) = zi(half-1,j)*Ig0 + zi(ni-1,j)*Ig1 + zi(0,j)*Ig2 + zi(half,j)*Ig3
      ii = 2
      do i = 0, half-2
        ti(ii) = zi(i,j)*Ih0 + zi(i+half,j)*Ih1 + zi(i+1,j)*Ih2 + zi(i+half+1,j)*Ih3
        ii = ii + 1
        ti(ii) = zi(i,j)*Ig0 + zi(i+half,j)*Ig1 + zi(i+1,j)*Ig2 + zi(i+half+1,j)*Ig3
        ii = ii + 1
      enddo
      zo(0:ni-1,j) = ti(0:ni-1)
    enddo
  endif

  if(alongy) then
    half = nj/2
    do i0 = 0, ni-1, strip
      in = min(i0+strip-1,ni-1)
      sl = in-i0
      jj = 2
      tj(0:sl,0) = zi(i0:in,half-1)*Ih0 + zi(i0:in,nj-1)*Ih1 + zi(i0:in,0)*Ih2 + zi(i0:in,half)*Ih3
      tj(0:sl,1) = zi(i0:in,half-1)*Ig0 + zi(i0:in,nj-1)*Ig1 + zi(i0:in,0)*Ig2 + zi(i0:in,half)*Ig3
      do j = 0, half-2
        tj(0:sl,jj) = zi(i0:in,j)*Ih0 + zi(i0:in,j+half)*Ih1 + zi(i0:in,j+1)*Ih2 + zi(i0:in,j+half+1)*Ih3
        jj = jj + 1
        tj(0:sl,jj) = zi(i0:in,j)*Ig0 + zi(i0:in,j+half)*Ig1 + zi(i0:in,j+1)*Ig2 + zi(i0:in,j+half+1)*Ig3
        jj = jj + 1
      enddo
      zo(i0:in,:) = tj(0:sl,:)
    enddo  ! i0
  endif
  return
end subroutine idwt_d4
subroutine dwt_trim(z,ni,nj,val,quant,zero)
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
  real, intent(IN) :: val
  logical, intent(IN) :: quant, zero

  integer :: i, j, count, temp
  integer, dimension(0:7) :: cnts

  count = 0
  cnts = 0
  do j = 1, nj-1
  do i = 1, ni-1
    temp = z(i,j)/val + .5
    if(abs(temp)<8) cnts(abs(temp)) = cnts(abs(temp)) +1
    if(quant) then
      z(i,j) = temp * val
    endif
    if(abs(z(i,j)) < val .and. zero) then
      z(i,j) = 0.0
      count = count + 1
    endif
  enddo
  enddo
  print *, 'trimmed',count,' values out of',ni*nj,cnts
  return
end subroutine dwt_trim

subroutine dwt_quant(z,ni,nj,nbits)  ! shuffle zu into zs
  implicit none
  integer, intent(IN) :: ni, nj, nbits
  real, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z

  integer :: i, j, intval
  real :: high, low, span
  real *8 :: quant, fact, mult

  high  = maxval(z)
  low   = minval(z)
  span  = high -low
  if(span == 0.0) return   ! no quantization needed
  mult  = 2.0**nbits
  quant = span / mult  ! quantum size
  fact  = mult / span  ! 1 / quantum
  do j=0,nj-1
  do i=0,ni-1
    intval = nint( (z(i,j)) * fact)  ! number of quanta in value
    z(i,j) = intval * quant          ! quantized value
  enddo
  enddo
  return
end subroutine dwt_quant

subroutine dwt_diag(z,ni,nj,nbits,text)  ! shuffle zu into zs
  implicit none
  integer, intent(IN) :: ni, nj, nbits
  real, dimension(0:ni-1, 0:nj-1), intent(INOUT) :: z
  character(len=*) :: text

  integer :: i, j, npts
  integer, dimension(0:ni-1, 0:nj-1) :: ival
  real :: high, low, span, entropy
  real *8 :: quant, fact, mult, sum1, sum2, avg, var, ovnpts, log2
  integer, dimension(0:65600) :: itab

  log2 = log(2.0)
  log2 = 1.0/log2
!  print *,'log2=',log2*log(2.0)
  npts  = ni*ni
  ovnpts = npts ; ovnpts = 1.0 / ovnpts
  high  = maxval(z)
  low   = minval(z)
  sum1  = sum(z)
  sum2  = sum(z*z)
  avg   = sum1 / npts
  var   = max(0.d0,1.0d0*(sum2 + avg*avg*npts - 2*avg*sum1) / npts)
  var   = sqrt(var)
  span  = high -low
  if(span == 0.0) return   ! no quantization needed
  mult  = 2.0**nbits
  quant = span / mult  ! quantum size
  fact  = mult / span  ! 1 / quantum
  ival = nint((z-low) * fact)
  entropy = nbits
  if(nbits <= 16) then
    itab = 0
    do j=0,nj-1
    do i = 0,ni-1
      itab(ival(i,j)) = itab(ival(i,j)) + 1
    enddo
    enddo
    entropy = 0
    do i=0,65535
      if(itab(i) .ne.0) entropy = entropy - (itab(i)*ovnpts) * log(itab(i)*ovnpts) * log2
    enddo
  endif
  print 101,trim(text)//'min/max/range/avg/var/quantum :',low,high,span,avg,var,quant,entropy,1.0*nbits
!   print 101,trim(text)//'quantum               :'
100 format(A40,10G12.4)
101 format(A40,10G15.8)
  return
end subroutine dwt_diag

subroutine dwt_verif(z1,z2,ni,nj,text)
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(0:ni-1, 0:nj-1), intent(IN) :: z1
  real, dimension(0:ni-1, 0:nj-1), intent(IN) :: z2
  character(len=*), intent(IN) :: text

  integer :: i, j, npts
  real *8 :: sum1, sum2, var1, var2, sum2z1, sum2z2, moyz1, moyz2, sumz1, sumz2, sumz1z2, var12
  real, dimension(0:ni-1, 0:nj-1) :: zd,z

  npts = ni*nj

  sumz1 = sum(z1)
  sumz2 = sum(z2)
  sumz1z2=sum(z1*z2)
  moyz1 = sumz1/npts
  moyz2 = sumz2/npts
  sum2z1 = sum(z1*z1)
  sum2z2 = sum(z2*z2)
  var1 = max(0.d0,1.0d0*(sum2z1 + moyz1*moyz1*npts - 2*moyz1*sumz1) / npts)
  var1 = sqrt(var1)
  var2 = max(0.d0,1.0d0*(sum2z2 + moyz2*moyz2*npts - 2*moyz2*sumz2) / npts)
  var2 = sqrt(var2)
  var12 = sumz1z2/sqrt(sum2z1*sum2z2)

  zd = z1 - z2
  sum1 = sum(zd)
  sum2 = sum(zd*zd)
  print 100,trim(text)//' min z1, z2   ',minval(z1),minval(z2)
  print 100,trim(text)//' max z1, z2   ',maxval(z1),maxval(z2)
  print 100,trim(text)//' avg z1, z2   ',moyz1,moyz2
  print 100,trim(text)//' var z1, z2   ',var1,var2
  print 101,trim(text)//' cor z1, z2   ',var12
  print 100,trim(text)//' min/max error',minval(zd),maxval(zd)
  print 100,trim(text)//' bias         ',sum1/ni/nj
  print 100,trim(text)//' avg abs error',sum(abs(zd))/ni/nj
  print 100,trim(text)//' avg rms error',sqrt(sum2/ni/nj)
  return
100 format(A,5G12.4)
101 format(A,5G15.8)
!  z = abs(zd/z0)
!  print *,'haar min/max rel hor error',minval(z),maxval(z)
!  print *,'haar avg rel hor error',sum(z)/ni/nj
!  print *,'haar avg rel abs hor error',sum(abs(z))/ni/nj  
end subroutine dwt_verif
