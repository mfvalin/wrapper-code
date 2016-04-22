module globalstats
  implicit none
  integer, dimension(-1:2), save :: biases, nbiases
end module globalstats
program test_compress
  use ISO_C_BINDING
  use globalstats
  implicit none
  integer, external :: fnom, fstouv, fstnbr, fstinf, fstsui
  integer :: iun, status, nrec, key, ni, nj, nk, irec, ilev, ilen
  integer :: date,deet,npas,nbits,datyp,ip1,ip2,ip3,ig1,ig2,ig3,ig4
  integer :: swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=2) :: typvar
  character(len=4) :: nomvar, oldnam
  character(len=12) :: etiket
  character(len=1) :: grtyp
  real, dimension(:,:), pointer :: z=>NULL()
  integer, dimension(:,:), pointer :: iz=>NULL()
  character (len=128) :: filename
  integer :: i, j

  write(0,*)'======= compression algorithm test ======='
  iun=0
  call get_command_argument(1,filename,ilen,status)
  if(status .ne. 0) stop
  status = fnom(iun,trim(filename),'RND+STD+R/O+OLD',0)
!   call fstopi("MSGLVL",8,0)
  if(status < 0) goto 999
  status = fstouv(iun,'RND')
  if(status < 0) goto 999
  nrec = fstnbr(iun)
  irec = 0
  ilev = 0
  oldnam='    '
  write(0,*)nrec,' records found, unit=',iun
  biases = 0
  nbiases = 0

  key = fstinf(iun,ni,nj,nk,-1,'            ',-1,-1,-1,'  ','    ')
  do while(key >= 0)
    irec = irec + 1
    if(ni>1 .and. nj>1) then
      if(associated(z)) then
        deallocate(z)
      endif
      if(associated(iz)) then
        deallocate(iz)
      endif
      allocate(z(ni,nj))
      allocate(iz(ni,nj))
      call fstprm(key,date,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3,  &
                  typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                  swa,lng,dltf,ubc,extra1,extra2,extra3)
      ilev = ilev + 1
!      if(trim(nomvar) == 'TT') then
!       if(ni > 5 .and. nj > 5 .and. trim(nomvar) .eq. 'TT' .and. ip1 <= 1000) then
      if(ni > 5 .and. nj > 5 .and. ip1 <= 850 .and. ip1 > 0) then
!       if(ni > 5 .and. nj > 5 ) then
!         write(filename,100)nomvar,ip2,'.data'
100     format(A2,I4.4,A5)
!         print *,"'"//trim(filename)//"'"
        call fstluk(z,key,ni,nj,nk)
!         open(33,file=trim(filename),form='FORMATTED')
!         write(33,*)ni,nj
!         write(33,*)((z(i,j),i=1,ni),j=1,nj)
!         close(33)
        call test_quantizing(z,ni,nj,nomvar,0)
        call test_quantizing(z,ni,nj,nomvar,-1)
!         z = max(z,0.0)
  !      call test_compression(z,ni,nj,nomvar)
        call test_quantizing(z,ni,nj,nomvar,1)
        call test_quantizing(z,ni,nj,nomvar,2)
      endif
      if(nomvar .ne. oldnam) then
!         if(oldnam .ne. '    ') write(0,*)'NAME=',oldnam,' levels=',ilev
        oldnam = nomvar
        ilev = 0
      endif
    endif
    key = fstsui(iun,ni,nj,nk)
  enddo
  write(6,*)'number of records processed:',irec,' out of',nrec
  write(6,*)'biases=',biases
  write(6,*)'nbiases=',nbiases
  call fstfrm(iun)

  stop
999 continue
  write(0,*)'=== ERROR opening files ==='
  stop
end program test_compress

#define NSCORES 16

subroutine scores(s,fa,fb,ni,nj,toler,diag,msg,mode)
  use ISO_C_BINDING
  use globalstats
  implicit none
  integer, intent(IN) :: ni, nj, mode
  real, dimension(NSCORES), intent(OUT) :: s
  real, dimension(ni,nj) :: fa, fb
  real, intent(IN) :: toler             ! tolerance
  logical, intent(IN) :: diag
  character(len=*), intent(IN) :: msg

  real *8 :: suma, sumb, sumab, suma2, sumb2, bias, rms, abserr, errmax, errrel, one
  real *8 :: errrelmax, errrelavg, avga, avgb, stda, stdb, corr, slope, intercept
  integer :: i, j, n, isd
  real :: nonzero, fmin, fmax, span, sdl, sd
  integer, dimension(6) :: dig

  s = 0.0
  one = 1.0
  call s1scor(s(1),fa,fb,ni,nj,1,1,ni,nj,1)
  suma = 0.0
  sumb = 0.0
  sumab = 0.0
  suma2 = 0.0
  sumb2 = 0.0
  bias = 0.0
  abserr = 0.0
  rms = 0.0
  errmax = 0.0
  errrel = 0.0
  errrelavg = 0.0
  nonzero = 0
  sd = 0       ! significant digits
  n = ni*nj
  corr = 0.0
  slope = 0.0
  intercept = 0.0
  fmin = fa(1,1)
  fmax = fmin
  dig = 0
  do j = 1 , nj
  do i = 1 , ni
    suma = suma + fa(i,j)
    suma2 = suma2 + one*fa(i,j)*fa(i,j)
    sumb = sumb + fb(i,j)
    sumb2 = sumb2 + one*fb(i,j)*fb(i,j)
    sumab = sumab + one*fa(i,j)*fb(i,j)
    rms = rms + (fa(i,j)-fb(i,j)) * (fa(i,j)-fb(i,j))
    abserr = abserr + abs(fa(i,j)-fb(i,j))
    bias = bias +  (fb(i,j)-fa(i,j))
    errmax = max( errmax,  abs(fa(i,j)-fb(i,j)) )
!     if( fa(i,j)*fb(i,j) .ne. 0) then
    if( fa(i,j) .ne. 0) then
      nonzero = nonzero + 1
      errrel = abs(fa(i,j)-fb(i,j)) / abs(fa(i,j))
      errrel = max(errrel, .0000005)
      sdl = log10(errrel)
      isd = max(1, min(6,nint(-sdl)))
      dig(isd) = dig(isd)+1
      sd = sd - sdl
      errrelavg = errrelavg + errrel
      errrelmax = max( errrelmax , errrel )
    endif
    fmin = min(fmin,fa(i,j))
    fmax = max(fmax,fa(i,j))
  enddo
  enddo
  avga = suma / (ni*nj)
  stda = sqrt( max( 0.0 , (suma2/n) - ((suma/n)*(suma/n)) ) )
  avgb = sumb / (ni*nj)
  stdb = sqrt( max( 0.0 , (sumb2/n) - ((sumb/n)*(sumb/n)) ) )
  rms = sqrt( rms / (ni*nj) )
  if(rms <= 1.0E-20) rms = 1.0E-20
  bias = bias / (ni*nj)
!   if(abs(bias) <= 1.0E-20) bias = 1.0E-20
!   rms = max(abs(bias),rms)
  abserr = abserr / (ni*nj)
  if(nonzero > 0) then
    dig = nint(dig * 100.0 / nonzero)
    sd = sd / nonzero
    errrelavg = errrelavg / nonzero
  endif
  span = min( abs(fmax-avga) , abs(avga-fmin) )
  span = (fmax - fmin) / span
  corr = (sumab - n*avga*avgb) / (n*stda*stdb)
  slope = (sumab/n - avga*avgb) / ( suma2/n - avga*avga)
  intercept = avgb - slope * avga
  s(2) = rms           ! average squared error
  s(3) = bias          ! bias
  s(4) = abserr        ! average absolute error
  s(5) = errmax        ! max absolute error
  s(6) = errrelmax     ! max realtive error
  s(7) = errrelavg     ! average relative error
  s(8) = stda
  s(9) = stdb
  s(10) = sd
  s(11) = corr
  if(diag) then
    print 100,msg,' S1=',s(1),' Bias=',s(3),' RMS=',s(2),' B/R=',nint(100*abs(s(3)/s(2))), &
             ' Ae=',s(5),s(4),' Re=',s(6),s(7),' NZ=',nint(nonzero/(ni*nj)*100.0), &
             ' dSTD=',stda-stdb,' SD=',sd,' dC1=',1.0-corr,' Sp=',nint(span*.5), &
             ' A,B=',1.0-slope,intercept,dig
100 format(A,A,E8.2,A,E9.2,A,E8.2,A,I2,1H%,A,E8.2,E9.2,A,E8.2,E9.2,A,I3,1H%,A,E9.2,A,F3.1,A,E8.2,A,I5,A,E9.2,E10.2,6I3)
  endif

  nbiases(mode) = nbiases(mode) + 1
  if(bias > 0) then
    biases(mode) = biases(mode) + 1
  else
    biases(mode) = biases(mode) - 1
  endif
  return
end subroutine scores

subroutine un_quantize(z,iz,ni,nj,imode)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni,nj
  real, dimension(ni,nj), intent(OUT) :: z
  integer, dimension(ni,nj), intent(IN) :: iz
  integer, intent(IN) :: imode
  interface
    subroutine float_unpacker(field,header,data,nelm,nbits) bind(C,name='c_float_unpacker')
      import C_FLOAT, C_INT
      real(C_FLOAT), dimension(nelm), intent(OUT) :: field
      integer(C_INT), dimension(*), intent(IN) :: header,data
      integer(C_INT), value :: nelm
      integer(C_INT), intent(OUT) :: nbits
    end subroutine float_unpacker
  end interface

  integer :: mode, i, j, nbits
  real :: zz

  mode = imode
  if(mode <= 0) then        ! linear quantization
    if(mode == -1) then
      call float_unpacker(z,iz,iz(11,1),ni*nj,nbits)
      return
    endif
    do j=1,nj
    do i=1,ni
      zz = transfer(iz(i,j),zz)
      z(i,j) = zz
    enddo
    enddo
  else
    if(mode == 1) then      ! 5 bits exponent, 11 bits mantissa
      call pseudo_ieee_un_quantize
    else                    ! 4 bits exponent, 12 bits mantissa
      call pseudo_ieee_un_quantize
    endif
  endif
contains
  subroutine pseudo_ieee_un_quantize
    implicit none
    integer :: i, j
    real :: zz

    do j=1,nj
    do i=1,ni
      zz = transfer(iz(i,j),zz)
      z(i,j) = zz
    enddo
    enddo
  end subroutine pseudo_ieee_un_quantize
end subroutine un_quantize

subroutine quantize(z,iz,ni,nj,span,rrange,imode)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni,nj
  real, dimension(ni,nj), intent(IN) :: z
  integer, dimension(ni,nj), intent(OUT) :: iz
  real, intent(OUT) :: span, rrange
  integer, intent(IN) :: imode
  interface
    subroutine float_packer(field,nbits,header,data,nelm) bind(C,name='c_float_packer')
      import C_FLOAT, C_INT
      real(C_FLOAT), dimension(nelm), intent(IN) :: field
      integer(C_INT), dimension(*), intent(OUT) :: header,data
      integer(C_INT), value :: nelm
      integer(C_INT), value :: nbits
    end subroutine float_packer
  end interface

  integer :: i, j, mode, mask, lowexp, maxexp, iszero, power, nonzero
  real :: zmin, zmin0, zmax, meanval, fudge, rr, toler, r, zz

  mode = imode
  nonzero = 0
  call extrema
  if(nonzero < ni*nj/8 .and. mode > 0) then
    span = 0.0
    return
  endif
  if(mode <= 0) then        ! linear quantization
    if(rrange == 0.0) return
!     if(span < 16.0) return
    if(mode == -1) then
      call float_packer(z,16,iz,iz(11,1),ni*nj)
      return
    endif
    power = ishft(1,16)
    rr = rrange/(power - 2.01)
    toler = rr * .5
    if(abs(zmin) < toler) zmin = 0
    r = 1.0 / rr
    do j=1,nj
    do i=1,ni
      iz(i,j) = nint((z(i,j)-zmin)*r)
      zz = iz(i,j)*rr + zmin
      iz(i,j) = transfer(zz,iz(i,j))
    enddo
    enddo
  else
    if(rrange == 0.0) return
!     if(span < 16.0) return
    if( span > 96 .and. mode == 2) mode = 1
    if(mode == 1) then      ! 5 bits exponent, 11 bits mantissa
      mask = ishft(-1,12)
      lowexp = 31
      fudge = 1.0001762
    else                    ! 4 bits exponent, 12 bits mantissa
      mask = ishft(-1,11)
      lowexp = 15
      fudge = 1.000099
    endif
    call pseudo_ieee_quantize
  endif
  return
contains
  subroutine extrema
    implicit none
    integer :: i, j, maxexp1, maxexp2
    zmin = z(1,1)
    zmin0 = 1.0
    zmax = z(1,1)
    meanval = 0.0
    do j=1,nj
    do i=1,ni
      zmin = min(z(i,j),zmin)
      zmax = max(z(i,j),zmax)
      if(z(i,j) .ne. 0.0) then
        nonzero = nonzero + 1
        zmin0 = min(z(i,j),zmin0)
      endif
      meanval = meanval + z(i,j)
    enddo
    enddo
    rrange = zmax - zmin
    if(rrange == 0.0) return
    meanval = meanval / (ni*nj)  ! moyenne
    span = min(zmax-meanval,meanval-zmin)
    if(span == 0.0) span = rrange * .0001
    span = 0.5 * (rrange / span)
    maxexp1 = transfer(zmin,maxexp1)
    maxexp1 = iand(ishft(maxexp1,-23),Z'000000FF')
    maxexp2 = transfer(zmax,maxexp2)
    maxexp2 = iand(ishft(maxexp2,-23),Z'000000FF')
    maxexp = max(maxexp1,maxexp2)
  end subroutine extrema

  subroutine pseudo_ieee_quantize
    implicit none
    integer :: i, j, izz, ifactor
    real :: zz

    iszero = 0
    do j=1,nj
    do i=1,ni
      zz = z(i,j) * fudge
      izz = transfer(zz,izz)
      izz = iand(izz,mask)
      ifactor = iand(ishft(izz,-23),Z'000000FF')
      if(ifactor < maxexp-lowexp) then
        izz = 0
        iszero = iszero + 1
      endif
      iz(i,j) = izz
    enddo
    enddo
  end subroutine pseudo_ieee_quantize
end subroutine quantize

subroutine test_quantizing(z,ni,nj,nomvar,imode)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni,nj
  real, dimension(ni,nj), intent(IN) :: z
  character(len=4), intent(IN) :: nomvar
  integer, intent(IN) :: imode

  integer :: nbits
  integer, dimension(ni,nj) :: iz
  real, dimension(ni,nj) :: zz
  integer :: mode
  real, dimension(NSCORES) :: s
  real :: span,rrange, toler

!   if(mode .ne. 0) return
  mode = imode
  do nbits = 16, 16, 4
    call quantize(z,iz,ni,nj,span,rrange,mode)
    if(rrange == 0.0) return
!     if(span < 16.0) return
    if(mode > 0 .and. nint(span) < 16) return
    call un_quantize(zz,iz,ni,nj,mode)
!     if(mode == 0)print *,nomvar
    call scores(s,z,zz,ni,nj,toler,.true.,'    ',mode)
101 format(A,2X,9I8)
  enddo
  return
end subroutine test_quantizing

! subroutine test_compression(z,ni,nj,nomvar)
!   use ISO_C_BINDING
!   implicit none
!   integer, intent(IN) :: ni,nj
!   real, dimension(ni,nj), intent(IN) :: z
!   character(len=4), intent(IN) :: nomvar
!   interface
! !   int TO_jpeg2000(unsigned char *cin,int width,int height,int nbits,
! !                  int ltype, int ratio, int retry, char *outjpc, 
! !                  int jpclen)
!     function enc_jpeg(cin,width,height,nbits,ltype,ratio,retry,cout,jpclen) result(status) bind(C,name='TO_jpeg2000') 
!       import
!       implicit none
!       integer(C_INT), intent(IN), value :: width,height,nbits,ltype,retry,jpclen
!       real(C_FLOAT), intent(IN), value :: ratio
!       type(C_PTR), intent(IN), value :: cin, cout
!       integer(C_INT) :: status
!     end function enc_jpeg
! !   int FROM_jpeg2000(char *injpc,int bufsize,int *outfld)
!     function dec_jpeg(cin,bufsize,cout) result(status) bind(C,name='FROM_jpeg2000') 
!       import
!       implicit none
!       integer(C_INT), intent(IN), value :: bufsize
!       type(C_PTR), intent(IN), value :: cin, cout
!       integer(C_INT) :: status
!     end function dec_jpeg
! 
!   end interface
! 
!   integer, dimension(ni,nj) :: iz
!   real :: the_min, the_max, error
!   integer :: nbits, nbts, minv, nbits2, mn1, mx1, nbytes, njpeg
!   integer, external :: dwt_quantize, dwt_lorenzo, dwt_nbits, dwt_pack
!   integer(C_CHAR), dimension(10000000), target :: buf
!   integer(C_CHAR), dimension(10000000), target :: jbuf
! 
!   error = 0.0
!   nbts = 16
!   nbits = dwt_quantize(z,iz,ni*nj,nbts,error,the_min,the_max)
!   mn1 = minval(iz)
!   mx1 = maxval(iz)
!   minv = dwt_lorenzo(iz,ni,nj,.true.)
!   iz = iz - minv
!   nbits2=dwt_nbits(maxval(iz))
!   nbytes = dwt_pack(iz,ni*nj,buf,10000000,nbits)
!   njpeg = enc_jpeg(c_loc(buf),ni,nj,nbits2,0,1.0,0,c_loc(jbuf),10000000)
!   
!   write(0,100)nomvar,nbits,the_min,the_max,mn1,mx1,minval(iz),maxval(iz),nbits2,minv,nbytes,njpeg,ni*nj*4
! 100 format(1X,A4,1X,I3,2G12.5,10I8)
!   return
! end subroutine test_compression
