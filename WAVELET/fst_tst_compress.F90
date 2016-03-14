program test_compress
  use ISO_C_BINDING
  implicit none
  integer, external :: fnom, fstouv, fstnbr, fstinf, fstsui
  integer :: iun, status, nrec, key, ni, nj, nk, irec, ilev
  integer :: date,deet,npas,nbits,datyp,ip1,ip2,ip3,ig1,ig2,ig3,ig4
  integer :: swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=2) :: typvar
  character(len=4) :: nomvar, oldnam
  character(len=12) :: etiket
  character(len=1) :: grtyp
  real, dimension(:,:), pointer :: z=>NULL()
  integer, dimension(:,:), pointer :: iz=>NULL()

  write(0,*)'======= compression algorithm test ======='
  iun=0
  status = fnom(iun,'input.fst','RND+STD+R/O+OLD',0)
  call fstopi("MSGLVL",8,0)
  if(status < 0) goto 999
  status = fstouv(iun,'RND')
  if(status < 0) goto 999
  nrec = fstnbr(iun)
  irec = 0
  ilev = 0
  oldnam='    '
  write(0,*)nrec,' records found, unit=',iun

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
      if(ni > 1 .and. nj > 1 .and. trim(nomvar) .ne. 'KTkt') then
        call fstluk(z,key,ni,nj,nk)
        call test_quantizing(z,ni,nj,nomvar,0)
!         z = max(z,0.0)
  !      call test_compression(z,ni,nj,nomvar)
        call test_quantizing(z,ni,nj,nomvar,1)
      endif
      if(nomvar .ne. oldnam) then
!         if(oldnam .ne. '    ') write(0,*)'NAME=',oldnam,' levels=',ilev
        oldnam = nomvar
        ilev = 0
      endif
    endif
    key = fstsui(iun,ni,nj,nk)
  enddo
  write(0,*)'number of records processed:',irec,' out of',nrec
  call fstfrm(iun)

  stop
999 continue
  write(0,*)'=== ERROR opening files ==='
  stop
end program test_compress

#define NSCORES 16

subroutine scores(s,fa,fb,ni,nj,toler,diag,msg)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(NSCORES), intent(OUT) :: s
  real, dimension(ni,nj) :: fa, fb
  real, intent(IN) :: toler             ! tolerance
  logical, intent(IN) :: diag
  character(len=*), intent(IN) :: msg

  real *8 :: suma, sumb, sumab, suma2, sumb2, bias, rms, abserr, errmax, errrel
  real *8 :: errrelmax, errrelavg, avga, avgb, stda, stdb, sd, corr, slope, intercept
  integer :: i, j, n
  real :: nonzero, fmin, fmax, span

  s = 0.0
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
  do j = 1 , nj
  do i = 1 , ni
    suma = suma + fa(i,j)
    suma2 = suma2 + fa(i,j)*fa(i,j)
    sumb = sumb + fb(i,j)
    sumb2 = sumb2 + fb(i,j)*fb(i,j)
    sumab = sumab + fa(i,j)*fb(i,j)
    rms = rms + (fa(i,j)-fb(i,j)) * (fa(i,j)-fb(i,j))
    abserr = abserr + abs(fa(i,j)-fb(i,j))
    bias = bias +  (fb(i,j)-fa(i,j))
    errmax = max( errmax,  abs(fa(i,j)-fb(i,j)) )
    if( fa(i,j)*fb(i,j) .ne. 0) then
      nonzero = nonzero + 1
      errrel = abs(fa(i,j)-fb(i,j)) / abs(fa(i,j))
      errrel = max(errrel, .0000005)
      sd = sd - log10(errrel)
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
  bias = bias / (ni*nj)
  abserr = abserr / (ni*nj)
  if(nonzero > 0) then
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
             ' A,B=',slope-1.0,intercept
100 format(A,A,E8.2,A,E10.3,A,E9.3,A,I3,1H%,A,2E10.3,A,2E10.3,A,I3,1H%,A,E10.3,A,F4.1,A,E10.4,A,I5,A,E10.3,E11.3)
  endif

  return
end subroutine scores

subroutine test_quantizing(z,ni,nj,nomvar,mode)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni,nj
  real, dimension(ni,nj), intent(IN) :: z
  character(len=4), intent(IN) :: nomvar
  integer, intent(IN) :: mode

  real *8 :: bias, rms, relrms
  real :: errmax, znew, zmin, zmax, err, rrange, toler, oldzmin, r, rr, meanval, span, errrel, factor, zmin0, s1
  integer :: i, j, coded, nbits, power, imax, izero, nonzero, ifactor, iszero
  integer, dimension(0:16) :: population
  integer :: ipop
  integer, dimension(:,:), pointer :: iz
  type(C_PTR) :: piz
  integer :: mask
  real, dimension(ni,nj), target :: zz
  integer :: izz, maxexp1, maxexp2, maxexp
  real, dimension(NSCORES) :: s

!   if(mode .ne. 0) return
  zmin = z(1,1)
  zmin0 = 1.0
  zmax = z(1,1)
  meanval = 0.0
  do j=1,nj
  do i=1,ni
    zmin = min(z(i,j),zmin)
    zmax = max(z(i,j),zmax)
    if(z(i,j) .ne.0) zmin0 = min(z(i,j),zmin0)
    meanval = meanval + z(i,j)
  enddo
  enddo
  maxexp1 = transfer(zmin,maxexp1)
  maxexp1 = iand(ishft(maxexp1,-23),Z'000000FF')
  maxexp2 = transfer(zmax,maxexp2)
  maxexp2 = iand(ishft(maxexp2,-23),Z'000000FF')
  maxexp = max(maxexp1,maxexp2)
!   ifactor = 127+15-maxexp
!   ifactor = ishft(ifactor,23)
!   factor = transfer(ifactor,factor)
!   print 777," maxexp, factor, max, min, max', min' =",maxexp,factor,zmax,zmin, factor*zmax,factor*zmin
777 format(A,i6,6Z10.8)
  rrange = zmax - zmin
  if(rrange == 0) return
!   if(zmin == zmin0) return
  meanval = meanval / (ni*nj)  ! moyenne
  span = min(zmax-meanval,meanval-zmin)
  if(span == 0.0) span = rrange * .0001
  if ( rrange / span < 33 ) return
!   if ( rrange / span > 64) then
!     print *,nomvar,nint(rrange / span),zmin,zmin0,zmax,meanval
!     return
!   else
!     return
!   endif
!
  zz = z
  piz = C_LOC(zz(1,1))
  call c_f_pointer(piz,iz,[ni,nj])
  mask = ishft(-1,12)
!   mask = not(mask)
!   factor = 1.00001
!   zz = zz * factor
!   iz = iand(iz,mask)
!   mask = -1
  oldzmin = zmin
  if(rrange == 0) return
  do nbits = 16, 16, 4
    zmin = oldzmin
    power = ishft(1,nbits)
    rr = rrange/(power - 2.01)
    toler = rr * .5
    if(abs(zmin) < toler) zmin = 0
    r = 1.0 / rr
!     if(zmin > 0) zmin = zmin - mod(zmin,rr)
!     if(zmin < 0) zmin = zmin + mod(zmin,rr) - rr
   izero = -1
   if(zmin<0 .and. zmax > 0) then  ! make sure that 0 decodes back as zero
      zmin = -rr*nint((0.0-zmin)*r)  ! make zmin a multiple of rr
      izero = nint((0.0-zmin)*r)
      imax = nint((zmax-zmin)*r)
   endif
    errmax = 0.0
    bias = 0
    rms = 0
!     meanval = 0
    imax = nint((zmax-zmin)*r)
    population = 0
    relrms = 0
    nonzero = 0
    iszero = 0
    do j=1,nj
    do i=1,ni
!       meanval = meanval + z(i,j)
      coded = nint((z(i,j)-zmin)*r)
      ipop = ishft(coded,3-nbits)
      population(ipop) = population(ipop)+1
      if(coded==0) population(8) = population(8)+1
      if(mode == 0) znew = coded*rr + zmin
      if(mode .ne. 0) then
        zz(i,j) = z(i,j) * 1.000177
        izz = transfer(zz(i,j),izz)
        izz = iand(izz,mask)
        ifactor = iand(ishft(izz,-23),Z'000000FF')
        if(ifactor < maxexp-31) then
          izz = 0
          iszero = iszero + 1
        endif
        znew = transfer(izz,znew)
      endif
      zz(i,j) = znew
      err = znew - z(i,j)
      errmax = max(abs(err),errmax)
      bias = bias + err
      rms = rms + err * err
      errrel = 0
      if(z(i,j) .ne. 0 ) then
        nonzero = nonzero + 1
        errrel = ( abs(z(i,j)-znew) / abs(z(i,j)) )
      endif
      relrms = relrms + errrel
    enddo
    enddo
!     print *,'iszero =',iszero,ni*nj
    relrms = relrms / (nonzero)
    bias =  bias / (ni*nj)
    rms = sqrt(rms / (ni*nj))
!     meanval = meanval / (ni*nj)  ! moyenne
    span = max(min(zmax-meanval,meanval-zmin),rr)
    span = rrange/span
!     return
!     call s1scor(s1,z,zz,ni,nj,1,1,ni,nj,1)
    if(mode == 0)print 101,nomvar,population(8),population(0)-population(8),population(1:7)
    if(mode .ne. 0 .and. nint(span) < 32) return
    call scores(s,z,zz,ni,nj,toler,.true.,'    ')
!     print 100, nbits, toler*1.2 >= errmax, toler >= rms, toler,errmax,  bias,rms,relrms,rrange,s1,   zmax,zmin,meanval, &
!        nint(errmax/toler*100)*1.0, nint(rms/toler*100)*1.0, nint(bias/toler*1000000)*.0001, abs(bias/rms*100.0), &
!        (power-imax),nint(span),'  |'
100 format(I4,2X,2L1, 2G12.4,2H |, 5G12.4,2H |, 3G12.4,2H | ,4F9.2,2H | 2I6,A)
101 format(A,2X,9I8)
    if(imax > power - 2) then
      print *,'imax > power ',imax,izero,power,zmax,zmin,zmin-oldzmin,rr,(zmax-zmin)*r,(zmax-oldzmin)*r
      print *,zmin+(power-1)*rr,zmax,zmin+(power-1)*rr-zmax
      print *,oldzmin+(power-1)*rr,zmax,oldzmin+(power-1)*rr-zmax
      stop
    endif
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
