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
!  call fstopi("MSGLVL",8,0)
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
        call test_quantizing(z,ni,nj,nomvar)
        z = max(z,0.0)
  !      call test_compression(z,ni,nj,nomvar)
        call test_quantizing(z,ni,nj,nomvar)
      endif
      if(nomvar .ne. oldnam) then
        if(oldnam .ne. '    ') write(0,*)'NAME=',oldnam,' levels=',ilev
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

subroutine test_quantizing(z,ni,nj,nomvar)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni,nj
  real, dimension(ni,nj), intent(IN) :: z
  character(len=4), intent(IN) :: nomvar

  real *8 :: bias, rms
  real :: errmax, znew, zmin, zmax, err, rrange, toler, oldzmin, r, rr, meanval, span
  integer :: i, j, coded, nbits, power, imax, izero
  integer, dimension(0:16) :: population
  integer :: ipop

  zmin = z(1,1)
  zmax = z(1,1)
  do j=1,nj
  do i=1,ni
    zmin = min(z(i,j),zmin)
    zmax = max(z(i,j),zmax)
  enddo
  enddo
  rrange = zmax - zmin
  oldzmin = zmin
  if(rrange == 0) return
  do nbits = 8, 16, 4
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
    meanval = 0
    imax = nint((zmax-zmin)*r)
    population = 0
    do j=1,nj
    do i=1,ni
      meanval = meanval + z(i,j)
      coded = nint((z(i,j)-zmin)*r)
      ipop = ishft(coded,3-nbits)
      population(ipop) = population(ipop)+1
      if(coded==0) population(8) = population(8)+1
      znew = coded*rr + zmin
      err = znew - z(i,j)
      errmax = max(abs(err),errmax)
      bias = bias + err
      rms = rms + err * err
    enddo
    enddo
    bias =  bias / (ni*nj)
    rms = sqrt(rms / (ni*nj))
    meanval = meanval / (ni*nj)  ! moyenne
    span = max(min(zmax-meanval,meanval-zmin),rr)
    span = rrange/span
!     return
    print 100, nbits, toler*1.2 >= errmax, toler, errmax, bias, rms, rrange, zmax, zmin, meanval, &
       nint(errmax/toler*100)*1.0, nint(rms/toler*100)*1.0, nint(bias/toler*1000000)*.0001, &
       (power-imax),nint(span)
    print 101,population(8),population(0)-population(8),population(1:7)
100 format(I4,L5, 2G12.4,1H|, 3G12.4,1H|, 3G12.4, 3F8.2, 2i6)
101 format(9I8)
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
