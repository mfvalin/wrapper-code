module globalstats
  implicit none
  integer, dimension(-1:4), save :: biases, nbiases
end module globalstats

#define NSCORES 16

subroutine bilorentz(p, q, ni, nj) ! 2 level Lorenzo predictor
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT) :: q
  integer :: i, j
  real :: e
  do j = nj, 3, -1
    do i = ni, 3, -1
      e = 2.0*p(i-1,j) - p(i-2,j) + 2.0*p(i,j-1) -4.0*p(i-1,j-1) + 2.0*p(i-2,j-1) - p(i,j-2) + 2.0*p(i-1,j-2) - p(i-2,j-2)
      q(i,j) = p(i,j) - e
    enddo
    q(2,j) = p(2,j) - (p(1,j) + p(2,j-1) - p(1,j-1))  ! Lorenzo prediction
    q(1,j) = q(1,j) - q(1,j-1)                        ! delta
  enddo
  do i = ni, 2, -1                                    ! row 2
    q(i,2) = p(i,2) - (p(i-1,2) + p(i,1) - p(i-1,1))
  enddo
  q(1,2) = p(1,2) - p(1,1)
  do i = ni, 2, -1                                    ! row 1
    q(i,1) = p(i,1) - p(i-1,1)
  enddo
  q(1,1) = p(1,1)
  q(1,1) = 0
end

subroutine lorenzo(p, q, ni, nj) ! Lorenzo predictor
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT) :: q
  integer :: i, j
  do j = nj, 2, -1
    do i = ni, 2, -1
      q(i,j) = p(i,j) - (p(i-1,j) + p(i,j-1) - p(i-1,j-1))
    enddo
    q(1,j) = p(1,j) - p(1,j-1)
  enddo
  do i = ni, 2, -1                   ! row 1
    q(i,1) = p(i,1) - p(i-1,1)
  enddo
  q(1,1) = p(1,1)
  q(1,1) = 0
end

subroutine my_slope32(q, ni, nj)  ! replace values with 32x32 local average + x-y slope correction
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(INOUT) :: q
  real :: avg, xlo, xhi, ylo, yhi, dx, dy
  integer :: i, j, ii, jj
  do j = 1, nj-31, 32
  do i = 1, ni-31, 32
    avg = sum(q(i:i+31,j:j+31)) / 1024.0
!     xlo = sum(q(i   :i+15,j+8   :j+23  )/256.0)
!     xhi = sum(q(i+16:i+31,j+8   :j+23  )/256.0)
    xlo = sum(q(i   :i+15,j   :j+31  )/512.0)
    xhi = sum(q(i+16:i+31,j   :j+31  )/512.0)
    dx = (xhi-xlo) / 16.0
!     ylo = sum(q(i+8 :i+23  ,j   :j+15)/256.0)
!     yhi = sum(q(i+8 :i+23  ,j+16:j+31)/256.0)
    ylo = sum(q(i   :i+31  ,j   :j+15)/512.0)
    yhi = sum(q(i   :i+31  ,j+16:j+31)/512.0)
    dy = (yhi-ylo) / 16.0
    do jj = 0, 31
    do ii = 0, 31
      q(i+ii,j+jj) = avg + dx * (ii-15.5) + dy * (jj-15.5)
    enddo
    enddo
!     q(i:i+31,j:j+31) = avg
  enddo
  enddo
end

subroutine my_filter32(q, ni, nj)  ! replace values with 32x32 local average
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(INOUT) :: q
  real :: avg
  integer :: i, j
  do j = 1, nj-31, 32
  do i = 1, ni-31, 32
    avg = sum(q(i:i+31,j:j+31)) / 1024.0
    q(i:i+31,j:j+31) = avg
  enddo
  enddo
end

subroutine my_slope16(q, ni, nj)  ! replace values with 16x16 local average + x-y slope correction
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(INOUT) :: q
  real :: avg, xlo, xhi, ylo, yhi, dx, dy
  integer :: i, j, ii, jj
  do j = 1, nj-15, 16
  do i = 1, ni-15, 16
    avg = sum(q(i:i+15,j:j+15)) / 256.0
    xlo = sum(q(i   :i+7 ,j   :j+15  )/128.0)
    xhi = sum(q(i+8 :i+15,j   :j+15  )/128.0)
    dx = (xhi-xlo) / 8.0
    ylo = sum(q(i   :i+15  ,j   :j+7 )/128.0)
    yhi = sum(q(i   :i+15  ,j+8 :j+15)/128.0)
    dy = (yhi-ylo) / 16.0
    do jj = 0, 15
    do ii = 0, 15
      q(i+ii,j+jj) = avg + dx * (ii-7.5) + dy * (jj-7.5)
    enddo
    enddo
!     q(i:i+15,j:j+15) = avg
  enddo
  enddo
end

subroutine my_filter16(q, ni, nj)  ! replace values with 16x16 local average
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(INOUT) :: q
  real :: avg
  integer :: i, j
  do j = 1, nj-15, 16
  do i = 1, ni-15, 16
    avg = sum(q(i:i+15,j:j+15)) / 256.0
    q(i:i+15,j:j+15) = avg
  enddo
  enddo
end

subroutine my_filter8(q, ni, nj)  ! replace values with 8x8 local average
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(INOUT) :: q
  real :: avg
  integer :: i, j
  do j = 1, nj-7, 8
  do i = 1, ni-7, 8
    avg = sum(q(i:i+7,j:j+7)) / 64.0
    q(i:i+7,j:j+7) = avg
  enddo
  enddo
end

program test_compress
  use ISO_C_BINDING
  use globalstats
  implicit none
#include <misc_pack.hf>
  integer, external :: fnom, fstouv, fstnbr, fstinf, fstsui
  integer :: iun, status, nrec, key, ni, nj, nk, irec, ilev, ilen
  integer :: date,deet,npas,nbits,datyp,ip1,ip2,ip3,ig1,ig2,ig3,ig4
  integer :: swa,lng,dltf,ubc,extra1,extra2,extra3
  character(len=2) :: typvar
  character(len=4) :: nomvar, oldnam
  character(len=12) :: etiket
  character(len=1) :: grtyp
  real, dimension(:,:), pointer :: z=>NULL()
  real, dimension(:,:), pointer :: y=>NULL()
  real, dimension(:), pointer :: p=>NULL()
  real, dimension(:), pointer :: q=>NULL()
  integer :: sizep
  integer, dimension(:,:), pointer :: iz=>NULL()
  character (len=128) :: filename, varname, str_ip
  character (len=128) :: outfile
  integer :: iunout
  integer :: i, j, the_kind
  real, dimension(NSCORES) :: s
  real :: p1
  real :: maxvalue, minvalue, minabsvalue, quantum
  integer :: nmiss
  type(PackHeader) :: phead
  interface
    subroutine ieee_clip(f, n, nbits) bind(C,name='ieee_clip')
      import :: C_FLOAT, C_INT
      implicit none
      real(C_FLOAT), dimension(*), intent(INOUT) :: f
      integer(C_INT), intent(IN), value :: n, nbits
    end subroutine ieee_clip
  end interface

  write(0,*)'======= compression algorithm test ======='
  iun=0
  iunout = 0
  call get_command_argument(1,filename,ilen,status)
  if(status .ne. 0) stop
  call get_command_argument(2,varname,ilen,status)
  if(status .ne. 0) stop
  call get_command_argument(3,outfile,ilen,status)
  if(status == 0) then
    print *,"OUT = '",trim(outfile),"'"
    status = fnom(iunout,trim(outfile),'RND+STD',0)
    status = fstouv(iunout,'RND')
    print *,'iunout =',iunout,' status =',status
!     goto 888
  endif
  status = fnom(iun,trim(filename),'RND+STD+R/O+OLD',0)
  call fstopi("MSGLVL",8,0)
  if(status < 0) goto 999
  status = fstouv(iun,'RND')
  if(status < 0) goto 999
  nrec = fstnbr(iun)
  irec = 0
  ilev = 0
  oldnam='    '
  write(0,*)nrec,' records found, unit=',iun
  write(6,*)"name      min            max            avg           std"
  biases = 0
  nbiases = 0
  sizep = 0;

  key = fstinf(iun,ni,nj,nk,-1,'            ',-1,-1,-1,'  ','    ')
  do while(key >= 0)
    if(ni>10 .and. nj>10) then
      if(ni*nj*nk > sizep) then
        if(associated(p)) deallocate(p)
        if(associated(q)) deallocate(q)
        sizep = ni*nj*nk
        allocate(p(sizep))
        allocate(q(sizep))
      endif
      call fstprm(key,date,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3,  &
                  typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                  swa,lng,dltf,ubc,extra1,extra2,extra3)
!       varname(1:2) = nomvar(1:2)
      if(nomvar(1:2) == varname(1:2)) then
        irec = irec + 1
        call fstluk(p,key,ni,nj,nk)
        q(1:ni*nj) = p(1:ni*nj)  !    aucun filtre
        nmiss = float_info(q, ni, ni, nj, maxvalue, minvalue, minabsvalue)
        write(6,*) 'maxvalue, minvalue, minabsvalue :',maxvalue, minvalue, minabsvalue
        quantum = 0.0
        if(nomvar(1:2) == 'TT') quantum = 0.02  ! .02 degree C
        if(nomvar(1:2) == 'TD') quantum = 0.02
        if(nomvar(1:2) == 'ES') quantum = 0.02
        if(nomvar(1:2) == 'GZ') quantum = 0.1   ! .1 dam
        if(nomvar(1:2) == 'ZZ') quantum = 0.1   ! .1 m
        if(nomvar(1:2) == 'UU') quantum = 0.1
        if(nomvar(1:2) == 'VV') quantum = 0.1
!         quantum = 0.0
        call float_quantize_prep(12, phead, maxvalue, minvalue, quantum)
        print *,'NBITS from header =',phead%nbits,' quantum =',phead%quantum
        q(1:ni*nj) = p(1:ni*nj)  !    aucun filtre
!         call ieee_clip(q, ni*nj, 8)
!         print *,'writing ',typvar//nomvar//etiket//grtyp,date,deet,npas,ni,nj,nk
!         print *,'writing ',ig1,ig2,ig3,ig4
!         call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
!                     typvar,nomvar,'NEWFIELD',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        q(1:ni*nj) = p(1:ni*nj)  !    aucun filtre
        call my_filter8(q, ni, nj)
        call fstecr(p-q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'DIFF08',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        q(1:ni*nj) = p(1:ni*nj)  !    aucun filtre
        call my_filter16(q, ni, nj)
        call fstecr(p-q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'DIFF16',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        q(1:ni*nj) = p(1:ni*nj)  !    aucun filtre
        call my_slope16(q, ni, nj)
        call fstecr(p-q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'SLOPE16',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        q(1:ni*nj) = p(1:ni*nj)  !    aucun filtre
        call my_filter32(q, ni, nj)
        call fstecr(p-q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'DIFF32',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        q(1:ni*nj) = p(1:ni*nj)  !    aucun filtre
        call my_slope32(q, ni, nj)
        call fstecr(p-q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'SLOPE32',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call lorenzo(p, q, ni, nj)
        q(1) = 0.0
        call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'LORENZO',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call bilorentz(p, q, ni, nj)
        q(1) = 0.0
        call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'BILORENTZ',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        z(1:ni,1:nj) => p(1:ni*nj)
        y(1:ni,1:nj) => q(1:ni*nj)
!         y(1:ni,1:nj) => q(1:ni*nj)
!         write(6,*) nomvar, minval(z), maxval(z)
        call scores(s,z,y,ni,nj,1.00000,.true.,'diag:',1)
!         call scores(s,z,y,ni,nj,1.00000,.false.,'diag:',1)
        call CONVIP_plus( ip1, p1, the_kind, -1, str_ip, .false. )
        write(6,*) nomvar, s(12), s(13), s(16), s(8), p1
        if(irec == 25) exit
      endif
    endif
    key = fstsui(iun,ni,nj,nk)
  enddo
  write(6,*)'number of records processed:',irec,' out of',nrec
  write(6,*)'biases=',biases
  write(6,*)'nbiases=',nbiases
888 continue
  if(iun .ne. 0) call fstfrm(iun)
  if(iunout .ne. 0) call fstfrm(iunout)
  stop
999 continue
  write(0,*)'=== ERROR opening files ==='
  stop
end program test_compress

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

  real *8 :: suma, sumb, sumab, suma2, sumb2, bias, rms, abserr, errrel, one
  real :: errmax
  real *8 :: errrelmax, errrelavg, avga, avgb, stda, stdb, corr, slope, intercept
  integer :: i, j, n, isd
  real :: nonzero, fmina, fmaxa, span, sdl, sd, fminb, fmaxb
  integer, dimension(6) :: dig

  s = 0.0
  one = 1.0
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
  fmina = fa(1,1)
  fminb = fb(1,1)
  fmaxa = fmina
  fmaxb = fminb
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
      errrel = max(errrel, .000005_8)
      sdl = log10(max(errrel, .0001_8))
      isd = max(1, min(6,nint(-sdl)))
      dig(isd) = dig(isd)+1
      sd = sd - sdl
      errrelavg = errrelavg + errrel
      errrelmax = max( errrelmax , errrel )
    endif
    fmina = min(fmina,fa(i,j))
    fmaxa = max(fmaxa,fa(i,j))
    fminb = min(fminb,fb(i,j))
    fmaxb = max(fmaxb,fb(i,j))
  enddo
  enddo
  avga = suma / (ni*nj)
  stda = sqrt( max( 0.0_8 , (suma2/n) - ((suma/n)*(suma/n)) ) )
  avgb = sumb / (ni*nj)
  stdb = sqrt( max( 0.0_8 , (sumb2/n) - ((sumb/n)*(sumb/n)) ) )
  rms = sqrt( rms / (ni*nj) )
  if(rms <= 1.0E-20) rms = 1.0E-20
  bias = bias / (ni*nj)
  if(abs(bias) <= 1.0E-20) bias = 1.0E-20
  rms = max(abs(bias),rms)
  abserr = abserr / (ni*nj)
  if(nonzero > 0) then
    dig = nint(dig * 100.0 / nonzero)
    sd = sd / nonzero
    errrelavg = errrelavg / nonzero
  endif
  span = min( abs(fmaxa-avga) , abs(avga-fmina) )
  span = max(span, 1.0E-10)
  span = (fmaxa - fmina) / span
  if(stda*stdb .ne. 0) corr = (sumab - n*avga*avgb) / (n*stda*stdb)
  if( ( suma2/n - avga*avga) .ne. 0) slope = (sumab/n - avga*avgb) / ( suma2/n - avga*avga)
  intercept = avgb - slope * avga
  call s1scor(s(1),fa,fb,ni,nj,1,1,ni,nj,1)
!   s( 2) = rms           ! average squared error
  s( 2) = rms
  s( 2) = max(s(2), 1.0E-10)
  s( 3) = bias          ! bias
  s( 4) = abserr        ! average absolute error
  s( 5) = errmax        ! max absolute error
  s( 7) = errrelmax     ! max realtive error
  s( 6) = errrelavg     ! average relative error
  s( 8) = stda
  s( 9) = stdb
  s(10) = sd
  s(11) = corr
  s(12) = fmina
  s(13) = fmaxa
  s(14) = fminb
  s(15) = fmaxb
  s(16) = avga
  if(diag) then
    print 100,msg,' S1=',s(1),' Bias=',s(3),' RMS=',s(2),' B/R=',nint(100*abs(s(3)/s(2))), &
             ' Ae=',s(5),s(4),' Re=',s(7),s(6),' NZ=',nint(nonzero/(ni*nj)*100.0), &
             ' dSTD=',stda-stdb,' SD=',sd,' dC1=',1.0-corr,' Sp=',nint(span*.5), &
             ' A,B=',1.0-slope,intercept,dig
100 format(A,A,E8.2,A,E9.2,A,E8.2,A,I2,1H%,A,E8.2,E9.2,A,E8.2,E9.2,A,I3,1H%,A,E9.2,A,F3.1,A,E8.2,A,I5,A,E9.2,E10.2,6I4)
  endif

  nbiases(mode) = nbiases(mode) + 1
  if(bias > 0) then
    biases(mode) = biases(mode) + 1
  else
    biases(mode) = biases(mode) - 1
  endif
  return
end subroutine scores
! subroutine extrema
!   implicit none
!   integer :: i, j, maxexp1, maxexp2
!   zmin = z(1,1)
!   zmin0 = 1.0
!   zmax = z(1,1)
!   meanval = 0.0
!   do j=1,nj
!   do i=1,ni
!     zmin = min(z(i,j),zmin)
!     zmax = max(z(i,j),zmax)
!     if(z(i,j) .ne. 0.0) then
!       nonzero = nonzero + 1
!       zmin0 = min(z(i,j),zmin0)
!     endif
!     meanval = meanval + z(i,j)
!   enddo
!   enddo
!   rrange = zmax - zmin
!   if(rrange == 0.0) return
!   meanval = meanval / (ni*nj)  ! moyenne
!   span = min(zmax-meanval,meanval-zmin)
!   if(span == 0.0) span = rrange * .0001
!   span = 0.5 * (rrange / span)
!   maxexp1 = transfer(zmin,maxexp1)
!   maxexp1 = iand(ishft(maxexp1,-23),Z'000000FF')
!   maxexp2 = transfer(zmax,maxexp2)
!   maxexp2 = iand(ishft(maxexp2,-23),Z'000000FF')
!   maxexp = max(maxexp1,maxexp2)
! end subroutine extrema
