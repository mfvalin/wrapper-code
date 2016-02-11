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
  status = fnom(iun,'test_file.fst','RND+STD+R/O+OLD',0)
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
      call fstluk(z,key,ni,nj,nk)
      ilev = ilev + 1
      call fstprm(key,date,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3,  &
                  typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                  swa,lng,dltf,ubc,extra1,extra2,extra3)
      if(nomvar .ne. oldnam) then
        if(oldnam .ne. '    ') write(0,*)'NAME=',oldnam,' levels=',ilev
        oldnam = nomvar
        ilev = 0
      endif
      call test_compression(z,ni,nj,nomvar)
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
subroutine test_compression(z,ni,nj,nomvar)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni,nj
  real, dimension(ni,nj), intent(IN) :: z
  character(len=4), intent(IN) :: nomvar
  interface
!   int TO_jpeg2000(unsigned char *cin,int width,int height,int nbits,
!                  int ltype, int ratio, int retry, char *outjpc, 
!                  int jpclen)
    function enc_jpeg(cin,width,height,nbits,ltype,ratio,retry,cout,jpclen) result(status) bind(C,name='TO_jpeg2000') 
      import
      implicit none
      integer(C_INT), intent(IN), value :: width,height,nbits,ltype,retry,jpclen
      real(C_FLOAT), intent(IN), value :: ratio
      type(C_PTR), intent(IN), value :: cin, cout
      integer(C_INT) :: status
    end function enc_jpeg
!   int FROM_jpeg2000(char *injpc,int bufsize,int *outfld)
    function dec_jpeg(cin,bufsize,cout) result(status) bind(C,name='FROM_jpeg2000') 
      import
      implicit none
      integer(C_INT), intent(IN), value :: bufsize
      type(C_PTR), intent(IN), value :: cin, cout
      integer(C_INT) :: status
    end function dec_jpeg

  end interface

  integer, dimension(ni,nj) :: iz
  real :: the_min, the_max, error
  integer :: nbits, nbts, minv, nbits2, mn1, mx1, nbytes, njpeg
  integer, external :: dwt_quantize, dwt_lorenzo, dwt_nbits, dwt_pack
  integer(C_CHAR), dimension(10000000), target :: buf
  integer(C_CHAR), dimension(10000000), target :: jbuf

  error = 0.0
  nbts = 16
  nbits = dwt_quantize(z,iz,ni*nj,nbts,error,the_min,the_max)
  mn1 = minval(iz)
  mx1 = maxval(iz)
  minv = dwt_lorenzo(iz,ni,nj,.true.)
  iz = iz - minv
  nbits2=dwt_nbits(maxval(iz))
  nbytes = dwt_pack(iz,ni*nj,buf,10000000,nbits)
  njpeg = enc_jpeg(c_loc(buf),ni,nj,nbits2,0,1.0,0,c_loc(jbuf),10000000)
  
  write(0,100)nomvar,nbits,the_min,the_max,mn1,mx1,minval(iz),maxval(iz),nbits2,minv,nbytes,njpeg,ni*nj*4
100 format(1X,A4,1X,I3,2G12.5,10I8)
  return
end subroutine test_compression
