#define IN_FORTRAN_CODE__

!=========================== UNUSED CODE ===========================
#if 0
module analyze_data_mod
  use ISO_C_BINDING
  implicit none
#include <misc_pack.hf>

contains

  function array_stats_1(zi, ni, lni, nj, quantum) result(bits)
    implicit none
    integer, intent(IN), value :: ni, lni, nj
    real, dimension(*), intent(IN), target :: zi
    integer, dimension(:,:), pointer :: bits

    real, dimension(lni,64) :: z
    pointer(pz, z)
    integer, dimension(:), allocatable :: boundi, boundj
    integer, dimension(64,64) :: q
    real, dimension(64,64) :: zr
    integer :: i0, j0, i, j, ni0, nj0, nb0
    integer, dimension(2) :: t
    real, intent(IN) :: quantum

    ni0 = (ni+63)/64
    nj0 = (nj+63)/64
    allocate( bits(ni0, nj0), boundi(ni0+1), boundj(nj0+1) )
    bits = 0
    zr = 5.0
    zr(64,:) = 55.0

    boundi = [ ( (i*64)-63 , i = 1, size(boundi) ) ]
    boundi(ni0) = ni + 1 - mod(ni,64)
    boundi(ni0+1) = ni + 1
    print 1,boundi

    boundj = [ ( ((j*64)-63) , j = 1, size(boundj) ) ]
    boundj(nj0) = nj + 1 - mod(nj,64)
    boundj(nj0+1) = nj + 1
    print 1,boundj

    bits = 0
!     quantum = 1.0
    do j0 = 1, nj0
    do i0 = 1, ni0
      pz = loc( zi(boundi(i0)+(boundj(j0)-1)*lni) )
      bits(i0,j0) = float_quantize_simple(z, q, 64, 64, 64, 64, quantum, t)
    enddo
    enddo
    do j = nj0, 1, -1
      print 2, bits(:,j)
    enddo
1 format(50I5)
2 format(50I3)
  end function
end module
#endif


module globalstats
  use ISO_C_BINDING
  use analyze_data_mod
  implicit none
  integer, dimension(-1:4), save :: biases, nbiases
end module globalstats

#define NSCORES 16

subroutine smooth124(s, d, ni, nj)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: s
  real, dimension(ni,nj), intent(OUT)  :: d
#define IN_FORTRAN_CODE
#include <smooth124.h>
  call Fsmooth_124_2D(s, d, ni, ni, ni, nj)
  d = s - d
end subroutine smooth124

subroutine avg2f(p, q, r, r2, ni, nj)
  use ISO_C_BINDING
  implicit none
#define IN_FORTRAN_CODE
#include <average_2x2.h>
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT)  :: r, r2
  real, dimension((ni+1)/2,(nj+1)/2), intent(OUT)  :: q
  call average_2x2_2D_F32(p, q, ni, ni, nj)
  call expand_2x2_2D_F32(r, q, ni, ni, nj)
  r2 = r - p
end subroutine avg2f

subroutine avg2(p, q, r, r2, ni, nj)
  use ISO_C_BINDING
  implicit none
#define IN_FORTRAN_CODE
#include <average_2x2.h>
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT)  :: r, r2
  real, dimension((ni+1)/2,(nj+1)/2), intent(OUT)  :: q
  integer(C_INT32_t), dimension(ni,nj) :: ip
  integer(C_INT32_t), dimension(ni,nj) :: ir, ir2
  integer(C_INT32_t), dimension((ni+1)/2,(nj+1)/2) :: iq
!   call quant16i(p, ip, ni, nj)
  call quant12i(p, ip, ni, nj)
  call avgres_2x2_2D_I32(ip, iq, ir2, ni, ni, nj)
  call average_2x2_2D_I32(ip, iq, ni, ni, nj)
  call expand_2x2_2D_I32(ir, iq, ni, ni, nj)
  ir = ip - ir  ! residual after restore
  call population(iq, ((ni+1)/2)*((nj+1)/2), 'avg   2x2 apres')
  call population(ir, ni*nj, 'res   2x2 apres')
  call population(ir2, ni*nj, 'avres 2x2 apres')
  q = iq
  r = ir
  r2 = ir2
end subroutine avg2

subroutine quant12(p, q, ni, nj)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT) :: q
  integer(C_INT32_t), dimension(ni,nj) :: iq
  real :: pmax, pmin, fac
  pmax = maxval(p)
  pmin = minval(p)
  fac = 1.0 / (pmax-pmin)
  iq = (p - pmin) * fac * 4095 + .5
  q = iq
print *,'QUANT12, quantum =',(pmax-pmin)/4095.0
end subroutine quant12

subroutine quant16(p, q, ni, nj)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT) :: q
  integer(C_INT32_t), dimension(ni,nj) :: iq
  real :: pmax, pmin, fac
  pmax = maxval(p)
  pmin = minval(p)
  fac = 1.0 / (pmax-pmin)
  iq = (p - pmin) * fac * 65535 + .5
  q = iq
print *,'QUANT16, quantum =',(pmax-pmin)/65535.0
end subroutine quant16

subroutine quant12i(p, q, ni, nj)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  integer, dimension(ni,nj), intent(OUT) :: q
  real :: pmax, pmin, fac
  pmax = maxval(p)
  pmin = minval(p)
  fac = 1.0 / (pmax-pmin)
  q = (p - pmin) * fac * 4095 + .5
end subroutine quant12i

subroutine quant16i(p, q, ni, nj)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  integer, dimension(ni,nj), intent(OUT) :: q
  real :: pmax, pmin, fac
  pmax = maxval(p)
  pmin = minval(p)
  fac = 1.0 / (pmax-pmin)
  q = (p - pmin) * fac * 65535 + .5
end subroutine quant16i

subroutine tile(q, ni, nj, step)
  use ISO_C_BINDING
  implicit none
  integer, dimension(ni,nj), intent(INOUT) :: q
  integer(C_INT32_t), intent(IN) :: ni, nj, step
  integer :: i0, j0, mini, maxi
  do j0 = 1, nj-step, step
  do i0 = 1, ni-step, step
    mini = minval(q(i0:i0+step-1,j0:j0+step-1))
    maxi = maxval(q(i0:i0+step-1,j0:j0+step-1))
    q(i0:i0+step-1,j0:j0+step-1) = maxi - mini
!     if(mini >= 0) q(i0:i0+step-1,j0:j0+step-1) = maxi
!     if(maxi <  0) q(i0:i0+step-1,j0:j0+step-1) = mini
  enddo
  enddo
end subroutine tile

subroutine population(p, n, msg)
  use ISO_C_BINDING
  implicit none
#include <rmn/misc_operators.h>
  integer(C_INT32_t), intent(IN), dimension(n) :: p
  integer(C_INT32_t), intent(IN) :: n
  character(len=*) :: msg
  integer(C_INT32_t), dimension(34) :: pop
  integer :: i
  pop = 0
  call BitPop(p, pop, n)
  print 1, trim(msg), pop(1:21),BitEntropy(p,n,18,0)
  do i=2,34
    pop(i) = pop(i) + pop(i-1)
  enddo
  print 1, 'cumul',pop(1:20),pop(34)
1 format(A15,21I8,2X,F6.2)
end subroutine population

subroutine wavelet(p, q, ni, nj)
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT) :: q
#include <rmn/misc_operators.h>
  interface
    subroutine FDWT53i_2D_split_inplace_n(x, ni, lni, nj, levels) bind(C, name='FDWT53i_2D_split_inplace_n')
      import :: C_INT32_T
      implicit none
      integer(C_INT32_t), intent(IN), value :: ni, lni, nj, levels
      integer(C_INT32_t), dimension(lni,nj), intent(INOUT) :: x
    end subroutine FDWT53i_2D_split_inplace_n
  end interface
  integer(C_INT32_t), dimension(ni,nj) :: iq
  integer(C_INT32_t), dimension(34) :: pop
  real :: pmax, pmin, fac
  integer :: nz, i, j
  call quant16i(p, iq, ni, nj)
  pmax = maxval(p)
  pmin = minval(p)
  call population(iq, ni*nj, 'wavelet avant')
!   print *,'entropy =',BitEntropy(iq, ni*nj, 16, 0)
  call FDWT53i_2D_split_inplace_n(iq, ni, ni, nj, 3)
  q = iq
  call population(iq, ni*nj, 'wavelet apres')
  call population(iq(1     :ni/2 ,1     :nj/2), ni*nj/4, 'wavelet corner')
  call population(iq(ni/2+1:ni   ,1     :nj/2), ni*nj/4, 'wavelet corner')
  call population(iq(1     :ni/2 ,nj/2+1:nj  ), ni*nj/4, 'wavelet corner')
  call population(iq(ni/2+1:ni   ,nj/2+1:nj  ), ni*nj/4, 'wavelet corner')
  call tile(iq, ni, nj, 4)
  call population(iq, ni*nj, 'wavelet tiled')
!   iq(1:(ni+7)/8 , 1:(nj+7)/8) = 0
!   call population(iq, ni*nj, 'wavelet q0')
!   q = iq
end subroutine wavelet

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

subroutine lorenzo12(p, q, ni, nj) ! Lorenzo predictor
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT) :: q
  integer :: i, j
  integer(C_INT32_t), dimension(34) :: pop
  integer(C_INT32_t), dimension(ni,nj) :: iq
#include <rmn/misc_operators.h>
  call quant12i(p, iq, ni, nj)
  call population(iq, ni*nj, 'lorenzo12 avant')
!   print *,'entropy =',BitEntropy(iq, ni*nj, 16, 0)
  do j = nj, 2, -1
    do i = ni, 2, -1
      q(i,j) = iq(i,j) - (iq(i-1,j) + iq(i,j-1) - iq(i-1,j-1))
    enddo
    q(1,j) = iq(1,j) - iq(1,j-1)
  enddo
  do i = ni, 2, -1                   ! row 1
    q(i,1) = iq(i,1) - iq(i-1,1)
  enddo
  q(1,1) = iq(1,1)
  iq = q
  call population(iq, ni*nj, 'lorenzo12 apres')
!   call population(iq(1     :ni/2 ,1     :nj/2), ni*nj/4, 'lorenzo corner')
!   call population(iq(ni/2+1:ni   ,1     :nj/2), ni*nj/4, 'lorenzo corner')
!   call population(iq(1     :ni/2 ,nj/2+1:nj  ), ni*nj/4, 'lorenzo corner')
!   call population(iq(ni/2+1:ni   ,nj/2+1:nj  ), ni*nj/4, 'lorenzo corner')
!   print *,'entropy =',BitEntropy(iq, ni*nj, 16, 0)
  call tile(iq, ni, nj, 4)
  call population(iq, ni*nj, 'lorenzo12 tiled')
end

subroutine lorenzo(p, q, ni, nj) ! Lorenzo predictor
  use ISO_C_BINDING
  implicit none
  integer, intent(IN) :: ni, nj
  real, dimension(ni,nj), intent(IN)  :: p
  real, dimension(ni,nj), intent(OUT) :: q
  integer :: i, j
  integer(C_INT32_t), dimension(34) :: pop
  integer(C_INT32_t), dimension(ni,nj) :: iq
#include <rmn/misc_operators.h>
  call quant16i(p, iq, ni, nj)
  call population(iq, ni*nj, 'lorenzo avant')
!   print *,'entropy =',BitEntropy(iq, ni*nj, 16, 0)
  do j = nj, 2, -1
    do i = ni, 2, -1
      q(i,j) = iq(i,j) - (iq(i-1,j) + iq(i,j-1) - iq(i-1,j-1))
    enddo
    q(1,j) = iq(1,j) - iq(1,j-1)
  enddo
  do i = ni, 2, -1                   ! row 1
    q(i,1) = iq(i,1) - iq(i-1,1)
  enddo
  q(1,1) = iq(1,1)
  iq = q
  call population(iq, ni*nj, 'lorenzo apres')
!   call population(iq(1     :ni/2 ,1     :nj/2), ni*nj/4, 'lorenzo corner')
!   call population(iq(ni/2+1:ni   ,1     :nj/2), ni*nj/4, 'lorenzo corner')
!   call population(iq(1     :ni/2 ,nj/2+1:nj  ), ni*nj/4, 'lorenzo corner')
!   call population(iq(ni/2+1:ni   ,nj/2+1:nj  ), ni*nj/4, 'lorenzo corner')
!   print *,'entropy =',BitEntropy(iq, ni*nj, 16, 0)
  call tile(iq, ni, nj, 4)
  call population(iq, ni*nj, 'lorenzo tiled')
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
! #include <misc_pack.hf>
  integer, external :: fnom, fstouv, fstnbr, fstinf, fstsui
  integer :: iun, status, nrec, key, nk, irec, ilev, ilen
  integer, target :: ni, nj, ninj
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
  real, dimension(:), pointer :: r=>NULL()
  real, dimension(:), pointer :: r2=>NULL()
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
  interface
    function c_creat(name, mode) result(fd) bind(C,name='creat')
      import :: C_CHAR, C_INT
      implicit none
      character(C_CHAR), dimension(*), intent(IN) :: name
      integer(C_INT), intent(IN), value :: mode
      integer(C_INT) :: fd
    end function
    function c_write(fd, buf, cnt) result(nc) bind(C, name='write')
      import :: C_SIZE_T, C_PTR, C_INT
      implicit none
      integer(C_INT), intent(IN), value :: fd
      type(C_PTR), intent(IN), value :: buf
      integer(C_SIZE_T), intent(IN), value :: cnt
      integer(C_SIZE_T) :: nc
    end function
    function c_read(fd, buf, cnt) result(nc) bind(C, name='read')
      import :: C_SIZE_T, C_PTR, C_INT
      implicit none
      integer(C_INT), intent(IN), value :: fd
      type(C_PTR), intent(IN), value :: buf
      integer(C_SIZE_T), intent(IN), value :: cnt
      integer(C_SIZE_T) :: nc
    end function
    function c_close(fd) result(status) bind(C,name='close')
      import :: C_INT
      implicit none
      integer(C_INT), intent(IN), value :: fd
      integer(C_INT) :: status
    end function
  end interface
  integer, dimension(:,:), pointer :: bits0
  integer :: fd, fdstatus, fdmode, ipkind
  real :: ipvalue
  integer(C_SIZE_T) :: nc
  character(len=128) :: c_fname, ipstring

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
        if(associated(r)) deallocate(r)
        if(associated(r2)) deallocate(r2)
        sizep = ni*nj*nk
        allocate(p(sizep))
        allocate(q(sizep))
        allocate(r(sizep))
        allocate(r2(sizep))
      endif
      call fstprm(key,date,deet,npas,ni,nj,nk,nbits,datyp,ip1,ip2,ip3,  &
                  typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
                  swa,lng,dltf,ubc,extra1,extra2,extra3)
!       varname(1:2) = nomvar(1:2)
      if(nomvar(1:2) == varname(1:2)) then
        irec = irec + 1
        call fstluk(p,key,ni,nj,nk)
! ==========================================================================
        call CONVIP_plus( ip1, ipvalue, ipkind, -1, ipstring, .true. )  ! convert ip1
        do i = 1, 20
          if(ipstring(i:i) == ' ') then
            ipstring(i:i) = '_'
            exit
          endif
        enddo
        write(c_fname,1111)'RAW/',trim(nomvar),'_'//trim(ipstring)   ! create file name
1111 format(A,A,A,F10.10)
        fdmode = INT(O'777')
        fd = c_creat(trim(c_fname)//achar(0), fdmode) ! create raw file
        if(fd > 0) then
          ninj = 2         ! 2D data
          nc = 4
          nc = c_write(fd, C_LOC(ninj), nc)
          nc = 4
          nc = c_write(fd, C_LOC(ni), nc)
          nc = 4
          nc = c_write(fd, C_LOC(nj), nc)
          nc = 4 * ni * nj
          nc = c_write(fd, C_LOC(p(1)), nc)
          print *,'INFO, wrote',nc+16,' bytes into ',trim(c_fname)
          ninj = ni * nj ;
          nc = 4
          nc = c_write(fd, C_LOC(ninj), nc)
          fdstatus = c_close(fd)
        else
          print *,'ERROR creating '//trim(c_fname)
        endif
! ==========================================================================
        q(1:ni*nj) = p(1:ni*nj)  !    aucun filtre
        nmiss = float_info(q, ni, ni, nj, maxvalue, minvalue, minabsvalue)
        write(6,*) 'maxvalue, minvalue, minabsvalue :',maxvalue, minvalue, minabsvalue
        quantum = 0.0
        if(nomvar(1:2) == 'TT') quantum = 0.01  ! .01 degree C
        if(nomvar(1:2) == 'TD') quantum = 0.02
        if(nomvar(1:2) == 'ES') quantum = 0.02
        if(nomvar(1:2) == 'GZ') quantum = 0.1   ! .1 dam
        if(nomvar(1:2) == 'ZZ') quantum = 0.1   ! .1 m
        if(nomvar(1:2) == 'UU') quantum = 0.1
        if(nomvar(1:2) == 'VV') quantum = 0.1
        if(nomvar(1:2) == 'WW') quantum = 0.01
!         quantum = 0.0

bits0 => array_stats_1(p, ni, ni, nj, quantum)

#if 0
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
        call lorenzo12(p, q, ni, nj)
        q(1) = 0.0
        call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'LORENZO12',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call lorenzo(p, q, ni, nj)
        q(1) = 0.0
        call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'LORENZO',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call smooth124(p, r2, ni, nj)
        call fstecr(r2, r2, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'SMOOTH124',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call bilorentz(p, q, ni, nj)
        q(1) = 0.0
        call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'BILORENTZ',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call wavelet(p, q, ni, nj)
        call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'WAVELET',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call quant12(p, q, ni, nj)
        call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'QUANT12',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call quant16(p, q, ni, nj)
        call fstecr(q, q, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'QUANT16',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call avg2(p, q, r, r2, ni, nj)
        call fstecr(q, q, -32, iunout, date,deet,npas,ni/2,nj/2,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'AVG2x2','X',0,0,0,0, 5, .false. )
        call fstecr(r, r, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'RES2X2',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call fstecr(r2, r2, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'AVGRES',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        call avg2f(p, q, r, r2, ni, nj)
        call fstecr(q, q, -32, iunout, date,deet,npas,ni/2,nj/2,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'AVG2x2F','X',0,0,0,0, 5, .false. )
        call fstecr(r2, r2, -32, iunout, date,deet,npas,ni,nj,nk,ip1,ip2,ip3,  &
                    typvar,nomvar,'AVGRESF',grtyp,ig1,ig2,ig3,ig4, 5, .false. )
        z(1:ni,1:nj) => p(1:ni*nj)
        y(1:ni,1:nj) => q(1:ni*nj)
!         y(1:ni,1:nj) => q(1:ni*nj)
!         write(6,*) nomvar, minval(z), maxval(z)
        call scores(s,z,y,ni,nj,1.00000,.true.,'diag:',1)
!         call scores(s,z,y,ni,nj,1.00000,.false.,'diag:',1)
        call CONVIP_plus( ip1, p1, the_kind, -1, str_ip, .false. )
        write(6,*) nomvar, s(12), s(13), s(16), s(8), p1
#endif
        if(irec == 10) exit
!         if(irec == 25) exit
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
