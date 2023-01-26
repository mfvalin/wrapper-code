! Copyright (C) 2022  Environnement et Changement climatique Canada
! 
! This is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
! 
! This software is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! Author:
!     M. Valin,   Environnement et Changement climatique Canada, 2022
!
module analyze_data_mod
  use ISO_C_BINDING
  use lorenzo_mod
  use rmn_misc_operators
  use rmn_zfpx
  implicit none
#include <rmn/misc_pack.hf>
#define IN_FORTRAN_CODE
! #include <misc_operators.h>
#include <average_2x2.h>
#include <misc_analyze.h>
#if 0
  interface
    ! void AnalyzeCompressionErrors(float *fa, float *fb, int np, float small, char *str)
    subroutine AnalyzeCompressionErrors(fa, fb, np, small, str) BIND(C, name='AnalyzeCompressionErrors')
      import ::C_FLOAT, C_INT32_T, C_CHAR
      implicit none
      real(C_FLOAT), dimension(*), intent(IN) :: fa, fb
      integer(C_INT32_T), intent(IN), value :: np
      real(C_FLOAT), intent(IN), value :: small
      character(C_CHAR), dimension(*), intent(IN) :: str
    end subroutine
  end interface
#endif
contains

  function array_stats_1(zi, ni, lni, nj, quantum) result(bits)
    implicit none
    integer, intent(IN), value :: ni, lni, nj
    real, dimension(*), intent(IN), target :: zi
    integer, dimension(:,:), pointer :: bits

    real, dimension(lni,64) :: z
    real, dimension(ni,nj) :: zzi
    real, dimension(:,:), pointer :: pzzi
    pointer(pz, z)
    integer, dimension(:), allocatable :: boundi, boundj
    integer, dimension(64,64) :: q, ql, tnb, tnb2, tnbo, residu, tnbr, qr
    integer, dimension(32,32) :: avg2x2
    real, dimension(64,64) :: zr
    integer :: i0, j0, i, j, ni0, nj0, nb0, nbits, i1, j1, nbits0, bsz, maxbits, totsaved, errors
    integer, dimension(2) :: t
    real, intent(IN) :: quantum
    integer :: nil, njl   ! local ni and nj (<= 64)
    integer, dimension(:,:), pointer :: bitsl
    integer, dimension(32) :: tbl, tblo
    integer, dimension(lni,nj) :: qqq
    integer, dimension(64,64) :: qzfp
    type(C_PTR) :: pzfp
    integer :: Ssize
    integer, dimension(34) :: bitpop32

    print *,'Calling ZfpCompressReal'
    i = ZfpCompress_set_debug(1)
    pzfp = ZfpCompressReal(zi, ni, nj, 1, quantum, 0, Ssize)
    print *,'Ssize =',Ssize
    i = ZfpExpandReal(zzi, ni, nj, 1, pzfp, Ssize)
    pzzi(1:ni,i:nj) => zi(1:ni*nj)
    print *,'maxerr = ',maxval(abs(pzzi-zzi))

    ni0 = (ni+63)/64
    nj0 = (nj+63)/64
    allocate( bits(ni0, nj0), boundi(ni0+1), boundj(nj0+1) )
    allocate( bitsl(ni0, nj0))   ! lorenzo
    bits = 0
    bitsl = 0
    zr = 5.0
    zr(64,:) = 55.0
    i = float_quantize_simple(zi, qqq, ni, ni, ni, nj, quantum, t)
    print *, 'min, max, quantum =', minval(zi(1:ni*nj)), maxval(zi(1:ni*nj)), quantum
    print *, 'minq, maxq, maxbits =', t(1), t(2), i

    boundi = [ ( (i*64)-63 , i = 1, size(boundi) ) ]
!     boundi(ni0) = ni + 1 - mod(ni,64)
    boundi(ni0+1) = ni + 1
!     print 1,boundi
!     print 1,(boundi(i+1)-boundi(i) , i = 1, ni0)

    boundj = [ ( ((j*64)-63) , j = 1, size(boundj) ) ]
!     boundj(nj0) = nj + 1 - mod(nj,64)
    boundj(nj0+1) = nj + 1
!     print 1,boundj
!     print 1,(boundj(j+1)-boundj(j) , j = 1, nj0)

    bits = 0
!     quantum = 1.0
    bitpop32 = 0 ! set bit population to 0
    do j0 = 1, nj0
      njl = boundj(j0+1) - boundj(j0)                ! dimension along j of local tile
      do i0 = 1, ni0
        pz = loc( zi(boundi(i0)+(boundj(j0)-1)*lni) )  ! base address of local tile
  !       function float_quantize_simple(z, q, ni, lniz, lniq, nj, quantum, t) result(nbits)
        nil = boundi(i0+1) - boundi(i0)                ! dimension along i of local tile
        bits(i0,j0) = float_quantize_simple(z, q, nil, 64, 64, njl, quantum, t)
        ql = 0
        call lorenzopredict(q, ql, nil, 64, 64, njl)
!         i = to_zigzag(ql, ql, 4096)
        do j = 1, njl
        do i = 1, nil
          tnb(i,j)  = BitsNeeded_u32(to_zigzag_32(ql(i,j)))
          bitpop32(tnb(i,j)) = bitpop32(tnb(i,j)) + 1
          tnbo(i,j) = BitsNeeded_u32(to_zigzag_32(q(i,j)))
        enddo
        enddo
        if(i0 == 9 .and. j0 == 7)then  ! tile(9,7)
          qzfp = 0
!           do j = 1, 64
!           do i = 1, 64
!             q(i,j) = i**2 - 2*i*j - j*j ! + 1 * (rand() - .5)
!           enddo
!           enddo
          call lorenzopredict(q, ql, nil, 64, 64, njl)
!           call lorenzopredict2_f(q, qzfp, 64, 64, 64, 64)
!           call zfpx_gather_64_64(q, 64, qzfp, 4)
!           do j = 57, 1, -8
!             print 32, '<',j,q(1:16,j),q(17:63:3,j)
!             print 32, '+',j,qzfp(1:16,j),qzfp(17:63:3,j)
! !             print 32, '+',j,qzfp(1:16,j),qzfp(17:32,j)
! !             print 32, '|',j,qzfp(33:48,j),qzfp(49:64,j)
!             print 32, '-',j,ql(1:16,j),ql(17:63:3,j)
!           enddo
          qr=999
          call lorenzounpredict(qr, ql, nil, 64, 64, njl)
          errors = sum(qr-q)
          print *,'lorenzo predict<->unpredict total error =', errors
!           call avgres_2x2_2D_I32(q, avg2x2, residu, 64, 64, 64)
!           do j = 1, 64
!           do i = 1, 64
!               tnbr(i,j) = BitsNeeded_u32(to_zigzag_32(residu(i,j)))
!           enddo
!           enddo
!           totsaved = 0
!           do j1 = 57, 1, -8
!             do i1 = 1, 56, 8
!               tbl = 0
!               tblo = 0
!               tbl(17) = maxval(tnb(i1:i1+7,j1:j1+7))
!               tbl(16) = maxval(tnbr(i1:i1+7,j1:j1+7))
!               tbl(15) = maxval(tnbo(i1:i1+7,j1:j1+7))
!               maxbits = ((tbl(17)+1) / 2)
!               tnb2 = 1 + (tnb -1) / maxbits
!               do j = 0, 7
!               do i = 0, 7
!                 tbl(tnb2(i1+i,j1+j)) = tbl(tnb2(i1+i,j1+j)) + 1
!                 tblo(tnbo(i1+i,j1+j)) = tblo(tnbo(i1+i,j1+j)) + 1
!                 tbl(18) = tbl(18)+1
!               enddo
!               enddo
!               print 2,tbl(1:18), maxbits, tbl(17)*64 - 64 - tbl(1)*maxbits - sum(tbl(2:4))*tbl(17)  !, tblo(1:16)
!               totsaved = totsaved + (tbl(17)*64 - 64 - tbl(1)*maxbits - sum(tbl(2:4))*tbl(17))
!             enddo
!             print *,""
!           enddo
!           print *,"totsaved = ", totsaved/64
          do j = 64, 1, -1
!             print 3,tnb(:,j)
!             print 3,tnbo(:,j)
!             print 3,tnbo(:,j) - tnb(:,j)
!             print 3,tnb(:,j)-tnbr(:,j)
!             if(mod(j,8) == 1) print *,""
          enddo
!           print *,'tnbo-tnbr, tnbo-tnb, tnbr-tnb', sum(tnbo-tnbr), sum(tnbo-tnb), sum(tnbr-tnb)
        endif
        nbits = 0
        nbits0 = 0
        bsz = 8
        do j1 = 1, njl, bsz
        do i1 = 1, nil, bsz
          maxbits = maxval(tnb(i1:min(i1+bsz-1,nil),j1:min(j1+bsz-1,njl)))
          do j = j1, min(j1+bsz-1,njl)
          do i = i1, min(i1+bsz-1,nil)
            if(tnb(i,j) > (maxbits+1)/2) then
              nbits = nbits + maxbits
            else
              nbits = nbits + (maxbits+1)/2
            endif
          enddo
          enddo
          nbits = nbits + bsz*bsz
!           nbits = nbits + maxbits * bsz*bsz
          nbits0 = nbits0 + 1
        enddo
        enddo
        bitsl(i0,j0) = nbits + nbits0 * 32
      enddo
    enddo
    print *,"=== needed bits per point ==="
    print '(34I8)', bitpop32(1:16)
    print *,"quantized", sum(bits(1:ni0,1:nj0)) *1.0 / (ni0*nj0)
    print *,"lorenzo  ",sum(bitsl(1:ni0,1:nj0)) *1.0 / (ni*nj)
!     do j = nj0, 1, -1
!       print 2, bits(:,j)
!     enddo
!     do j = nj0, 1, -1
!       print 2, bitsl(:,j)
!     enddo
1 format(64I5)
2 format(64I4)
3 format(8(8I2,1X))
4 format(8(8Z2,1X))
32 format(A,I3,3I5,13I5,3x,16I5)
64 format(A,I3,3I5,13I5,3x,16I5,/,10x,16I5,3x,16I5)
  end function
end module

subroutine analyze_self_test
  use analyze_data_mod
  implicit none
  integer, parameter :: NI = 2048
  integer, parameter :: NJ = 1024
  real, dimension(NI, NJ) :: zi, zo
  integer, dimension(NI, NJ) :: q
  integer, dimension(:,:), pointer :: bits0
  integer :: i, j
  integer, dimension(2) :: t
  real, dimension(2) :: tf
  real :: x, y, v, quantum, r
! float x = 2.0 * i / nx;
! float y = 2.0 * j / ny;
! float z = 2.0 * k / nz;
! float v = (x * x + y * y + z * z) + .002f ;

  quantum = .001
  do j = 1, NJ
  do i = 1, NI
    x = 2.0 * i / ni
    y = 2.0 * j / nj
    v = x*x + y*y + .002
    call random_number(r)
    zi(i,j) = v + (r - .5) * quantum * 5
!     zi(i,j) = sqrt((i-ni*.5)**2 + (j-nj*.5)**2) - 1.0
!     zi(i,j) = ((i-ni*.5)**2 + (j-nj*.5)**2)
  enddo
  enddo
  bits0 => array_stats_1(zi, NI, NI, NJ, 0.001)
  print *,'size of bits0 =', size(bits0, 1), 'x',  size(bits0, 2)
  i = float_quantize_simple(zi, q, ni, ni, ni, nj, quantum, t)
  call float_unquantize_simple(zo, q, ni, ni, ni, nj, quantum, tf)
  call AnalyzeCompressionErrors(zi, zo, ni*nj, 0.0, "self"//achar(0))
  print *,'================================================================'
end
