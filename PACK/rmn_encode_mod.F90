! Copyright (C) 2023  Environnement et Changement climatique Canada
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
!     M. Valin,   Environnement et Changement climatique Canada, 2023
!
! interface to the C bit stream encoder/decoder
!
module rmn_encode
  use ISO_C_BINDING
  use rmn_stream_pack
  implicit none

  interface
!   uint32_t stream_get_block_8x8(uint32_t * restrict src, int lni, uint32_t * restrict block);
    function stream_get_block_8x8(src, lni, blk) result(nbits) bind(C, name='')
      import :: C_INT32_T
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN) :: src
      integer(C_INT32_T), intent(IN), value :: lni
      integer(C_INT32_T), dimension(*), intent(OUT) :: blk
      integer(C_INT32_T) :: nbits
    end function
!   void stream_encode_init(bitstream *bstream, void *buffer, size_t bufsize);
    subroutine stream_encode_init(stream, buffer, bufsize) bind(C, name='stream_encode_init')
      import :: C_PTR, C_SIZE_T
      implicit none
      type(C_PTR), intent(IN), value :: stream, buffer
      integer(C_SIZE_T), intent(IN), value :: bufsize
    end subroutine
!   uint32_t stream_encode_ublock(uint32_t *src, int nx, int ny, int nbits, bitstream *bstream);
    function stream_encode_ublock(src, nx, ny, nbits, stream) result(total_size) bind(C, name='stream_encode_ublock')
      import :: C_INT32_T, C_PTR
      implicit none
      integer(C_INT32_T), dimension(*), intent(IN) :: src
      integer(C_INT32_T), intent(IN), value :: nx, ny, nbits
      type(C_PTR), intent(IN), value :: stream
      integer(C_INT32_T) :: total_size
    end function
!   uint32_t stream_decode_ublock(uint32_t *dst, bitstream *bstream);
    function stream_decode_ublock(dst, stream) result(dims) bind(C, name='stream_decode_ublock')
      import :: C_INT32_T, C_PTR
      implicit none
      integer(C_INT32_T), dimension(*), intent(OUT) :: dst
      type(C_PTR), intent(IN), value :: stream
      integer(C_INT32_T) :: dims
    end function
  end interface
end module
