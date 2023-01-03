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
! interface to the C bit stream packer/unpacker
!
module rmn_stream_pack
  use ISO_C_BINDING
  implicit none

  interface
    ! STATIC inline void  LeStreamInit(bitstream *p, uint32_t *buffer){
    subroutine LeStreamInit(stream, buffer) bind(C, name='LeStreamInit')
      import :: C_PTR
      implicit none
      type(C_PTR), intent(IN), value :: stream, buffer
    end subroutine
    subroutine BeStreamInit(stream, buffer) bind(C, name='BeStreamInit')
      import :: C_PTR
      implicit none
      type(C_PTR), intent(IN), value :: stream, buffer
    end subroutine
    ! void  LeStreamInsert(bitstream *p, uint32_t *w32, int nbits, int nitems);
    subroutine LeStreamInsert(stream, what, nbits, nitems) bind(C, name='LeStreamInsert')
      import :: C_PTR, C_INT32_T
      implicit none
      type(C_PTR), intent(IN), value :: stream, what
      integer(C_INT32_T), intent(IN), value :: nbits, nitems
    end subroutine
    subroutine BeStreamInsert(stream, what, nbits, nitems) bind(C, name='BeStreamInsert')
      import :: C_PTR, C_INT32_T
      implicit none
      type(C_PTR), intent(IN), value :: stream, what
      integer(C_INT32_T), intent(IN), value :: nbits, nitems
    end subroutine
    ! void  LeStreamXtract(bitstream *p, uint32_t *w32, int nbits, int n);
    subroutine LeStreamXtract(stream, what, nbits, nitems) bind(C, name='LeStreamXtract')
      import :: C_PTR, C_INT32_T
      implicit none
      type(C_PTR), intent(IN), value :: stream, what
      integer(C_INT32_T), intent(IN), value :: nbits, nitems
    end subroutine
    subroutine BeStreamXtract(stream, what, nbits, nitems) bind(C, name='BeStreamXtract')
      import :: C_PTR, C_INT32_T
      implicit none
      type(C_PTR), intent(IN), value :: stream, what
      integer(C_INT32_T), intent(IN), value :: nbits, nitems
    end subroutine
    ! STATIC inline void  LeStreamFlush(bitstream *p){
    subroutine LeStreamFlush(stream) bind(C, name='LeStreamFlush')
      import :: C_PTR
      implicit none
      type(C_PTR), intent(IN), value :: stream
    end subroutine
    subroutine BeStreamFlush(stream) bind(C, name='BeStreamFlush')
      import :: C_PTR
      implicit none
      type(C_PTR), intent(IN), value :: stream
    end subroutine
    ! STATIC inline void  LeStreamRewind(bitstream *p){
    subroutine LeStreamRewind(stream) bind(C, name='LeStreamRewind')
      import :: C_PTR
      implicit none
      type(C_PTR), intent(IN), value :: stream
    end subroutine
    subroutine BeStreamRewind(stream) bind(C, name='BeStreamRewind')
      import :: C_PTR
      implicit none
      type(C_PTR), intent(IN), value :: stream
    end subroutine
  end interface

  private :: stream_create_le, stream_create_be
  interface stream_create
    module procedure stream_create_le
    module procedure stream_create_be
  end interface
  private :: stream_get_le, stream_get_be
  interface stream_get
    module procedure stream_get_le
    module procedure stream_get_be
  end interface
  private :: stream_put_le, stream_put_be
  interface stream_put
    module procedure stream_put_le
    module procedure stream_put_be
  end interface
  private :: stream_flush_le, stream_flush_be
  interface stream_flush
    module procedure stream_flush_le
    module procedure stream_flush_be
  end interface
  private :: stream_rewind_le, stream_rewind_be
  interface stream_rewind
    module procedure stream_rewind_le
    module procedure stream_rewind_be
  end interface

  private :: pack32_put_init, pack32_put_fast, pack32_put_check, pack32_put_flush

  ! the followinf type layout MUST MATCH the struct bitstream typedef from bi_endian_pack.h
  type, bind(C) :: c_bitstream    ! generic C bit stream (struct bitstream)
    integer(C_INT64_T) :: accum   ! 64 bit unsigned bit accumulator
    integer(C_INT32_T) :: insert  ! # of bits used in accumulator (0 <= insert <= 64)
    integer(C_INT32_T) :: xtract  ! # of bits extractable from accumulator (0 <= xtract <= 64)
    type(C_PTR)        :: stream  ! pointer into packed stream (both insert and extract mode)
    type(C_PTR)        :: start   ! pointer to start of stream data storage
    type(C_PTR)        :: stop    ! pointer to end of stream data storage (1 byte beyond stream buffer end)
  end type

  type :: le_stream               ! little endian type stream
    type(c_bitstream) :: s
  contains
    procedure, PASS :: create => stream_create_le
    procedure, PASS :: flush  => stream_flush_le
    procedure, PASS :: rewind => stream_rewind_le
    procedure, PASS :: put    => stream_put_le
    procedure, PASS :: get    => stream_get_le
  end type

  type :: be_stream               ! big endian type stream
    type(c_bitstream) :: s
  contains
    procedure, PASS :: create => stream_create_be
    procedure, PASS :: flush  => stream_flush_be
    procedure, PASS :: rewind => stream_rewind_be
    procedure, PASS :: put    => stream_put_be
    procedure, PASS :: get    => stream_get_be
  end type
contains

! create a little endian stream (usually for writing)
subroutine stream_create_le(self, mem)
  implicit none
  class(le_stream), intent(INOUT) :: self
  integer(C_INT32_T), dimension(*), intent(IN), target :: mem
#if ! defined(__GFORTRAN__) && ! defined(__INTEL_COMPILER)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: mem
!$PRAGMA IGNORE_TKR mem
!DIR$ IGNORE_TKR mem
!IBM* IGNORE_TKR mem
#endif
  type(C_PTR) :: t, m
  t = transfer( LOC(self%s) , t)
  m = C_LOC(mem)
  call BeStreamInit(t, m)
end subroutine

! create a big endian stream (usually for writing)
subroutine stream_create_be(self, mem)
  implicit none
  class(be_stream), intent(INOUT) :: self
  integer(C_INT32_T), dimension(*), intent(IN), target :: mem
#if ! defined(__GFORTRAN__) && ! defined(__INTEL_COMPILER)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: mem
!$PRAGMA IGNORE_TKR mem
!DIR$ IGNORE_TKR mem
!IBM* IGNORE_TKR mem
#endif
  type(C_PTR) :: t, m
  t = transfer( LOC(self%s) , t)
  m = C_LOC(mem)
  call BeStreamInit(t, m)
end subroutine

! rewind a little endian stream for reading
subroutine stream_rewind_le(self)
  implicit none
  class(le_stream), intent(INOUT) :: self
  type(C_PTR) :: t
  t = transfer( LOC(self%s) , t)
  call LeStreamRewind(t)
end subroutine

! rewind a big endian stream for reading
subroutine stream_rewind_be(self)
  implicit none
  class(be_stream), intent(INOUT) :: self
  type(C_PTR) :: t
  t = transfer( LOC(self%s) , t)
  call BeStreamRewind(t)
end subroutine

! flush write leftovers in a little endian stream
subroutine stream_flush_le(self)
  implicit none
  class(le_stream), intent(INOUT) :: self
  type(C_PTR) :: t
  t = transfer( LOC(self%s) , t)
  call LeStreamFlush(t)
end subroutine

! flush write leftovers in a big endian stream
subroutine stream_flush_be(self)
  implicit none
  class(be_stream), intent(INOUT) :: self
  type(C_PTR) :: t
  t = transfer( LOC(self%s) , t)
  call BeStreamFlush(t)
end subroutine

! extract data from a little endian stream
subroutine stream_get_le(self, what, nbits, nitems)
  implicit none
  class(le_stream), intent(INOUT) :: self
  integer(C_INT32_T), dimension(*), intent(OUT), target  :: what
#if ! defined(__GFORTRAN__) && ! defined(__INTEL_COMPILER)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: what
!$PRAGMA IGNORE_TKR what
!DIR$ IGNORE_TKR what
!IBM* IGNORE_TKR what
#endif
  integer(C_INT32_T), intent(IN), value :: nbits, nitems
  type(C_PTR) :: t
  t = transfer( LOC(self%s) , t)
  call LeStreamXtract(t, C_LOC(what), nbits, nitems)
end subroutine

! extract data from a big endian stream
subroutine stream_get_be(self, what, nbits, nitems)
  implicit none
  class(be_stream), intent(INOUT) :: self
  integer(C_INT32_T), dimension(*), intent(OUT), target  :: what
#if ! defined(__GFORTRAN__) && ! defined(__INTEL_COMPILER)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: what
!$PRAGMA IGNORE_TKR what
!DIR$ IGNORE_TKR what
!IBM* IGNORE_TKR what
#endif
  integer(C_INT32_T), intent(IN), value :: nbits, nitems
  type(C_PTR) :: t
  t = transfer( LOC(self%s) , t)
  call BeStreamXtract(t, C_LOC(what), nbits, nitems)
end subroutine

! write data into a little endian stream
subroutine stream_put_le(self, what, nbits, nitems)
  implicit none
  class(le_stream), intent(INOUT) :: self
  integer(C_INT32_T), dimension(*), intent(IN), target  :: what
#if ! defined(__GFORTRAN__) && ! defined(__INTEL_COMPILER)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: what
!$PRAGMA IGNORE_TKR what
!DIR$ IGNORE_TKR what
!IBM* IGNORE_TKR what
#endif
  integer(C_INT32_T), intent(IN), value :: nbits, nitems
  type(C_PTR) :: t
  t = transfer( LOC(self%s) , t)
  call LeStreamInsert(t, C_LOC(what), nbits, nitems)
end subroutine

! write data into a big endian stream
subroutine stream_put_be(self, what, nbits, nitems)
  implicit none
  class(be_stream), intent(INOUT) :: self
  integer(C_INT32_T), dimension(*), intent(IN), target  :: what
#if ! defined(__GFORTRAN__) && ! defined(__INTEL_COMPILER)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: what
!$PRAGMA IGNORE_TKR what
!DIR$ IGNORE_TKR what
!IBM* IGNORE_TKR what
#endif
  integer(C_INT32_T), intent(IN), value :: nbits, nitems
  type(C_PTR) :: t
  t = transfer( LOC(self%s) , t)
  call BeStreamInsert(t, C_LOC(what), nbits, nitems)
end subroutine

subroutine pack32_put_init(accum, cnt)
  implicit none
  integer(kind=8), intent(INOUT) :: accum
  integer(kind=4), intent(INOUT) :: cnt
  accum = 0
  cnt = 64
end subroutine

subroutine pack32_put_fast(accum, cnt, token, nbits)
  implicit none
  integer(kind=8), intent(INOUT) :: accum
  integer(kind=4), intent(INOUT) :: cnt
  integer(kind=4), intent(IN)    :: token, nbits
  integer(kind=8) :: token64, mask64
  token64 = token
!   mask64 = maskr(nbits)
  mask64 = -1
  mask64 = rshift(mask64, 64-nbits)
  cnt    = cnt - nbits
  accum = ior( lshift(accum,nbits), iand(token64, mask64) )
end subroutine

subroutine pack32_put_check(sp, indx, accum, cnt)
  implicit none
  integer(kind=4), dimension(*), intent(INOUT) :: sp
  integer(kind=8), intent(INOUT) :: accum
  integer(kind=4), intent(INOUT) :: cnt, indx
  if(cnt < 32) then
    sp(indx) = accum
    accum = rshift(accum, 32)
    indx = indx + 1
    cnt = cnt + 32
  endif
end subroutine

subroutine pack32_put(sp, indx, accum, cnt, token, nbits)
  implicit none
  integer(kind=4), dimension(*), intent(INOUT) :: sp
  integer(kind=8), intent(INOUT) :: accum
  integer(kind=4), intent(INOUT) :: cnt, indx
  integer(kind=4), intent(IN) :: token, nbits
  call pack32_put_check(sp, indx, accum, cnt)
  call pack32_put_fast(accum, cnt, token, nbits)
end subroutine

subroutine pack32_put_flush(sp, indx, accum, cnt)
  implicit none
  integer(kind=4), dimension(*), intent(INOUT) :: sp
  integer(kind=8), intent(INOUT) :: accum
  integer(kind=4), intent(INOUT) :: cnt, indx
  if(cnt < 32) call pack32_put_check(sp, indx, accum, cnt)
  if(cnt < 64) then
    sp(indx) = lshift(accum, cnt - 32)
    indx = indx + 1
    cnt = 64
  endif
end subroutine

subroutine pack_stream(stream, indx, what, nbits, nitems)
  implicit none
  integer(kind=4), dimension(*), intent(INOUT) :: stream
  integer(kind=4), intent(INOUT) :: indx
  integer(kind=4), dimension(*), intent(IN) :: what
  integer(kind=4), intent(IN) :: nbits, nitems
  integer :: i, cnt
  integer(kind=8) :: accum
  indx = 0
  call pack32_put_init(accum, cnt)
  do i = 1, abs(nitems)
    call pack32_put(stream, indx, accum, cnt, what(i), nbits)
  enddo
  call pack32_put_flush(stream, indx, accum, cnt)
end subroutine
end module

#define NPTS 4097
program test_fortran_streams
  use rmn_stream_pack
  implicit none
  interface
    function c_elapsed_us() result(t) bind(C,name='elapsed_us')  ! elapsed microseconds
      import C_INT64_T
      implicit none
      integer(C_INT64_T) :: t
    end function
    function c_elapsed_cycles() result(t) bind(C,name='elapsed_cycles')  ! elapsed timer ticks
      import C_INT64_T
      implicit none
      integer(C_INT64_T) :: t
    end function
    function c_cycles_counter_freq() result(t) bind(C,name='cycles_counter_freq')  ! timer tick frequency
      import C_INT64_T
      implicit none
      integer(C_INT64_T) :: t
    end function
    function c_cycles_to_ns(cycles) result(ns) bind(C, name='cycles_to_ns')
      import :: C_DOUBLE, C_INT64_T
      implicit none
      integer(C_INT64_T), intent(IN), value :: cycles
      real(C_DOUBLE) :: ns
    end function
end interface
  type(le_stream) :: lestream
  type(be_stream) :: bestream
  integer, dimension(1,NPTS,1) :: unpacked_data
  integer, dimension(1,NPTS,1) :: packed_data_le, packed_data_be
  integer, dimension(1,NPTS,1) :: restored_data_le, restored_data_be
  integer :: i, errorsle, errorsbe
  integer(C_INT64_T) :: t0, t1, t2, t3

  do i = 1, NPTS
    unpacked_data(1,i,1) = i
  enddo
  print *,'=== original data ==='
  print 1,unpacked_data(1,1:8,1)

  packed_data_le = -1
  packed_data_be = -1
  restored_data_le = 0
  restored_data_be = 0
  call lestream % create(packed_data_le)
  call bestream % create(packed_data_be)
  t0 = c_elapsed_cycles()
  call lestream % put(unpacked_data, 16, NPTS)
  t0 = c_elapsed_cycles() - t0
  t1 = c_elapsed_cycles()
  call bestream % put(unpacked_data, 16, NPTS)
  t1 = c_elapsed_cycles() - t1
  call lestream % flush()
  call bestream % flush()
  call lestream % rewind()
  call bestream % rewind()
  t2 = c_elapsed_cycles()
  call lestream % get(restored_data_le, 16, NPTS)
  t2 = c_elapsed_cycles() - t2
  t3 = c_elapsed_cycles()
  call bestream % get(restored_data_be, 16, NPTS)
  t3 = c_elapsed_cycles() - t3
  errorsle = 0
  errorsbe = 0
  do i = 1, NPTS
    if(restored_data_le(1,i,1) .ne. unpacked_data(1,i,1)) errorsle = errorsle +1
    if(restored_data_be(1,i,1) .ne. unpacked_data(1,i,1)) errorsbe = errorsbe +1
  enddo
  print *,'npts =', NPTS, 'errors le/be =', errorsle, errorsbe
  print *,'t =', t0, t1, t2, t3
  print *,'=== packed data (le/be) ==='
  print 1,packed_data_le(1,1:8,1)
  print 1,packed_data_be(1,1:8,1)
  print *,'=== restored data (le/be) ==='
  print 1,restored_data_le(1,1:8,1), restored_data_le(1,NPTS-6:NPTS,1)
  print 1,restored_data_be(1,1:8,1), restored_data_be(1,NPTS-6:NPTS,1)
1 format(8z9.8,3x,8z9.8)
end program
