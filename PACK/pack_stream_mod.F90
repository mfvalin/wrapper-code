module pack_stream_mod
  use ISO_C_BINDING
  implicit none
#include "pack_stream.hf"
  public :: init, clone, create
  type, public :: pstream
    private
    type(c_pstream) :: c
  contains
  procedure :: init
  procedure :: clone
  procedure :: create
  procedure :: alloc
  GENERIC :: new => init, clone, create, alloc
  procedure :: unpack_u64, unpack_i64, unpack_u32, unpack_i32
  GENERIC :: unpack_u => unpack_u64, unpack_u32
  GENERIC :: unpack_i => unpack_i64, unpack_i32
  procedure :: pack_u64, pack_i64, pack_u32, pack_i32
  GENERIC :: pack_u => pack_u64, pack_u32
  GENERIC :: pack_i => pack_i64, pack_i32
  end type

  interface
  function pstream_unpack_u64(this, src, dest, n) result(ntoken) bind(C,name='pstream_unpack_u64')
    import :: C_PTR, C_INT, C_INT64_T
    implicit none
    type(C_PTR), intent(IN) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT64_T), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT) :: ntoken
  end function pstream_unpack_u64

  function pstream_unpack_i64(this, src, dest, n) result(ntoken) bind(C,name='pstream_unpack_i64')
    import :: C_PTR, C_INT, C_INT64_T
    implicit none
    type(C_PTR), intent(IN) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT64_T), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT) :: ntoken
  end function pstream_unpack_i64

  function pstream_unpack_u32(this, src, dest, n) result(ntoken) bind(C,name='pstream_unpack_u32')
    import :: C_PTR, C_INT
    implicit none
    type(C_PTR), intent(IN) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT) :: ntoken
  end function pstream_unpack_u32

  function pstream_unpack_i32(this, src, dest, n) result(ntoken) bind(C,name='pstream_unpack_i32')
    import :: C_PTR, C_INT
    implicit none
    type(C_PTR), intent(IN) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT) :: ntoken
  end function pstream_unpack_i32

  function pstream_pack_u64(this, src, dest, n, nbits, navail, mode) result(nused) bind(C,name='pstream_pack_u64')
    import :: C_PTR, C_INT, C_INT64_T
    implicit none
    type(C_PTR), intent(IN) :: this
    integer(C_INT64_T), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT), intent(IN), value :: nbits
    integer(C_INT), intent(IN), value :: navail, mode
    integer(C_INT) :: nused
  end function pstream_pack_u64

  function pstream_pack_i64(this, src, dest, n, nbits, navail, mode) result(nused) bind(C,name='pstream_pack_i64')
    import :: C_PTR, C_INT, C_INT64_T
    implicit none
    type(C_PTR), intent(IN) :: this
    integer(C_INT64_T), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT), intent(IN), value :: nbits
    integer(C_INT), intent(IN), value :: navail, mode
    integer(C_INT) :: nused
  end function pstream_pack_i64

  function pstream_pack_u32(this, src, dest, n, nbits, navail, mode) result(nused) bind(C,name='pstream_pack_u32')
    import :: C_PTR, C_INT, C_INT64_T
    implicit none
    type(C_PTR), intent(IN) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT), intent(IN), value :: nbits
    integer(C_INT), intent(IN), value :: navail, mode
    integer(C_INT) :: nused
  end function pstream_pack_u32

  function pstream_pack_i32(this, src, dest, n, nbits, navail, mode) result(nused) bind(C,name='pstream_pack_i32')
    import :: C_PTR, C_INT, C_INT64_T
    implicit none
    type(C_PTR), intent(IN) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT), intent(IN), value :: nbits
    integer(C_INT), intent(IN), value :: navail, mode
    integer(C_INT) :: nused
  end function pstream_pack_i32

  end interface

  contains

  subroutine init(this)
    implicit none
    class(pstream), intent(INOUT) :: this
    this%c%a  = 0
    this%c%s  = C_NULL_PTR
    this%c%p  = C_NULL_PTR
    this%c%m  = 0
    this%c%nf = 0
    this%c%na = 0
    this%c%ni = 0
    this%c%nb = 0
  end subroutine init

  subroutine create(this, addr, size)
    implicit none
    class(pstream), intent(INOUT) :: this
    type(C_PTR), intent(IN), value :: addr
    integer(C_INT), intent(IN), value :: size
    call init(this)
    this%c%s  = addr
    this%c%ni = size
  end subroutine create

  subroutine alloc(this, size)
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT), intent(IN), value :: size
    integer(C_SIZE_T) :: szt
    interface
    function malloc_c(sz) result(p) bind(C,name='malloc')
      import :: C_PTR, C_SIZE_T
      implicit none
      integer(C_SIZE_T), intent(IN), value :: sz
      type(C_PTR) :: p
    end function malloc_c
    end interface
    call init(this)
    szt = size
    this%c%s  = malloc_c(szt)
    this%c%ni = size
  end subroutine alloc

  subroutine clone(this, src)
    implicit none
    class(pstream), intent(INOUT) :: this
    type(pstream), intent(IN) :: src
    this%c%a  = src%c%a
    this%c%s  = src%c%s
    this%c%p  = src%c%p
    this%c%m  = src%c%m
    this%c%nf = src%c%nf
    this%c%na = src%c%na
    this%c%ni = src%c%ni
    this%c%nb = src%c%nb
  end subroutine clone

  function unpack_u64(this, src, dest, n) result(ntoken)  ! unpack into 64 bit unsigned integers
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT64_T), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT) :: ntoken

    ntoken = pstream_unpack_u64(this % c % p, src, dest, n)
  end function unpack_u64

  function unpack_i64(this, src, dest, n) result(ntoken)  ! unpack into 64 bit signed integers
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT64_T), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT) :: ntoken

    ntoken = pstream_unpack_i64(this % c % p, src, dest, n)
  end function unpack_i64

  function unpack_u32(this, src, dest, n) result(ntoken)  ! unpack into 32 bit unsigned integers
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT) :: ntoken

    ntoken = pstream_unpack_u32(this % c % p, src, dest, n)
  end function unpack_u32

  function unpack_i32(this, src, dest, n) result(ntoken)  ! unpack into 32 bit signed integers
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT) :: ntoken

    ntoken = pstream_unpack_i32(this % c % p, src, dest, n)
  end function unpack_i32

  function pack_u64(this, src, dest, n, nbits, navail, mode) result(nused)
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT64_T), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT), intent(IN), value :: nbits
    integer(C_INT), intent(IN), value :: navail, mode
    integer(C_INT) :: nused

    nused = pstream_pack_u64(this %c % p, src, dest, n, nbits, navail, mode)
  end function pack_u64

  function pack_i64(this, src, dest, n, nbits, navail, mode) result(nused)
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT64_T), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT), intent(IN), value :: nbits
    integer(C_INT), intent(IN), value :: navail, mode
    integer(C_INT) :: nused

    nused = pstream_pack_i64(this %c % p, src, dest, n, nbits, navail, mode)
  end function pack_i64

  function pack_u32(this, src, dest, n, nbits, navail, mode) result(nused)
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT), intent(IN), value :: nbits
    integer(C_INT), intent(IN), value :: navail, mode
    integer(C_INT) :: nused

    nused = pstream_pack_u32(this %c % p, src, dest, n, nbits, navail, mode)
  end function pack_u32

  function pack_i32(this, src, dest, n, nbits, navail, mode) result(nused)
    implicit none
    class(pstream), intent(INOUT) :: this
    integer(C_INT), dimension(*), intent(IN) :: src
    integer(C_INT), dimension(*), intent(OUT) :: dest
    integer(C_INT), intent(IN), value :: n
    integer(C_INT), intent(IN), value :: nbits
    integer(C_INT), intent(IN), value :: navail, mode
    integer(C_INT) :: nused

    nused = pstream_pack_i32(this %c % p, src, dest, n, nbits, navail, mode)
  end function pack_i32

end module pack_stream_mod
