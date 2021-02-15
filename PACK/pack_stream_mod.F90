module pack_stream_mod
  use ISO_C_BINDING
  implicit none
#include "pack_stream.hf"
  type :: pstream
    private
    type(c_pstream) :: c
  contains
  procedure init
  procedure clone
  procedure create
  GENERIC :: new => init, clone, create
  end type

  contains

  subroutine init(stream)
    implicit none
    class(pstream), intent(INOUT) :: stream
    stream%c%a  = 0
    stream%c%s  = C_NULL_PTR
    stream%c%p  = C_NULL_PTR
    stream%c%m  = 0
    stream%c%nf = 0
    stream%c%na = 0
    stream%c%ni = 0
    stream%c%nb = 0
  end subroutine init

  subroutine create(stream, addr, size)
    implicit none
    class(pstream), intent(INOUT) :: stream
    type(C_PTR), intent(IN), value :: addr
    integer(C_INT), intent(IN), value :: size
    call init(stream)
    stream%c%s  = addr
    stream%c%ni = size
  end subroutine create

  subroutine clone(stream, src)
    implicit none
    class(pstream), intent(INOUT) :: stream
    type(pstream), intent(IN) :: src
    stream%c%a  = src%c%a
    stream%c%s  = src%c%s
    stream%c%p  = src%c%p
    stream%c%m  = src%c%m
    stream%c%nf = src%c%nf
    stream%c%na = src%c%na
    stream%c%ni = src%c%ni
    stream%c%nb = src%c%nb
  end subroutine clone

end module pack_stream_mod
