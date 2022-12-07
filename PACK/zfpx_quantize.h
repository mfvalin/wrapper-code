/*
  Hopefully useful code
  Copyright (C) 2022  Recherche en Prevision Numerique

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * Neither the name of the copyright holder nor the names of its
    contributors may be used to endorse or promote products derived from
    this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  transforms and helper functions for use in floating point and integer compression
  interfaces to the zfpx functions

 */

#if defined IN_FORTRAN_CODE || defined(__GFORTRAN__)
interface
  subroutine zfpx_shuffle_2d(src, dst) bind(C, name='zfpx_shuffle_2d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(IN)  :: src
    integer(C_INT32_T), dimension(*), intent(OUT) :: dst
  end subroutine
  subroutine zfpx_unshuffle_2d(src, dst) bind(C, name='zfpx_unshuffle_2d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(IN)  :: src
    integer(C_INT32_T), dimension(*), intent(OUT) :: dst
  end subroutine
  subroutine zfpx_shuffle_3d(src, dst) bind(C, name='zfpx_shuffle_3d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(IN)  :: src
    integer(C_INT32_T), dimension(*), intent(OUT) :: dst
  end subroutine
  subroutine zfpx_unshuffle_3d(src, dst) bind(C, name='zfpx_unshuffle_3d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(IN)  :: src
    integer(C_INT32_T), dimension(*), intent(OUT) :: dst
  end subroutine
  subroutine zfpx_fwd_xform_1d(t) bind(C, name='zfpx_fwd_xform_1d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(INOUT)  :: t
  end subroutine
  subroutine zfpx_fwd_xform_2d(t) bind(C, name='zfpx_fwd_xform_2d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(INOUT)  :: t
  end subroutine
  subroutine zfpx_fwd_xform_3d(t) bind(C, name='zfpx_fwd_xform_3d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(INOUT)  :: t
  end subroutine
  subroutine zfpx_inv_xform_1d(t) bind(C, name='zfpx_inv_xform_1d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(INOUT)  :: t
  end subroutine
  subroutine zfpx_inv_xform_2d(t) bind(C, name='zfpx_inv_xform_2d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(INOUT)  :: t
  end subroutine
  subroutine zfpx_inv_xform_3d(t) bind(C, name='zfpx_inv_xform_3d')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(INOUT)  :: t
  end subroutine
  function int_to_negabinary(i) result(u) bind(C, name='int_to_negabinary')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: i
    integer(C_INT32_T) :: u
  end function
  subroutine v_int_to_negabinary(x, n) bind(C, name='v_int_to_negabinary')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(INOUT)  :: x
    integer(C_INT32_T), intent(IN), value :: n
  end subroutine
  function negabinary_to_int(u) result(i) bind(C, name='negabinary_to_int')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: u
    integer(C_INT32_T) :: i
  end function
  subroutine v_negabinary_to_int(x, n) bind(C, name='v_negabinary_to_int')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(INOUT)  :: x
    integer(C_INT32_T), intent(IN), value :: n
  end subroutine
  function get_ieee_emax(f, n) result(emax) bind(C, name='get_ieee_emax')
    import :: C_FLOAT, C_INT32_T
    implicit none
    real(C_FLOAT), dimension(*), intent(IN) :: f
    integer(C_INT32_T), intent(IN), value :: n
    integer(C_INT32_T) :: emax
  end function
  function zfpx_quantize(f, fi, n, emax) result(factor) bind(C, name='zfpx_quantize')
    import :: C_FLOAT, C_INT32_T
    implicit none
    real(C_FLOAT), dimension(*), intent(IN) :: f
    integer(C_INT32_T), dimension(*), intent(OUT)  :: fi
    integer(C_INT32_T), intent(IN), value :: n, emax
    real(C_FLOAT) :: factor
  end function
  subroutine zfpx_bit_plane_32_16(src, planes) bind(C, name='zfpx_bit_plane_32_16')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(IN)   :: src
    integer(C_INT32_T), dimension(*), intent(OUT)  :: planes
  end subroutine
  subroutine zfpx_bit_plane_32_64(src, planes) bind(C, name='zfpx_bit_plane_32_64')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(IN)   :: src
    integer(C_INT32_T), dimension(*), intent(OUT)  :: planes
  end subroutine
  subroutine zfpx_gather_64_64(f, lni, blocks, transform) bind(C, name='zfpx_gather_64_64')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(IN)   :: f
    integer(C_INT32_T), dimension(*), intent(OUT)  :: blocks
    integer(C_INT32_T), intent(IN), value :: lni, transform
  end subroutine
  subroutine zfpx_scatter_64_64(f, lni, blocks, transform) bind(C, name='zfpx_scatter_64_64')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), dimension(*), intent(OUT)  :: f
    integer(C_INT32_T), dimension(*), intent(IN)   :: blocks
    integer(C_INT32_T), intent(IN), value :: lni, transform
  end subroutine
end interface

#else

#if ! defined(ZFPX_QUANTIZE_EXTRAS)
#define ZFPX_QUANTIZE_EXTRAS

// reorder after 2 dimensional forward transform, explicit "no gather" version
void zfpx_shuffle_2d(int32_t *src, int32_t *dst);
// reorder before 2 dimensional inverse transform, explicit "no scatter" version
void zfpx_unshuffle_2d(int32_t *src, int32_t *dst);
// reorder after 3 dimensional forward transform, explicit "no gather" version
void zfpx_shuffle_3d(int32_t *src, int32_t *dst);
// reorder before 3 dimensional inverse transform, explicit "no scatter" version
void zfpx_unshuffle_3d(int32_t *src, int32_t *dst);

// reversible forward lifting transform, in place, 1 dimensional data (4)
void zfpx_fwd_xform_1d(int32_t *t);
// reversible forward lifting transform, in place, 2 dimensional data (4,4)
void zfpx_fwd_xform_2d(int32_t *t);
// reversible forward lifting transform, in place, 3 dimensional data (4,4,4)
void zfpx_fwd_xform_3d(int32_t *t);

// reversible inverse lifting transform, in place, 1 dimensional data (4)
void zfpx_inv_xform_1d(int32_t *t);
// reversible inverse lifting transform, in place, 2 dimensional data (4,4)
void zfpx_inv_xform_2d(int32_t *t0);
// reversible inverse lifting transform, in place, 3 dimensional data (4,4,4)
void zfpx_inv_xform_3d(int32_t *t0);

// transform a 64 x 64 block into 64 4x4x4 slices
void zfpx_gather_64_64(int32_t *f, int32_t lni, int32_t *blocks, int transform);
// inverse of zfpx_gather_64_64
void zfpx_scatter_64_64(int32_t *f, int32_t lni, int32_t *blocks, int transform);

// signed integer (2's complement) to negabinary (base -2) conversion
uint32_t int_to_negabinary(int32_t x);
void v_int_to_negabinary(int32_t x, int32_t n);
// negabinary (base -2) to signed integer (2's complement) conversion
int32_t negabinary_to_int(uint32_t x);
void v_negabinary_to_int(uint32_t x, int32_t n);

// get the IEEE exponent of the largest float (absolute value) in float array f
uint32_t get_ieee_emax(float *f, int n);

// quantize float array f by normalizing all exponents to the largest exponent emax
float zfpx_quantize(float *f, int32_t *fi, int n, uint32_t emax);

// 16 32 bit integers -> 32 bit planes (lower 16 bits in 64)
void zfpx_bit_plane_32_16(uint32_t *src, uint64_t *planes);
// 64 32 bit integers -> 32 bit planes (64 bits)
void zfpx_bit_plane_32_64(uint32_t *src, uint64_t *planes);

#endif

#endif
