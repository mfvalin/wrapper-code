/*
Copyright (C) 2022  Recherche en Prevision Numerique

This code is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation,
version 2.1 of the License.

This code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Library General Public License for more details.
*/
#if ! defined(IN_FORTRAN_CODE) && ! defined(__GFORTRAN__)

void average_2x1_I32(int32_t * restrict src,                           int32_t * restrict avg, uint32_t n) ;
void average_2x1_F32(float   * restrict src,                           float   * restrict avg, uint32_t n) ;
void average_2x2_I32(int32_t * restrict src1, int32_t * restrict src2, int32_t * restrict avg, uint32_t n) ;
void average_2x2_F32(float   * restrict src1, float   * restrict src2, float   * restrict avg, uint32_t n) ;
void average_2x2_2D_I32(int32_t * restrict src, int32_t * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj) ;
void average_2x2_2D_F32(float   * restrict src, float   * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj) ;

void avgres_2x2_2D_I32(int32_t * restrict src, int32_t * restrict avg, int32_t * restrict res, int ni, int lni, int nj) ;

void expand_2x1_row_along_x_I32(int32_t * restrict dst, int32_t * restrict avg, int n) ;
void expand_2x1_row_along_x_F32(float   * restrict dst, float   * restrict avg, int n) ;
void expand_2x2_2_rows_I32(int32_t * restrict row0, int32_t * restrict row1, int32_t * src0, int32_t * src1, int n) ;
void expand_2x2_2_rows_F32(float   * restrict row0, float   * restrict row1, float   * src0, float   * src1, int n) ;
void expand_2x2_2D_I32(int32_t * restrict dst, int32_t * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj) ;
void expand_2x2_2D_F32(float   * restrict dst, float   * restrict avg, uint32_t ni, uint32_t lni, uint32_t nj) ;

#else
  interface
    subroutine average_2x2_2D_I32(src, avg, ni, lni, nj) bind(C,name='average_2x2_2D_I32')
      import :: C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: ni, lni, nj
      integer(C_INT32_T), dimension(lni,nj), intent(IN) :: src
      integer(C_INT32_T), dimension(ni,nj), intent(OUT) :: avg
    end subroutine average_2x2_2D_I32
    subroutine average_2x2_2D_F32(src, avg, ni, lni, nj) bind(C,name='average_2x2_2D_F32')
      import :: C_INT32_T, C_FLOAT
      implicit none
      integer(C_INT32_T), intent(IN), value :: ni, lni, nj
      real(C_FLOAT), dimension(lni,nj), intent(IN) :: src
      real(C_FLOAT), dimension(ni,nj), intent(OUT) :: avg
    end subroutine average_2x2_2D_F32
    subroutine avgres_2x2_2D_I32(src, avg, res, ni, lni, nj) bind(C,name='avgres_2x2_2D_I32')
      import :: C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: ni, lni, nj
      integer(C_INT32_T), dimension(lni,nj), intent(IN) :: src
      integer(C_INT32_T), dimension((ni+1)/2,(nj+1)/2), intent(OUT) :: avg
      integer(C_INT32_T), dimension(lni,nj), intent(OUT) :: res
    end subroutine avgres_2x2_2D_I32
    subroutine avgres_2x2_2D_F32(src, avg, res, ni, lni, nj) bind(C,name='avgres_2x2_2D_F32')
      import :: C_INT32_T, C_FLOAT
      implicit none
      integer(C_INT32_T), intent(IN), value :: ni, lni, nj
      real(C_FLOAT), dimension(lni,nj), intent(IN) :: src
      real(C_FLOAT), dimension((ni+1)/2,(nj+1)/2), intent(OUT) :: avg
      real(C_FLOAT), dimension(lni,nj), intent(OUT) :: res
    end subroutine avgres_2x2_2D_F32
    subroutine expand_2x2_2D_I32(dst, avg, ni, lni, nj) bind(C,name='expand_2x2_2D_I32')
      import :: C_INT32_T
      implicit none
      integer(C_INT32_T), intent(IN), value :: ni, lni, nj
      integer(C_INT32_T), dimension(lni,nj), intent(OUT) :: dst
      integer(C_INT32_T), dimension(ni,nj), intent(IN) :: avg
    end subroutine expand_2x2_2D_I32
    subroutine expand_2x2_2D_F32(dst, avg, ni, lni, nj) bind(C,name='expand_2x2_2D_F32')
      import :: C_INT32_T, C_FLOAT
      implicit none
      integer(C_INT32_T), intent(IN), value :: ni, lni, nj
      real(C_FLOAT), dimension(lni,nj), intent(OUT) :: dst
      real(C_FLOAT), dimension(ni,nj), intent(IN) :: avg
    end subroutine expand_2x2_2D_F32
  end interface
#endif
