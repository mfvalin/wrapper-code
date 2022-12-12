/*
  Copyright (C) 2022  Environnement Canada

  This is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation,
  version 2.1 of the License.

  This software is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.
*/

#if defined(IN_FORTRAN_CODE)  || defined(__GFORTRAN__)
interface
  function ZfpCompress_set_debug(flag) result(old) bind(C, name='ZfpCompress_set_debug')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: flag
    integer(C_INT32_T) :: old
  end function
  function ZfpCompress_set_diag(flag) result(old) bind(C, name='ZfpCompress_set_diag')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: flag
    integer(C_INT32_T) :: old
  end function
  function get_zfp_codec_version() result(version) bind(C, name='get_zfp_codec_version')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T) :: version
  end function
  function get_ZFP_CODEC() result(version) bind(C, name='get_ZFP_CODEC')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T) :: version
  end function
  function zfp_codec_consistent() result(version) bind(C, name='zfp_codec_consistent')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T) :: version
  end function
  function ZfpCompressReal(array, nx, ny, nz, maxerror, precision, Ssize) result(Stream) bind(C, name='ZfpCompress')
    import :: C_INT32_T, C_FLOAT, C_PTR
    implicit none
    real(C_FLOAT), dimension(*), intent(IN) :: array
    integer(C_INT32_T), intent(IN), value :: nx, ny, nz, precision
    real(C_FLOAT), intent(IN), value :: maxerror
    integer(C_INT32_T), intent(OUT) :: Ssize
    type(C_PTR) :: Stream
  end function
  function ZfpExpandReal(array, nx, ny, nz, Stream, Ssize) result(status) bind(C, name='ZfpExpand')
    import :: C_INT32_T, C_FLOAT, C_PTR
    implicit none
    real(C_FLOAT), dimension(*), intent(OUT) :: array
    integer(C_INT32_T), intent(IN), value :: nx, ny, nz
    type(C_PTR), intent(IN), value :: Stream
    integer(C_INT32_T), intent(OUT) :: Ssize
    integer(C_INT32_T) :: status
  end function
  function ZfpArrayDims(Stream, d, s, Ssize) result(ndims) bind(C, name='ZfpExpand')
    import :: C_INT32_T, C_PTR
    implicit none
    type(C_PTR), intent(IN), value :: Stream
    integer(C_INT32_T), intent(OUT) :: d, s
    integer(C_INT32_T), intent(IN), value :: Ssize
    integer(C_INT32_T) :: ndims
  end function
end interface
#else

int ZfpCompress_set_debug(int flag);
int ZfpCompress_set_diag(int flag);
uint32_t get_zfp_codec_version();
uint32_t get_ZFP_CODEC();
int32_t zfp_codec_consistent();
void *ZfpCompress(void* array, int nx, int ny, int nz, float maxerror, int precision , int *Ssize);
int32_t ZfpExpand(void* array, int nx, int ny, int nz, void *Stream, int Ssize);
int32_t ZfpArrayDims(void *Stream, int *d, int *s, int Ssize);

#endif
