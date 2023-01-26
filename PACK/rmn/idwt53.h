
#if ! defined(IN_FORTRAN_CODE) && ! defined(__GFORTRAN__)
void FDWT53i_1D_inplace_full(int32_t *x, uint32_t n);
void FDWT53i_1D_split_full(int32_t *x, int32_t *e, int32_t *o, uint32_t n);
void FDWT53i_1D_inplace_split_full(int32_t *x, uint32_t n);
void IDWT53i_1D_inplace_full(int32_t *x, uint32_t n);
void IDWT53i_1D_split_full(int32_t *x, int32_t *e, int32_t *o, uint32_t n);
void IDWT53i_1D_inplace_split_full(int32_t *x, uint32_t n);
void IDWT53i_quadrants(int32_t *z, int ni, int lni, int nj, int32_t *min, int32_t *max);
void FDWT53i_2D_split_inplace(int32_t *x, int ni, int lni, int nj);
void FDWT53i_2D_split_inplace_n(int32_t *x, int ni, int lni, int nj, int levels);
void IDWT53i_2D_split_inplace(int32_t *x, int ni, int lni, int nj);
void IDWT53i_2D_split_inplace_n(int32_t *x, int ni, int lni, int nj, int levels);
void FDWT53i_8x8_3level(int32_t *x, int lni);
void IDWT53i_8x8_3level(int32_t *x, int lni);
#else
interface
  subroutine FDWT53i_2D_split_inplace_n(x, ni, lni, nj, levels) bind(C, name='FDWT53i_2D_split_inplace_n')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_t), intent(IN), value :: ni, lni, nj, levels
    integer(C_INT32_t), dimension(lni,nj), intent(INOUT) :: x
  end subroutine FDWT53i_2D_split_inplace_n
end interface
#endif
