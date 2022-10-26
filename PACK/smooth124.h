#if  ! defined(IN_FORTRAN_CODE) && ! defined(__GFORTRAN__)

#include <stdint.h>

void Ismooth_124_inplace(int32_t *f, int n) ;
void Fsmooth_124_inplace(float *f, int n) ;
void Fsmooth_124(float *s, float *d, int n) ;
void Fsmooth_124_2D_inplace(float *f, int ni, int lni, int nj) ;
void Fsmooth_124_2D(float *s, float *d, int ni, int lnis, int lnid, int nj) ;

#else
interface
  subroutine Ismooth_124_inplace(f, n) bind(C,name='Ismooth_124_inplace')
    import :: C_INT32_t
    implicit none
    integer(C_INT32_t), intent(IN), value :: n
    integer(C_INT32_t), dimension(*), intent(INOUT) :: f
  end subroutine Ismooth_124_inplace
  subroutine Fsmooth_124_inplace(f, n) bind(C,name='Fsmooth_124_inplace')
    import :: C_FLOAT, C_INT32_t
    implicit none
    integer(C_INT32_t), intent(IN), value :: n
    real(C_FLOAT), dimension(*), intent(INOUT) :: f
  end subroutine Fsmooth_124_inplace
  subroutine Fsmooth_124(s, d, n) bind(C,name='Fsmooth_124')
    import :: C_FLOAT, C_INT32_t
    implicit none
    integer(C_INT32_t), intent(IN), value :: n
    real(C_FLOAT), dimension(*), intent(IN) :: s
    real(C_FLOAT), dimension(*), intent(OUT) :: d
  end subroutine Fsmooth_124
  subroutine Fsmooth_124_2D_inplace(f, ni, lni, nj) bind(C,name='Fsmooth_124_2D_inplace')
    import :: C_FLOAT, C_INT32_t
    implicit none
    integer(C_INT32_t), intent(IN), value :: ni, lni, nj
    real(C_FLOAT), dimension(*), intent(INOUT) :: f
  end subroutine Fsmooth_124_2D_inplace
  subroutine Fsmooth_124_2D(s, d, ni, lnis, lnid, nj) bind(C,name='Fsmooth_124_2D')
    import :: C_FLOAT, C_INT32_t
    implicit none
    integer(C_INT32_t), intent(IN), value :: ni, lnis, lnid, nj
    real(C_FLOAT), dimension(lnis,nj), intent(IN) :: s
    real(C_FLOAT), dimension(lnis,nj), intent(OUT) :: d
  end subroutine Fsmooth_124_2D
end interface
#endif
