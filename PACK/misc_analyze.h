#if defined(IN_FORTRAN_CODE) || defined(__GFORTRAN__)

  interface
    subroutine AnalyzeCompressionErrors(fa, fb, np, small, str) bind(C, name='AnalyzeCompressionErrors')
      import :: C_FLOAT, C_INT32_T, C_CHAR
      implicit none
      real(C_FLOAT), dimension(*), intent(IN) :: fa, fb
#if ! defined(__GFORTRAN__) && ! defined(__INTEL_COMPILER)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: fa, fb
!$PRAGMA IGNORE_TKR fa, fb
!DIR$ IGNORE_TKR fa, fb
!IBM* IGNORE_TKR fa, fb
#endif
      integer(C_INT32_t), intent(IN), value :: np
      real(C_FLOAT), intent(IN), value :: small
      character(C_CHAR), dimension(*), intent(IN) :: str
    end subroutine

    subroutine FloatTestData_2d(z, ni, lni, nj, alpha, delta, noise, dkind, minf, maxf) bind(C, name='FloatTestData_2d')
      import :: C_FLOAT, C_INT32_T
      implicit none
      real(C_FLOAT), dimension(*), intent(OUT) :: z
#if ! defined(__GFORTRAN__) && ! defined(__INTEL_COMPILER)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: z
!$PRAGMA IGNORE_TKR z
!DIR$ IGNORE_TKR z
!IBM* IGNORE_TKR z
#endif
      integer(C_INT32_t), intent(IN), value :: ni, lni, nj, dkind
      real(C_FLOAT), intent(IN), value :: alpha, delta, noise
      real(C_FLOAT), intent(OUT) :: minf, maxf
    end subroutine
  end interface

#else

void AnalyzeCompressionErrors(float *fa, float *fb, int np, float small, char *str);
void FloatTestData_2d(float *z, int ni, int lni, int nj, float alpha, float delta, float noise, int kind, float *min, float *max);

#endif
