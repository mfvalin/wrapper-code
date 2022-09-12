/*
 *
 * This software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation,
 * version 2.1 of the License.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 */

#if ! defined(IN_FORTRAN_CODE)

void LorenzoPredict_c(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);
void LorenzoUnpredict_c(int32_t *orig, int32_t *diff, int ni, int lnio, int lnid, int nj);

#else

interface
  subroutine lorenzopredict(orig, diff, ni, lnio, lnid, nj) bind(C,name='LorenzoPredict_c')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
    integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
    integer(C_INT32_T), dimension(lnid,nj), intent(OUT) :: diff
  end subroutine lorenzopredict_c
  subroutine lorenzounpredict(orig, diff, ni, lnio, lnid, nj) bind(C,name='LorenzoUnpredict_c')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
    integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
    integer(C_INT32_T), dimension(lnid,nj), intent(OUT) :: diff
  end subroutine lorenzounpredict_c
end interface

#endif
