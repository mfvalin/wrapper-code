module lorenzo_mod
  use ISO_C_BINDING
  implicit none
! C interfaces
interface
  ! the C version using SIMD is faster that the Fortran version
  subroutine lorenzopredict(orig, diff, ni, lnio, lnid, nj) bind(C,name='LorenzoPredict_c')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
    integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
    integer(C_INT32_T), dimension(lnid,nj), intent(OUT) :: diff
  end subroutine lorenzopredict
  subroutine lorenzopredict_c(orig, diff, ni, lnio, lnid, nj) bind(C,name='LorenzoPredict_c')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
    integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
    integer(C_INT32_T), dimension(lnid,nj), intent(OUT) :: diff
  end subroutine lorenzopredict_c

  ! the C version using SIMD is faster that the Fortran version
  subroutine lorenzopredictinplace(orig, ni, lnio, nj) bind(C,name='LorenzoPredictInplace_c')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: ni, lnio, nj
    integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
  end subroutine lorenzopredictinplace
  subroutine lorenzopredictinplace_c(orig, ni, lnio, nj) bind(C,name='LorenzoPredictInplace_c')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: ni, lnio, nj
    integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
  end subroutine lorenzopredictinplace_c

  ! the C version seems slower that the Fortran version
  subroutine lorenzounpredict_c(orig, diff, ni, lnio, lnid, nj) bind(C,name='LorenzoUnpredict_c')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
    integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
    integer(C_INT32_T), dimension(lnid,nj), intent(OUT) :: diff
  end subroutine lorenzounpredict_c

  ! the C version seems slower that the Fortran version
  subroutine lorenzounpredictinplace_c(orig, ni, lnio, nj) bind(C,name='LorenzoUnpredictInplace_c')
    import :: C_INT32_T
    implicit none
    integer(C_INT32_T), intent(IN), value :: ni, lnio, nj
    integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
  end subroutine lorenzounpredictinplace_c
end interface

contains

! 3 point Lorenzo predictor (forward derivative)
! upon exit, diff contains original value - predicted value
! with explicit SIMD primitives, the C version (LorenzoPredict_c) outperforms the Fortran version
subroutine lorenzopredict_f(orig, diff, ni, lnio, lnid, nj) bind(C,name='LorenzoPredict_f')
  implicit none
  integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
  integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
  integer(C_INT32_T), dimension(lnid,nj), intent(OUT) :: diff
  integer :: i, j

  diff(1,1) = orig(1,1)
  do i = 2, ni
    diff(i,1) = orig(i,1) - orig(i-1,1)         ! botom row, 1D prediction along i
  enddo
  do j = 2 , nj                                 ! all rows except bottom row
    diff(1,j) = orig(1,j) - orig(1,j-1)         ! special case for leftmost column, 1D prediction along j
    do i = 2, ni
      diff(i,j) = orig(i,j) - ( orig(i-1,j) + orig(i,j-1) - orig(i-1,j-1) )  ! value - predicted value
    enddo
  enddo
end subroutine lorenzopredict_f

#if 0
! 5 point Lorenzo predictor (centered derivatives)
subroutine lorenzopredict2_f(orig, diff, ni, lnio, lnid, nj) bind(C,name='LorenzoPredict2_f')
  implicit none
  integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
  integer(C_INT32_T), dimension(lnio,nj), intent(IN)  :: orig
  integer(C_INT32_T), dimension(lnid,nj), intent(OUT) :: diff
  integer :: i, j

  diff(1,1) = orig(1,1)                ! bottom row no prediction for first point
  diff(2,1) = orig(2,1) - orig(1,1)    ! bottom row 1 point prediction for second point
  diff(1,2) = orig(1,2) - orig(1,1)    ! first point of second row, 1 point prediction
  diff(2,2) = orig(2,2) - ( orig(1,2) + orig(2,1) - orig(1,1) )  ! second point row 2, 3 point prediction
  do i = 3, ni
    diff(i,1) = orig(i,1) - (orig(i-1,1)*2 - orig(i-2,1))                   ! bottom row, 2 point prediction
    diff(i,2) = orig(i,2) - ( orig(i-1,2) + (orig(i,1) - orig(i-2,1))/2 )   ! second row, 3 point prediction
  enddo
  do j = 3 , nj                                 ! all rows except 2 bottom rows
    diff(1,j) = orig(1,j) - (orig(1,j-1)*2 - orig(1,j-2) )         ! first column, 2 point prediction along j
    diff(2,j) = orig(2,j) - ( orig(2,j-1) + (orig(1,j) - orig(1,j-2))/2 )  ! second column, 3 point prediction
    do i = 3, ni
      diff(i,j) = orig(i,j) - &
                ( orig(i-1,j-1) + ( orig(i,j-1) - orig(i-2,j-1) + orig(i-1,j) - orig(i-1,j-2) )/2 )  ! 5 point prediction
    enddo
  enddo
end subroutine lorenzopredict2_f
#endif

! Lorenzo restoration
! upon exit, orig contains values restored from prediction in diff
! the Fortran version seems to have a performance >= C version
subroutine lorenzounpredict_f(orig, diff, ni, lnio, lnid, nj) bind(C,name='LorenzoUnpredict_f')
  implicit none
  integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
  integer(C_INT32_T), dimension(lnio,nj), intent(OUT) :: orig
  integer(C_INT32_T), dimension(lnid,nj), intent(IN)  :: diff
  integer :: i, j

  orig(1,1) = diff(1,1)
  do i = 2, ni
    orig(i,1) = diff(i,1) + orig(i-1,1)         ! botom row, 1D prediction along i
  enddo
  do j = 2 , nj                                 ! all rows except bottom row
    orig(1,j) = diff(1,j) + orig(1,j-1)         ! special case for leftmost column, 1D prediction along j
    do i = 2, ni
      orig(i,j) = diff(i,j) + ( orig(i-1,j) + orig(i,j-1) - orig(i-1,j-1) ) ! value + predicted value
    enddo
  enddo
end subroutine lorenzounpredict_f

! the Fortran version seems to have a performance >= C version
subroutine lorenzounpredict(orig, diff, ni, lnio, lnid, nj)
  implicit none
  integer(C_INT32_T), intent(IN), value :: ni, lnio, lnid, nj
  integer(C_INT32_T), dimension(lnio,nj), intent(OUT) :: orig
  integer(C_INT32_T), dimension(lnid,nj), intent(IN)  :: diff
  call lorenzounpredict_f(orig, diff, ni, lnio, lnid, nj)
end subroutine lorenzounpredict

! in place Lorenzo prediction
! upon exit, f contains original value - predicted alue
! with explicit SIMD primitives, the C version (LorenzoPredictInplace_c) outperforms the Fortran version
subroutine lorenzopredictinplace_f(f, ni, lni, nj) bind(C,name='LorenzoPredictInplace_f')
  implicit none
  integer(C_INT32_T), intent(IN), value :: ni, lni, nj
  integer(C_INT32_T), dimension(lni,nj), intent(INOUT) :: f
  integer :: i, j

  do j = nj, 2, -1                                            ! all rows except bottom row
    do i = ni, 2, -1
      f(i,j) = f(i,j) - ( f(i-1,j) + f(i,j-1) - f(i-1,j-1) )  ! value - predicted value
    enddo
    f(1,j) = f(1,j) - f(1,j-1)                                ! special case for leftmost column, 1D prediction along j
  enddo
  do i = ni, 2, -1
    f(i,1) = f(i,1) - f(i-1,1)                                ! botom row, 1D prediction along i
  enddo
end subroutine lorenzopredictinplace_f

! in place Lorenzo restore
! upon entry, f contains original value - predicted alue
! upon exit, f contains restored original values
! the Fortran version seems to have a performance >= C version
subroutine lorenzounpredictinplace_f(f, ni, lni, nj) bind(C,name='LorenzoUnpredictInplace_f')
  implicit none
  integer(C_INT32_T), intent(IN), value :: ni, lni, nj
  integer(C_INT32_T), dimension(lni,nj), intent(INOUT) :: f
  integer :: i, j

  do i = 2, ni                                                ! bottom row, 1D prediction along i
    f(i,1) = f(i,1) + f(i-1,1)                                ! (original - precicted) + predicted
  enddo
  do j = 2, nj
    f(1,j) = f(1,j) + f(1,j-1)                                ! special case for leftmost column, 1D prediction along j
    do i = 2,ni
      f(i,j) = f(i,j)+ ( f(i-1,j) + f(i,j-1) - f(i-1,j-1) )   ! (original - precicted) + predicted
    enddo
  enddo
end subroutine lorenzounpredictinplace_f

! the Fortran version seems to have a performance >= C version
subroutine lorenzounpredictinplace(f, ni, lni, nj)
  implicit none
  integer(C_INT32_T), intent(IN), value :: ni, lni, nj
  integer(C_INT32_T), dimension(lni,nj), intent(INOUT) :: f
  call lorenzounpredictinplace_f(f, ni, lni, nj)
end subroutine lorenzounpredictinplace

end module lorenzo_mod
