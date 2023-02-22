! Copyright (C) 2022  Environnement et Changement climatique Canada
! 
! This is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation,
! version 2.1 of the License.
! 
! This software is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! Author:
!     M. Valin,   Environnement et Changement climatique Canada, 2022
!
! interface to fine grained timer functions (time base counter)
! interface to wall clock (time of day)
!
module rmn_timers
  use ISO_C_BINDING
  implicit none

  private :: timer_t
  integer(C_INT64_t)            :: timer_t
!$OMP threadprivate(timer_t)

  interface
    function elapsed_us() result(t) bind(C,name='elapsed_us')  ! elapsed microseconds
      import C_INT64_T
      implicit none
      integer(C_INT64_T) :: t
    end function
    function elapsed_cycles() result(t) bind(C,name='elapsed_cycles')  ! elapsed timer ticks
      import C_INT64_T
      implicit none
      integer(C_INT64_T) :: t
    end function
    function elapsed_cycles_raw() result(t) bind(C,name='elapsed_cycles_raw')  ! elapsed timer ticks, no fencing
      import C_INT64_T
      implicit none
      integer(C_INT64_T) :: t
    end function
    function cycles_counter_freq() result(t) bind(C,name='cycles_counter_freq')  ! timer tick frequency
      import C_INT64_T
      implicit none
      integer(C_INT64_T) :: t
    end function
    function cycles_to_ns(cycles) result(ns) bind(C, name='cycles_to_ns')
      import :: C_DOUBLE, C_INT64_T
      implicit none
      integer(C_INT64_T), intent(IN), value :: cycles
      real(C_DOUBLE) :: ns
    end function
end interface
contains
  subroutine timer_start                       ! start timer
    implicit none
    timer_t = elapsed_cycles()
  end subroutine
  subroutine timer_stop                        ! stop timer
    implicit none
    timer_t = elapsed_cycles() - timer_t
  end subroutine
  function timer_get_cycles() result(cycles)   ! read timer (in cycles)
    implicit none
    integer(C_INT64_T) :: cycles
    cycles = timer_t
  end function
  function timer_get_ns() result(ns)           ! read timer (in ns)
    implicit none
    real(C_DOUBLE) :: ns
    ns = cycles_to_ns(timer_t)
  end function
  function timer_cycles() result(cycles)       ! read time base counter in ticks
    implicit none
    integer(C_INT64_T) :: cycles
    cycles = elapsed_cycles()
  end function
  function timer_cycles_raw() result(cycles)   ! read time base counter in ticks (no fencing)
    implicit none
    integer(C_INT64_T) :: cycles
    cycles = elapsed_cycles_raw()
  end function
  function timer_us() result(us)               ! read wall clock time in microseconds
    implicit none
    integer(C_INT64_T) :: us
    us = elapsed_us()
  end function
  function timer_freq() result(freq)           ! get time base counter frequency in Hertz
    implicit none
    integer(C_INT64_T) :: freq
    freq = cycles_counter_freq()
  end function
end module
