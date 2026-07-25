! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
program test_maps
  use fnonlinear, only : dp, tent_sim, logistic_sim, henon_sim, ikeda_sim, lorenz_sim, rossler_sim
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(dp), allocatable :: x(:), xy(:, :), trajectory(:, :)
  real(dp) :: times(2)
  integer :: status

  call tent_sim(5, 0, 1.9_dp, x, status, start=0.2_dp)
  call check(status == 0, "tent status")
  call check(maxval(abs(x - [0.2_dp, 0.38_dp, 0.722_dp, 0.5282_dp, 0.89642_dp])) < 1.0e-13_dp, &
    "tent recurrence")

  call logistic_sim(5, 0, 3.7_dp, x, status, start=0.2_dp)
  call check(status == 0, "logistic status")
  call check(maxval(abs(x - [0.2_dp, 0.592_dp, 0.8936832_dp, 0.35155009073971194_dp, &
    0.8434617104302653_dp])) < 1.0e-13_dp, "logistic recurrence")

  call henon_sim(5, 0, 1.4_dp, 0.3_dp, xy, status, start=[0.1_dp, 0.2_dp])
  call check(status == 0 .and. all(shape(xy) == [5, 2]), "Henon status and shape")
  call check(maxval(abs(xy(:, 1) - [0.1_dp, 1.046_dp, -0.5017624_dp, &
    0.9613282915247359_dp, -0.44434163772021446_dp])) < 1.0e-13_dp, "Henon x")
  call check(maxval(abs(xy(:, 2) - [0.2_dp, 0.1_dp, 1.046_dp, -0.5017624_dp, &
    0.9613282915247359_dp])) < 1.0e-13_dp, "Henon y")

  call ikeda_sim(5, 0, 0.4_dp, 6.0_dp, 0.9_dp, xy, status, start=[0.1_dp, 0.2_dp])
  call check(status == 0 .and. all(shape(xy) == [5, 2]), "Ikeda status and shape")
  call check(abs(xy(2, 1) - 0.9025912530985065_dp) < 1.0e-13_dp, "Ikeda real")
  call check(abs(xy(2, 2) - 0.17610092568490598_dp) < 1.0e-13_dp, "Ikeda imaginary")
  call check(all(ieee_is_finite(xy)), "Ikeda finite")

  times = [0.0_dp, 0.01_dp]
  call lorenz_sim(times, [16.0_dp, 45.92_dp, 4.0_dp], [-14.0_dp, -13.0_dp, 47.0_dp], &
    trajectory, status)
  call check(status == 0 .and. all(shape(trajectory) == [2, 4]), "Lorenz status and shape")
  call check(maxval(abs(trajectory(2, 2:4) - [-13.8311843488_dp, -12.7265657036_dp, &
    46.9116056197_dp])) < 1.0e-8_dp, "Lorenz RK4")

  call rossler_sim(times, [0.2_dp, 0.2_dp, 8.0_dp], [-1.894_dp, -9.92_dp, 0.025_dp], &
    trajectory, status)
  call check(status == 0, "Rossler status")
  call check(maxval(abs(trajectory(2, 2:4) - [-1.7948553734_dp, -9.9583230796_dp, &
    0.0245610102_dp])) < 5.0e-9_dp, "Rossler RK4")

  print '(a)', "Chaotic-map and RK4 tests passed."
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // message
      error stop 1
    end if
  end subroutine check
end program test_maps
