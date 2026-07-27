! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a translation of tseriesChaos and is distributed
! under the GNU General Public License version 2 only.
program test_core
  use tserieschaos, only : dp, delay_embed, delay_embed_lags, delay_embed_matrix, &
    correlation_integral, correlation_dimension_curve, average_mutual_information, &
    recurrence_distance_matrix, space_time_separation
  implicit none
  real(dp), allocatable :: embedded(:, :), eps(:), curve(:, :), ami(:), recurrence(:, :), isolines(:, :)
  real(dp) :: x(6), c2
  real(dp) :: xm(5, 2)
  integer :: status, i

  x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
  call delay_embed(x, 3, 2, embedded, status)
  call check(status == 0, "delay_embed status")
  call check(all(shape(embedded) == [2, 3]), "delay_embed shape")
  call check(maxval(abs(embedded - reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [2, 3]))) < 1.0e-14_dp, &
    "delay_embed values")

  call delay_embed_lags(x, [0, 1, 3], embedded, status)
  call check(status == 0 .and. all(shape(embedded) == [3, 3]), "delay_embed_lags shape")
  call check(maxval(abs(embedded(1, :) - [1.0_dp, 2.0_dp, 4.0_dp])) < 1.0e-14_dp, "delay_embed_lags values")

  do i = 1, 5
    xm(i, 1) = real(i, dp)
    xm(i, 2) = 10.0_dp + real(i, dp)
  end do
  call delay_embed_matrix(xm, [0, 2], embedded, status)
  call check(status == 0 .and. all(shape(embedded) == [3, 4]), "delay_embed_matrix shape")
  call check(maxval(abs(embedded(1, :) - [1.0_dp, 11.0_dp, 3.0_dp, 13.0_dp])) < 1.0e-14_dp, &
    "delay_embed_matrix values")

  call correlation_integral([0.0_dp, 1.0_dp, 2.0_dp], 1, 1, 1, 1.1_dp, c2, status)
  call check(status == 0, "correlation_integral status")
  call check(abs(c2 - 2.0_dp / 3.0_dp) < 1.0e-14_dp, "correlation_integral exact")

  call correlation_dimension_curve([0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], 2, 1, 1, 0.1_dp, 8, eps, curve, status)
  call check(status == 0, "correlation_dimension_curve status")
  call check(all(shape(curve) == [8, 2]), "correlation_dimension_curve shape")
  call check(all(eps(2:) > eps(:size(eps)-1)), "epsilon grid increasing")
  call check(all(curve(2:, :) >= curve(:size(curve, 1)-1, :)), "correlation curves monotone")
  call check(abs(curve(8, 1) - 1.0_dp) < 1.0e-12_dp, "dimension-one terminal curve")
  call check(abs(curve(8, 2) - 1.0_dp / 3.0_dp) < 1.0e-12_dp, "dimension-two terminal curve")

  call average_mutual_information([0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], 2, 2, ami, status)
  call check(status == 0, "average_mutual_information status")
  call check(lbound(ami, 1) == 0 .and. ubound(ami, 1) == 2, "AMI bounds")
  call check(abs(ami(0) - log(2.0_dp)) < 1.0e-14_dp, "AMI lag zero")
  call check(abs(ami(1) + 0.6_dp * log(0.6_dp) + 0.4_dp * log(0.4_dp)) < 1.0e-14_dp, "AMI lag one")

  call recurrence_distance_matrix([0.0_dp, 1.0_dp, 2.0_dp], 1, 1, recurrence, status)
  call check(status == 0, "recurrence status")
  call check(maxval(abs(recurrence - reshape([0.0_dp, 0.5_dp, 1.0_dp, &
                                               0.5_dp, 0.0_dp, 0.5_dp, &
                                               1.0_dp, 0.5_dp, 0.0_dp], [3, 3]))) < 1.0e-14_dp, &
    "recurrence values")

  call space_time_separation([0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp], &
    2, 1, 1, 4, isolines, status)
  call check(status == 0 .and. all(shape(isolines) == [10, 4]), "space_time_separation shape")
  call check(all(isolines(2:, :) >= isolines(:9, :)), "space-time quantiles monotone")
  call check(all(isolines >= 0.0_dp), "space-time distances nonnegative")

  print '(a)', "Core embedding and nonlinear-statistic tests passed."
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // message
      error stop 1
    end if
  end subroutine check
end program test_core
