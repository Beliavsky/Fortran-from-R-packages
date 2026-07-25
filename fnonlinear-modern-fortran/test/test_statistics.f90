! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2005 Antonio Fabio Di Narzo
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
program test_statistics
  use fnonlinear, only : dp, delay_embed, delay_embed_lags, delay_embed_matrix, &
    correlation_integral, correlation_dimension_curve, mutual_information_curve, &
    recurrence_distance_matrix, recurrence_matrix, space_time_separation
  implicit none
  real(dp), allocatable :: embedded(:, :), eps(:), curve(:, :), ami(:), distances(:, :)
  real(dp), allocatable :: isolines(:, :)
  logical, allocatable :: recurrence(:, :)
  real(dp) :: x(6), c2
  real(dp) :: xm(5, 2)
  integer :: status, i

  x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
  call delay_embed(x, 3, 2, embedded, status)
  call check(status == 0 .and. all(shape(embedded) == [2, 3]), "delay embedding")
  call check(maxval(abs(embedded - reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
    5.0_dp, 6.0_dp], [2, 3]))) < 1.0e-14_dp, "delay embedding values")

  call delay_embed_lags(x, [0, 1, 3], embedded, status)
  call check(status == 0 .and. all(shape(embedded) == [3, 3]), "explicit lags")
  call check(maxval(abs(embedded(1, :) - [1.0_dp, 2.0_dp, 4.0_dp])) < 1.0e-14_dp, &
    "explicit lag values")

  do i = 1, 5
    xm(i, 1) = real(i, dp)
    xm(i, 2) = 10.0_dp + real(i, dp)
  end do
  call delay_embed_matrix(xm, [0, 2], embedded, status)
  call check(status == 0 .and. all(shape(embedded) == [3, 4]), "matrix embedding")

  call correlation_integral([0.0_dp, 1.0_dp, 2.0_dp], 1, 1, 1, 1.1_dp, c2, status)
  call check(status == 0 .and. abs(c2 - 2.0_dp / 3.0_dp) < 1.0e-14_dp, &
    "correlation integral")

  call correlation_dimension_curve([0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp], &
    2, 1, 1, 0.1_dp, 8, eps, curve, status)
  call check(status == 0 .and. all(shape(curve) == [8, 2]), "dimension curve")
  call check(all(curve(2:, :) >= curve(:size(curve, 1) - 1, :)), "dimension monotonicity")

  call mutual_information_curve([0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], &
    2, 2, ami, status)
  call check(status == 0 .and. lbound(ami, 1) == 0 .and. ubound(ami, 1) == 2, &
    "mutual information bounds")
  call check(abs(ami(0) - log(2.0_dp)) < 1.0e-14_dp, "mutual information lag zero")

  call recurrence_distance_matrix([0.0_dp, 1.0_dp, 2.0_dp], 1, 1, distances, status)
  call check(status == 0 .and. all(shape(distances) == [3, 3]), "recurrence distances")
  call check(maxval(abs(distances - transpose(distances))) < 1.0e-14_dp, &
    "recurrence distance symmetry")

  call recurrence_matrix([0.0_dp, 1.0_dp, 2.0_dp], 1, 1, 3, 1.1_dp, recurrence, status)
  call check(status == 0 .and. all(shape(recurrence) == [3, 3]), "recurrence matrix")
  call check(all(recurrence .eqv. transpose(recurrence)), "recurrence symmetry")
  call check(all([(recurrence(i, i), i = 1, 3)]), "recurrence diagonal")
  call check(recurrence(1, 2) .and. .not. recurrence(1, 3), "recurrence threshold")

  call space_time_separation([0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, &
    1.0_dp, 0.0_dp, -1.0_dp], 2, 1, 1, 4, isolines, status)
  call check(status == 0 .and. all(shape(isolines) == [10, 4]), "space-time separation")
  call check(all(isolines(2:, :) >= isolines(:9, :)), "space-time quantiles")

  print '(a)', "Embedding, entropy, recurrence, and dimension tests passed."
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAILED: " // message
      error stop 1
    end if
  end subroutine check
end program test_statistics
