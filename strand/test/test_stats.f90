! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program test_stats
  use strand
  implicit none
  real(dp), allocatable :: z(:), adj(:)
  real(dp) :: x(3), y(8), factors(8, 2), dd
  integer :: groups(8, 1)

  x = [10.0_dp, 1.0_dp, 100.0_dp]
  z = rank_normal(x)
  call assert_close(z(1), 0.0_dp, 1.0e-12_dp)
  call assert_close(z(2), -1.0_dp, 1.0e-12_dp)
  call assert_close(z(3), 1.0_dp, 1.0e-12_dp)

  y = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp]
  groups(:, 1) = [1, 1, 1, 1, 2, 2, 2, 2]
  z = normalize_grouped(y, groups, loops=2)
  call assert_close(sum(z(1:4)) / 4.0_dp, 0.0_dp, 1.0e-12_dp)
  call assert_close(sum(z(5:8)) / 4.0_dp, 0.0_dp, 1.0e-12_dp)

  factors(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
  factors(:, 2) = [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 4.0_dp]
  adj = adjust_numeric(y, factors, loops=2)
  call assert_true(abs(sum(adj)) < 1.0e-10_dp)
  call assert_vector_close(adj, [-0.5336900676524256_dp, 0.5336900676524259_dp, &
    0.1731072048738221_dp, 0.9475089173720049_dp, -0.9475089173720044_dp, &
    -1.5124270755347036_dp, 1.5124270755347022_dp, -0.1731072048738215_dp], 1.0e-12_dp)

  dd = maximum_drawdown([0.10_dp, -0.05_dp, -0.10_dp, 0.20_dp])
  call assert_close(dd, -0.15_dp, 1.0e-14_dp)
  call assert_close(factor_exposure([100.0_dp, -50.0_dp], [1.0_dp, 2.0_dp], 1000.0_dp), &
    0.0_dp, 1.0e-14_dp)
  call assert_close(category_exposure([100.0_dp, -50.0_dp, 25.0_dp], [1, 2, 1], 1, 1000.0_dp), &
    0.125_dp, 1.0e-14_dp)

  print '(a)', 'test_stats: PASS'
contains
  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print *, 'mismatch:', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_vector_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    if (size(actual) /= size(expected)) error stop 1
    if (maxval(abs(actual - expected)) > tolerance) then
      print *, 'vector mismatch:', maxval(abs(actual - expected))
      error stop 1
    end if
  end subroutine assert_vector_close
  subroutine assert_true(value)
    logical, intent(in) :: value
    if (.not. value) error stop 1
  end subroutine assert_true
end program test_stats
