! SPDX-License-Identifier: GPL-3.0-only
program test_simplex
  use portfoliooptim, only : dp, lp_result, solve_lp
  implicit none
  real(dp) :: c(2), a(4, 2), b(4)
  real(dp) :: ai(2, 1), bi(2), ci(1)
  type(lp_result) :: result

  c = [-3.0_dp, -2.0_dp]
  a(1, :) = [1.0_dp, 1.0_dp]
  a(2, :) = [1.0_dp, 0.0_dp]
  a(3, :) = [0.0_dp, 1.0_dp]
  a(4, :) = [-1.0_dp, 0.0_dp]
  b = [4.0_dp, 2.0_dp, 3.0_dp, -1.0_dp]
  result = solve_lp(c, a, b)
  call assert_true(result%optimal)
  call assert_close(result%x(1), 2.0_dp, 1.0e-11_dp)
  call assert_close(result%x(2), 2.0_dp, 1.0e-11_dp)
  call assert_close(result%objective, -10.0_dp, 1.0e-11_dp)

  ci = 1.0_dp
  ai(:, 1) = [1.0_dp, -1.0_dp]
  bi = [0.0_dp, -1.0_dp]
  result = solve_lp(ci, ai, bi)
  call assert_true(.not. result%feasible)
  call assert_true(.not. result%optimal)
  print '(a)', 'test_simplex: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print *, 'mismatch:', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_simplex
