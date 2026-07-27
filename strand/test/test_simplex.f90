! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program test_simplex
  use strand_kinds, only : dp
  use strand_simplex, only : solve_bounded_lp
  use strand_types, only : lp_result
  implicit none
  type(lp_result) :: result
  real(dp) :: objective(2), a(3, 2), rhs(3), lower(2), upper(2)
  integer :: sense(3)

  objective = [-3.0_dp, -2.0_dp]
  a = reshape([1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [3, 2])
  rhs = [4.0_dp, 2.0_dp, 3.0_dp]
  sense = [-1, -1, -1]
  lower = 0.0_dp
  upper = huge(1.0_dp)
  result = solve_bounded_lp(objective, a, sense, rhs, lower, upper)
  call assert_true(result%optimal)
  call assert_close(result%x(1), 2.0_dp, 1.0e-10_dp)
  call assert_close(result%x(2), 2.0_dp, 1.0e-10_dp)
  call assert_close(result%objective, -10.0_dp, 1.0e-10_dp)

  objective = 0.0_dp
  a = 0.0_dp
  a(1, 1) = 1.0_dp
  rhs = 0.0_dp
  rhs(1) = 2.0_dp
  sense = -1
  sense(1) = 1
  upper = [1.0_dp, 1.0_dp]
  result = solve_bounded_lp(objective, a(1:1, :), sense(1:1), rhs(1:1), lower, upper)
  call assert_true(.not. result%optimal)
  call assert_true(.not. result%feasible)

  print '(a)', 'test_simplex: PASS'
contains
  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print *, 'mismatch:', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(value)
    logical, intent(in) :: value
    if (.not. value) error stop 1
  end subroutine assert_true
end program test_simplex
