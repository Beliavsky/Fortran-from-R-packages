! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program test_data
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use strand
  implicit none
  type(cross_section_state) :: first, second
  type(cross_section_stats) :: stats
  type(portfolio_state) :: portfolio
  real(dp) :: nan
  integer, allocatable :: shares(:, :)

  nan = ieee_value(0.0_dp, ieee_quiet_nan)
  first = update_cross_section(new_id=[1, 2], new_values=reshape([10.0_dp, 20.0_dp, 100.0_dp, 200.0_dp], [2, 2]), &
    date=1, replace_columns=[1], replace_values=[0.0_dp])
  second = update_cross_section(first, [2, 3], reshape([21.0_dp, 30.0_dp, 210.0_dp, nan], [2, 2]), 2, &
    carry_columns=[2], carry_values=[0.0_dp], replace_columns=[2], replace_values=[-1.0_dp])
  call assert_true(size(second%id) == 3)
  call assert_true(second%id(3) == 1)
  call assert_true(second%carry_forward(3))
  call assert_close(second%values(3, 2), 0.0_dp, 0.0_dp)
  call assert_close(second%values(2, 2), -1.0_dp, 0.0_dp)
  stats = compute_cross_section_stats(second, first)
  call assert_true(stats%input_rows == 3)
  call assert_true(stats%carry_forward_rows == 1)

  call initialize_portfolio(portfolio, 2, 2, reshape([10, 20, -5, 0], [2, 2]))
  call apply_adjustment_ratio(portfolio, [0.5_dp, 1.0_dp])
  shares = consolidated_shares(portfolio)
  call assert_true(all(shares(:, 1) == [20, 20]))
  call assert_true(all(shares(:, 2) == [-10, 0]))

  print '(a)', 'test_data: PASS'
contains
  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) error stop 1
  end subroutine assert_close
  subroutine assert_true(value)
    logical, intent(in) :: value
    if (.not. value) error stop 1
  end subroutine assert_true
end program test_data
