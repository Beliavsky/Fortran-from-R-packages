! SPDX-License-Identifier: GPL-3.0-or-later
program test_errors
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use intrinsicfrp, only: dp, vector_result, screening_result, rank_test_result
  use intrinsicfrp, only: oracle_result, oracle_control, tfrp, gkr_factor_screening
  use intrinsicfrp, only: iterative_kleibergen_paap_2006_beta_rank_test, oracle_tfrp
  implicit none
  real(dp) :: r(5, 3), f(4, 2), r2(8, 3), f2(8, 3), penalties(1)
  type(vector_result) :: v
  type(screening_result) :: s
  type(rank_test_result) :: rank
  type(oracle_result) :: o
  type(oracle_control) :: ctl
  integer :: i

  r = 0.0_dp
  f = 0.0_dp
  call tfrp(r, f, v)
  call check(v%status /= 0, 'dimension error')

  r2 = 0.0_dp
  f2 = 0.0_dp
  r2(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
  call tfrp(r2, f2, v)
  call check(v%status /= 0, 'nonfinite error')

  r2 = reshape([(0.01_dp * real(i, dp), i = 1, 24)], shape(r2))
  f2 = reshape([(0.02_dp * real(i, dp), i = 1, 24)], shape(f2))
  call gkr_factor_screening(r2, f2, s, target_level=1.2_dp)
  call check(s%status /= 0, 'screening level error')

  call iterative_kleibergen_paap_2006_beta_rank_test(r2, f2, rank)
  call check(rank%status /= 0, 'rank dimension error')

  penalties = 0.0_dp
  ctl = oracle_control()
  ctl%tuning_type = 'x'
  call oracle_tfrp(r2, f2(:, 1:2), penalties, o, ctl)
  call check(o%status /= 0, 'oracle tuning error')

  print '(a)', 'test_errors: PASS'

contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,1x,a)') 'FAIL:', trim(label)
      error stop 1
    end if
  end subroutine check
end program test_errors
