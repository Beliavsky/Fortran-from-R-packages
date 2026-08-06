! SPDX-License-Identifier: GPL-2.0-or-later
program test_missing_filter
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
  use fkf_module
  use test_support
  implicit none

  type(fkf_result) :: fit
  real(dp) :: nanv
  real(dp) :: a0(2), p0(2, 2), dt(2, 1), ct(2, 1)
  real(dp) :: tt(2, 2, 1), zt(2, 2, 1), hht(2, 2, 1), ggt(2, 2, 1), y(2, 4)
  real(dp) :: v, f, k1, k2

  nanv = ieee_value(0.0_dp, ieee_quiet_nan)
  a0 = [0.0_dp, 1.0_dp]
  p0 = reshape([2.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
  dt = 0.0_dp
  ct = 0.0_dp
  tt(:, :, 1) = reshape([1.0_dp, 0.0_dp, 0.2_dp, 1.0_dp], [2, 2])
  zt(:, :, 1) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
  hht(:, :, 1) = reshape([0.1_dp, 0.0_dp, 0.0_dp, 0.2_dp], [2, 2])
  ggt(:, :, 1) = reshape([0.3_dp, 0.0_dp, 0.0_dp, 0.4_dp], [2, 2])
  y(:, 1) = [1.0_dp, nanv]
  y(:, 2) = [nanv, 0.5_dp]
  y(:, 3) = [nanv, nanv]
  y(:, 4) = [1.5_dp, 0.2_dp]

  call fkf(a0, p0, dt, ct, tt, zt, hht, ggt, y, fit, .true.)
  call assert_true(fit%status == fkf_success, 'missing-data filter status')

  v = 1.0_dp
  f = 2.3_dp
  k1 = 2.0_dp / f
  k2 = 0.5_dp / f
  call assert_close(fit%att(1, 1), k1 * v, 1.0e-13_dp, 'partial update state 1')
  call assert_close(fit%att(2, 1), 1.0_dp + k2 * v, 1.0e-13_dp, 'partial update state 2')
  call assert_true(ieee_is_nan(fit%vt(2, 1)), 'missing innovation remains NaN')
  call assert_true(ieee_is_nan(fit%kt(1, 2, 1)), 'missing gain column remains NaN')
  call assert_close(fit%att(1, 3), fit%at(1, 3), 1.0e-13_dp, 'all-missing state unchanged')
  call assert_close(fit%ptt(2, 2, 3), fit%pt(2, 2, 3), 1.0e-13_dp, 'all-missing covariance unchanged')
  call assert_true(all(ieee_is_nan(fit%vt(:, 3))), 'all-missing innovations are NaN')
  call finish_test('partial and complete missing observations')
end program test_missing_filter
