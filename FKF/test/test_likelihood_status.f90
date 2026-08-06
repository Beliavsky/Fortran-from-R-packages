! SPDX-License-Identifier: GPL-2.0-or-later
program test_likelihood_status
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use fkf_module
  use test_support
  implicit none

  type(fkf_result) :: source_fit, corrected_fit, bad_fit
  real(dp) :: nanv, difference
  real(dp) :: a0(1), p0(1, 1), dt(1, 1), ct(1, 1)
  real(dp) :: tt(1, 1, 1), zt(1, 1, 1), hht(1, 1, 1), ggt(1, 1, 1), y(1, 3)

  nanv = ieee_value(0.0_dp, ieee_quiet_nan)
  a0 = 0.0_dp
  p0 = 1.0_dp
  dt = 0.0_dp
  ct = 0.0_dp
  tt = 1.0_dp
  zt = 1.0_dp
  hht = 0.1_dp
  ggt = 0.2_dp
  y(1, :) = [0.5_dp, nanv, -0.2_dp]

  call fkf(a0, p0, dt, ct, tt, zt, hht, ggt, y, source_fit)
  call fkf(a0, p0, dt, ct, tt, zt, hht, ggt, y, corrected_fit, .true.)
  difference = source_fit%log_likelihood - corrected_fit%log_likelihood
  call assert_close(difference, -0.5_dp * log(2.0_dp * acos(-1.0_dp)), 1.0e-13_dp, &
    'source missing-data likelihood constant')

  p0 = 0.0_dp
  hht = 0.0_dp
  ggt = 0.0_dp
  y(1, :) = [0.0_dp, 0.0_dp, 0.0_dp]
  call fkf(a0, p0, dt, ct, tt, zt, hht, ggt, y, bad_fit, .true.)
  call assert_true(bad_fit%status == fkf_non_pos_def, 'singular innovation covariance status')
  call assert_true(bad_fit%failure_time == 1, 'singular failure time')
  call finish_test('likelihood compatibility and failure status')
end program test_likelihood_status
