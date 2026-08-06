! SPDX-License-Identifier: GPL-2.0-or-later
program test_scalar_filter
  use fkf_module
  use test_support
  implicit none

  type(fkf_result) :: fit
  real(dp) :: a0(1), p0(1, 1), dt(1, 1), ct(1, 1)
  real(dp) :: tt(1, 1, 1), zt(1, 1, 1), hht(1, 1, 1), ggt(1, 1, 1), y(1, 2)

  a0 = 0.0_dp
  p0 = 1.0_dp
  dt = 0.0_dp
  ct = 0.0_dp
  tt = 1.0_dp
  zt = 1.0_dp
  hht = 0.1_dp
  ggt = 0.2_dp
  y(1, :) = [1.0_dp, 2.0_dp]

  call fkf(a0, p0, dt, ct, tt, zt, hht, ggt, y, fit, .true.)
  call assert_true(fit%status == fkf_success, 'scalar filter status')
  call assert_close(fit%att(1, 1), 5.0_dp / 6.0_dp, 1.0e-13_dp, 'att 1')
  call assert_close(fit%ptt(1, 1, 1), 1.0_dp / 6.0_dp, 1.0e-13_dp, 'ptt 1')
  call assert_close(fit%att(1, 2), 1.5_dp, 1.0e-13_dp, 'att 2')
  call assert_close(fit%ptt(1, 1, 2), 4.0_dp / 35.0_dp, 1.0e-13_dp, 'ptt 2')
  call assert_close(fit%at(1, 3), 1.5_dp, 1.0e-13_dp, 'prediction 3')
  call assert_close(fit%pt(1, 1, 3), 3.0_dp / 14.0_dp, 1.0e-13_dp, 'prediction variance 3')
  call assert_close(fit%log_likelihood, -3.422967818782874_dp, 1.0e-13_dp, 'log likelihood')
  call finish_test('scalar filter fixtures')
end program test_scalar_filter
