! SPDX-License-Identifier: GPL-2.0-or-later
program test_time_varying
  use fkf_module
  use test_support
  implicit none

  type(fkf_result) :: fit_tv, fit_manual
  real(dp) :: a0(1), p0(1, 1), dt(1, 3), ct(1, 3)
  real(dp) :: tt(1, 1, 3), zt(1, 1, 3), hht(1, 1, 3), ggt(1, 1, 3), y(1, 3)
  real(dp) :: a, p, v, f, k, att, ptt, ll
  integer :: i

  a0 = 0.2_dp
  p0 = 0.7_dp
  dt(1, :) = [0.1_dp, -0.2_dp, 0.0_dp]
  ct(1, :) = [0.0_dp, 0.3_dp, -0.1_dp]
  tt(1, 1, :) = [0.8_dp, 1.1_dp, 0.9_dp]
  zt(1, 1, :) = [1.0_dp, 0.5_dp, 1.5_dp]
  hht(1, 1, :) = [0.05_dp, 0.08_dp, 0.02_dp]
  ggt(1, 1, :) = [0.2_dp, 0.4_dp, 0.3_dp]
  y(1, :) = [0.4_dp, -0.1_dp, 0.8_dp]

  call fkf(a0, p0, dt, ct, tt, zt, hht, ggt, y, fit_tv, .true.)
  call assert_true(fit_tv%status == fkf_success, 'time-varying status')

  a = a0(1)
  p = p0(1, 1)
  ll = 0.0_dp
  do i = 1, 3
    v = y(1, i) - ct(1, i) - zt(1, 1, i) * a
    f = zt(1, 1, i)**2 * p + ggt(1, 1, i)
    k = p * zt(1, 1, i) / f
    att = a + k * v
    ptt = p - p * zt(1, 1, i) * k
    ll = ll - 0.5_dp * (log(2.0_dp * acos(-1.0_dp)) + log(f) + v * v / f)
    call assert_close(fit_tv%att(1, i), att, 1.0e-13_dp, 'time-varying state')
    call assert_close(fit_tv%ptt(1, 1, i), ptt, 1.0e-13_dp, 'time-varying variance')
    a = dt(1, i) + tt(1, 1, i) * att
    p = tt(1, 1, i)**2 * ptt + hht(1, 1, i)
  end do
  call assert_close(fit_tv%log_likelihood, ll, 1.0e-13_dp, 'time-varying likelihood')

  call fkf(a0, p0, dt(:, 1:1), ct(:, 1:1), tt(:, :, 1:1), zt(:, :, 1:1), &
    hht(:, :, 1:1), ggt(:, :, 1:1), y, fit_manual, .true.)
  call assert_true(fit_manual%status == fkf_success, 'constant-parameter extent')
  call finish_test('time-varying and constant parameter extents')
end program test_time_varying
