! SPDX-License-Identifier: GPL-2.0-or-later
program test_likelihood_forecast
   use garchsk, only : dp, forecast_result, garchsk_lik, gjrsk_lik, garchsk_fcst, gjrsk_fcst
   use test_support, only : assert_close, assert_true
   implicit none
   real(dp), parameter :: data(8) = [0.012_dp, -0.008_dp, 0.015_dp, -0.020_dp, &
      0.006_dp, 0.011_dp, -0.004_dp, 0.009_dp]
   real(dp), parameter :: pg(10) = [0.1_dp, 1.0e-4_dp, 0.08_dp, 0.85_dp, 0.0_dp, &
      0.03_dp, 0.60_dp, 0.60_dp, 0.05_dp, 0.75_dp]
   real(dp), parameter :: pj(13) = [0.1_dp, 1.0e-4_dp, 0.06_dp, 0.04_dp, 0.84_dp, &
      0.0_dp, 0.02_dp, -0.01_dp, 0.60_dp, 0.60_dp, 0.04_dp, 0.02_dp, 0.74_dp]
   type(forecast_result) :: fcst

   call assert_close(garchsk_lik(pg, data), -25.27094803375209_dp, 1.0e-11_dp)
   call assert_close(gjrsk_lik(pj, data), -25.254488186577028_dp, 1.0e-11_dp)

   fcst = garchsk_fcst(pg, data, 4)
   call assert_true(fcst%success, fcst%message)
   call assert_close(fcst%mu(1), 0.0009_dp)
   call assert_close(fcst%mu(2), 0.00009_dp)
   call assert_close(fcst%h(2), pg(2) + sum(pg(3:4))*fcst%h(1))

   fcst = gjrsk_fcst(pj, data, 4)
   call assert_true(fcst%success, fcst%message)
   call assert_close(fcst%mu(1), 0.0009_dp)
   call assert_close(fcst%mu(2), 0.00009_dp)
   call assert_close(fcst%h(2), pj(2) + sum(pj(3:5))*fcst%h(1))
   print '(a)', 'test_likelihood_forecast: PASS'
end program test_likelihood_forecast
