! SPDX-License-Identifier: GPL-2.0-or-later
program garchsk_demo
   use garchsk, only : dp, moment_path, forecast_result, garchsk_construct, garchsk_fcst, garchsk_lik
   implicit none
   real(dp), parameter :: data(8) = [0.012_dp, -0.008_dp, 0.015_dp, -0.020_dp, &
      0.006_dp, 0.011_dp, -0.004_dp, 0.009_dp]
   real(dp), parameter :: params(10) = [0.1_dp, 1.0e-4_dp, 0.08_dp, 0.85_dp, 0.0_dp, &
      0.03_dp, 0.60_dp, 0.60_dp, 0.05_dp, 0.75_dp]
   type(moment_path) :: path
   type(forecast_result) :: fcst

   path = garchsk_construct(params, data)
   fcst = garchsk_fcst(params, data, 5)
   print '(a,es14.6)', 'negative log likelihood: ', garchsk_lik(params, data)
   print '(a,es14.6)', 'last variance:           ', path%h(size(data))
   print '(a,5(1x,es12.4))', 'variance forecasts:', fcst%h
end program garchsk_demo
