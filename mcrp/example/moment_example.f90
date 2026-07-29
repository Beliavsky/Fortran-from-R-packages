! SPDX-License-Identifier: GPL-3.0-only
program moment_example
   use mcrp_module
   implicit none

   real(dp) :: r(6, 3), w(3)

   r = transpose(reshape([ &
      0.012_dp, -0.004_dp,  0.009_dp, &
     -0.006_dp,  0.011_dp,  0.003_dp, &
      0.018_dp,  0.002_dp, -0.005_dp, &
     -0.009_dp, -0.007_dp,  0.014_dp, &
      0.004_dp,  0.015_dp, -0.002_dp, &
      0.021_dp, -0.003_dp,  0.006_dp], [3, 6]))
   w = [0.40_dp, 0.35_dp, 0.25_dp]

   write(*, '(a, es14.6)') 'Portfolio variance: ', port_risk(r, w)
   write(*, '(a, es14.6)') 'Portfolio skewness: ', port_skew(r, w)
   write(*, '(a, es14.6)') 'Portfolio kurtosis: ', port_kurt(r, w)
   write(*, '(a, *(f10.6, 1x))') 'Variance contributions: ', &
      port_risk_contrib(r, w)
   write(*, '(a, *(f10.6, 1x))') 'Skewness contributions: ', &
      port_skew_contrib(r, w)
   write(*, '(a, *(f10.6, 1x))') 'Kurtosis contributions: ', &
      port_kurt_contrib(r, w)
end program moment_example
