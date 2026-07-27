! SPDX-License-Identifier: GPL-2.0-or-later
program basic_garchsk
   use garchsk, only : dp, moment_path, garchsk_construct, skewness, kurtosis
   implicit none
   real(dp), parameter :: x(6) = [0.01_dp, -0.015_dp, 0.008_dp, 0.012_dp, -0.006_dp, 0.004_dp]
   real(dp), parameter :: p(10) = [0.05_dp, 1.0e-4_dp, 0.06_dp, 0.88_dp, 0.0_dp, &
      0.02_dp, 0.60_dp, 0.60_dp, 0.04_dp, 0.75_dp]
   type(moment_path) :: path
   path = garchsk_construct(p, x)
   print '(a,f10.6)', 'sample skewness: ', skewness(x)
   print '(a,f10.6)', 'sample kurtosis: ', kurtosis(x)
   print '(a,es12.4)', 'last conditional variance: ', path%h(size(x))
end program basic_garchsk
