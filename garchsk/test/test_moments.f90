! SPDX-License-Identifier: GPL-2.0-or-later
program test_moments
   use garchsk, only : dp, skewness, kurtosis, moment_path, garchsk_construct, gjrsk_construct
   use test_support, only : assert_close, assert_true
   implicit none
   real(dp), parameter :: data(8) = [0.012_dp, -0.008_dp, 0.015_dp, -0.020_dp, &
      0.006_dp, 0.011_dp, -0.004_dp, 0.009_dp]
   real(dp), parameter :: pg(10) = [0.1_dp, 1.0e-4_dp, 0.08_dp, 0.85_dp, 0.0_dp, &
      0.03_dp, 0.60_dp, 0.60_dp, 0.05_dp, 0.75_dp]
   real(dp), parameter :: pj(13) = [0.1_dp, 1.0e-4_dp, 0.06_dp, 0.04_dp, 0.84_dp, &
      0.0_dp, 0.02_dp, -0.01_dp, 0.60_dp, 0.60_dp, 0.04_dp, 0.02_dp, 0.74_dp]
   type(moment_path) :: path

   call assert_close(skewness(data), -0.661244062378995_dp)
   call assert_close(kurtosis(data), 1.809605363019897_dp)

   path = garchsk_construct(pg, data)
   call assert_true(path%success, path%message)
   call assert_close(path%mu(2), 0.0012_dp)
   call assert_close(path%h(2), 0.000232330357142857_dp)
   call assert_close(path%skewness(8), -0.02215805346063205_dp)
   call assert_close(path%kurtosis(8), 2.370485600352341_dp)

   path = gjrsk_construct(pj, data)
   call assert_true(path%success, path%message)
   call assert_close(path%h(2), 0.0002290984375_dp)
   call assert_close(path%skewness(8), -0.017815637173621182_dp)
   call assert_close(path%kurtosis(8), 2.3000442319521515_dp)
   print '(a)', 'test_moments: PASS'
end program test_moments
