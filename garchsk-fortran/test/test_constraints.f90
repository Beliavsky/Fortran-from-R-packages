! SPDX-License-Identifier: GPL-2.0-or-later
program test_constraints
   use garchsk, only : dp, garchsk_ineqfun, gjrsk_ineqfun, garchsk_parameters_valid, &
      gjrsk_parameters_valid, garchsk_lik
   use test_support, only : assert_close, assert_true
   implicit none
   real(dp), parameter :: data(4) = [0.01_dp, -0.02_dp, 0.015_dp, 0.005_dp]
   real(dp) :: pg(10), pj(13)
   real(dp), allocatable :: values(:)

   pg = [0.1_dp, 1.0e-4_dp, 0.08_dp, 0.85_dp, 0.0_dp, 0.03_dp, 0.60_dp, &
      0.60_dp, 0.05_dp, 0.75_dp]
   pj = [0.1_dp, 1.0e-4_dp, 0.06_dp, 0.04_dp, 0.84_dp, 0.0_dp, 0.02_dp, &
      -0.01_dp, 0.60_dp, 0.60_dp, 0.04_dp, 0.02_dp, 0.74_dp]
   call assert_true(garchsk_parameters_valid(pg), 'valid GARCHSK parameters rejected')
   call assert_true(gjrsk_parameters_valid(pj), 'valid GJRSK parameters rejected')
   values = garchsk_ineqfun(pg)
   call assert_true(size(values) == 12, 'wrong GARCHSK constraint count')
   call assert_close(values(10), 0.93_dp)
   values = gjrsk_ineqfun(pj)
   call assert_true(size(values) == 15, 'wrong GJRSK constraint count')
   call assert_close(values(13), 0.94_dp)
   pg(4) = 0.95_dp
   call assert_true(.not. garchsk_parameters_valid(pg), 'invalid persistence accepted')
   call assert_true(garchsk_lik(pg, data) >= 1.0e12_dp, 'invalid model did not receive penalty')
   print '(a)', 'test_constraints: PASS'
end program test_constraints
