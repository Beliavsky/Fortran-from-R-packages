! Modern Fortran translation of cmprsk computational code.
! Upstream cmprsk Copyright (C) 2000 Robert Gray.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_gray_stratified
   use r_kinds, only : dp
   use cmprsk, only : gray_test_result, gray_test, cmprsk_success
   implicit none

   integer, parameter :: n = 18
   real(dp), parameter :: tol = 5.0e-13_dp
   real(dp) :: time(n)
   integer :: event_code(n)
   integer :: group(n)
   integer :: strata(n)
   integer :: status
   type(gray_test_result) :: gray

   time = [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 5.0_dp, 6.0_dp, &
           7.0_dp, 8.0_dp, 9.0_dp, 9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp, 12.0_dp, 13.0_dp]
   event_code = [1, 0, 2, 1, 2, 0, 1, 2, 0, 1, 2, 1, 0, 2, 1, 0, 2, 1]
   group = [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]
   strata = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2]

   call gray_test(time, event_code, group, strata, 0.7_dp, gray, status)
   if (status /= cmprsk_success) error stop 'gray_test status'
   if (maxval(abs(gray%score - [2.26522893145907744_dp, -2.06251094138504332_dp])) > tol) error stop 'gray score'
   if (abs(gray%covariance(1, 1) - 1.05801884838041493_dp) > tol) error stop 'gray variance 11'
   if (abs(gray%covariance(1, 2) + 0.566366058199916589_dp) > tol) error stop 'gray variance 12'
   if (abs(gray%covariance(2, 2) - 1.22532612550999609_dp) > tol) error stop 'gray variance 22'
   print *, 'test_gray_stratified: PASS'
end program test_gray_stratified
