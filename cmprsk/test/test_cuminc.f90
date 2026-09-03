! Modern Fortran translation of cmprsk computational code.
! Upstream cmprsk Copyright (C) 2000 Robert Gray.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_cuminc
   use r_kinds, only : dp
   use cmprsk, only : cuminc_curve, gray_test_result, cumulative_incidence, gray_test, cmprsk_success
   implicit none

   integer, parameter :: n = 12
   real(dp), parameter :: tol = 5.0e-13_dp
   real(dp) :: time(n)
   integer :: event_code(n)
   integer :: failed_any(n)
   integer :: failed_cause(n)
   integer :: group(n)
   integer :: strata(n)
   integer :: status
   type(cuminc_curve) :: curve
   type(gray_test_result) :: gray

   time = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 8.0_dp]
   event_code = [1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2]
   group = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
   strata = 1
   failed_any = merge(1, 0, event_code > 0)
   failed_cause = merge(1, 0, event_code == 1)

   call cumulative_incidence(time, failed_any, failed_cause, curve, status)
   if (status /= cmprsk_success) error stop 'cumulative_incidence status'
   if (size(curve%time) /= 10) error stop 'cumulative_incidence size'
   if (abs(curve%estimate(10) - 0.41071428571428564_dp) > tol) error stop 'cumulative_incidence estimate'
   if (abs(curve%variance(10) - 0.034832977493981176_dp) > tol) error stop 'cumulative_incidence variance'

   call gray_test(time, event_code, group, strata, 0.2_dp, gray, status)
   if (status /= cmprsk_success) error stop 'gray_test status'
   if (abs(gray%score(1) - 0.029891559491206698_dp) > tol) error stop 'gray_test score'
   if (abs(gray%covariance(1, 1) - 1.0147922836339014_dp) > tol) error stop 'gray_test covariance'
   print *, 'test_cuminc: PASS'
end program test_cuminc
