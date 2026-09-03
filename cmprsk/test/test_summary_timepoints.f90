! Modern Fortran translation of cmprsk computational code.
! Upstream cmprsk Copyright (C) 2000 Robert Gray.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_summary_timepoints
   use r_kinds, only : dp
   use cmprsk, only : crr_result, crr_summary_result, cuminc_result, fit_crr, fit_cuminc, &
                      summarize_crr, cuminc_timepoints, cmprsk_success
   implicit none

   integer, parameter :: n = 12
   real(dp), allocatable :: estimate(:, :)
   logical, allocatable :: present_value(:, :)
   real(dp) :: requested(2)
   real(dp) :: time(n)
   real(dp), allocatable :: variance(:, :)
   real(dp) :: x(n, 1)
   integer :: cause(n)
   integer :: group(n)
   integer :: i
   integer :: status
   type(crr_result) :: regression
   type(crr_summary_result) :: summary
   type(cuminc_result) :: incidence

   time = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 8.0_dp]
   cause = [1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2]
   group = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
   do i = 1, n
      x(i, 1) = real(i - 1, dp)/real(n - 1, dp)
   end do

   call fit_crr(time, cause, x, regression, status, censor_group=group)
   if (status /= cmprsk_success) error stop 'fit_crr failed'
   call summarize_crr(regression, summary, status)
   if (status /= cmprsk_success) error stop 'summarize_crr failed'
   if (summary%standard_error(1) <= 0.0_dp) error stop 'summary standard error'
   if (abs(summary%relative_risk(1) - exp(regression%coefficients(1))) > 1.0e-13_dp) error stop 'summary RR'
   if (abs(summary%likelihood_ratio + 2.0_dp*(regression%loglik_null - regression%loglik)) > 1.0e-13_dp) &
      error stop 'summary likelihood ratio'

   call fit_cuminc(time, cause, incidence, status, group=group)
   if (status /= cmprsk_success) error stop 'fit_cuminc failed'
   requested = [3.0_dp, 8.0_dp]
   call cuminc_timepoints(incidence, requested, estimate, variance, present_value)
   if (size(estimate, 1) /= 4 .or. size(estimate, 2) /= 2) error stop 'timepoints dimensions'
   if (.not. all(present_value)) error stop 'timepoints support'
   print *, 'test_summary_timepoints: PASS'
end program test_summary_timepoints
