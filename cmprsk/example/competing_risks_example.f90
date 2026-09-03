! Modern Fortran translation of cmprsk computational code.
! Upstream cmprsk Copyright (C) 2000 Robert Gray.
! SPDX-License-Identifier: GPL-2.0-or-later
program competing_risks_example
   use r_kinds, only : dp
   use cmprsk, only : crr_result, cuminc_result, fit_crr, fit_cuminc, cmprsk_success
   implicit none

   integer, parameter :: n = 12
   real(dp) :: time(n)
   real(dp) :: x(n, 1)
   integer :: cause(n)
   integer :: group(n)
   integer :: i
   integer :: status
   type(crr_result) :: regression
   type(cuminc_result) :: incidence

   time = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 8.0_dp]
   cause = [1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2]
   group = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
   do i = 1, n
      x(i, 1) = real(i - 1, dp)/real(n - 1, dp)
   end do

   call fit_cuminc(time, cause, incidence, status, group=group)
   if (status /= cmprsk_success) error stop 'fit_cuminc failed'
   call fit_crr(time, cause, x, regression, status, censor_group=group)
   if (status /= cmprsk_success) error stop 'fit_crr failed'

   print '(a,2f10.5)', 'Gray statistic, p-value: ', incidence%tests(1)%statistic, incidence%tests(1)%p_value
   print '(a,f10.5)', 'Fine-Gray coefficient: ', regression%coefficients(1)
end program competing_risks_example
