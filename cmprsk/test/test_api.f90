! Modern Fortran translation of cmprsk computational code.
! Upstream cmprsk Copyright (C) 2000 Robert Gray.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_api
   use r_kinds, only : dp
   use cmprsk, only : cuminc_result, fit_cuminc, curve_timepoints, cmprsk_success
   implicit none

   integer, parameter :: n = 12
   real(dp) :: estimate(3)
   real(dp) :: requested(3)
   real(dp) :: time(n)
   real(dp) :: variance(3)
   integer :: event_code(n)
   integer :: group(n)
   integer :: status
   logical :: present_value(3)
   type(cuminc_result) :: fit

   time = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 8.0_dp]
   event_code = [1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2]
   group = [10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20]
   call fit_cuminc(time, event_code, fit, status, group=group, rho=0.2_dp)
   if (status /= cmprsk_success) error stop 'fit_cuminc status'
   if (size(fit%estimates) /= 4 .or. size(fit%tests) /= 2) error stop 'fit_cuminc sizes'
   if (fit%estimates(1)%group_label /= 10 .or. fit%estimates(1)%cause_label /= 1) error stop 'fit_cuminc labels'
   if (abs(fit%tests(1)%score(1) - 0.029891559491206698_dp) > 1.0e-12_dp) error stop 'fit_cuminc Gray score'

   requested = [0.5_dp, 5.0_dp, 8.0_dp]
   call curve_timepoints(fit%estimates(1)%curve, requested, estimate, variance, present_value)
   if (.not. all(present_value)) error stop 'timepoints presence'
   if (estimate(1) /= 0.0_dp) error stop 'timepoints before first failure'
   print *, 'test_api: PASS'
end program test_api
