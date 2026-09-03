! Modern Fortran translation of cmprsk computational code.
! Upstream cmprsk Copyright (C) 2000 Robert Gray.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_crr_fit
   use r_kinds, only : dp
   use cmprsk, only : crr_result, fit_crr, predict_crr, cmprsk_success
   implicit none

   integer, parameter :: n = 12
   real(dp), parameter :: tol = 2.0e-9_dp
   real(dp), allocatable :: cif(:, :)
   real(dp) :: cov1(n, 2)
   real(dp) :: cov2(n, 1)
   real(dp) :: profile1(2, 2)
   real(dp) :: profile2(2, 1)
   real(dp) :: tf(4, 1)
   real(dp) :: time(n)
   integer :: event_code(n)
   integer :: group(n)
   integer :: i
   integer :: order(n)
   integer :: status
   type(crr_result) :: fit

   time = [5.0_dp, 1.0_dp, 8.0_dp, 3.0_dp, 6.0_dp, 1.0_dp, 7.0_dp, 4.0_dp, 2.0_dp, 8.0_dp, 5.0_dp, 3.0_dp]
   event_code = [1, 1, 0, 1, 2, 0, 1, 2, 2, 2, 0, 0]
   group = [1, 1, 1, 2, 1, 2, 2, 2, 1, 2, 2, 1]
   order = [7, 1, 11, 4, 9, 2, 10, 6, 3, 12, 8, 5]
   do i = 1, n
      cov1(i, 1) = (-1.0_dp)**order(i)*0.2_dp*real(order(i), dp)
      cov1(i, 2) = real(mod(order(i), 3), dp)*0.35_dp - 0.2_dp
      cov2(i, 1) = 0.1_dp*real(order(i), dp) - 0.4_dp
   end do
   tf(:, 1) = [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]

   call fit_crr(time, event_code, cov1, fit, status, cov2=cov2, time_functions=tf, &
                censor_group=group, gtol=1.0e-8_dp, maxiter=20)
   if (status /= cmprsk_success .or. .not. fit%converged) error stop 'fit_crr status'
   if (maxval(abs(fit%coefficients - [0.18003069023334850_dp, 1.6361297654789915_dp, &
                                       -1.2506868951378158_dp])) > tol) error stop 'fit_crr coefficients'
   if (abs(fit%loglik + 8.0022561401616557_dp) > tol) error stop 'fit_crr loglik'
   if (abs(fit%loglik_null + 8.6524231406763423_dp) > tol) error stop 'fit_crr null loglik'
   if (abs(fit%variance(2, 2) - 1.5881258806342764_dp) > tol) error stop 'fit_crr variance'
   if (abs(fit%baseline_jump(4) - 0.45144666233904740_dp) > tol) error stop 'fit_crr baseline'

   profile1 = reshape([0.0_dp, 0.0_dp, 0.5_dp, -0.2_dp], [2, 2])
   profile2(:, 1) = [0.0_dp, 0.3_dp]
   call predict_crr(fit, profile1, cif, status, cov2=profile2)
   if (status /= cmprsk_success) error stop 'predict_crr status'
   if (abs(cif(4, 1) - 0.84765396525219661_dp) > tol) error stop 'predict_crr first profile'
   if (abs(cif(4, 2) - 0.28142290550861304_dp) > tol) error stop 'predict_crr second profile'
   print *, 'test_crr_fit: PASS'
end program test_crr_fit
