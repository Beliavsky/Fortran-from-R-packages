! Modern Fortran translation of cmprsk computational code.
! Upstream cmprsk Copyright (C) 2000 Robert Gray.
! SPDX-License-Identifier: GPL-2.0-or-later
program test_crr_kernels
   use r_kinds, only : dp
   use cmprsk_censoring, only : censoring_survival_left
   use cmprsk_crr_kernels, only : crr_objective_score_info, crr_variance_kernel, &
                                  crr_score_residuals_kernel, crr_baseline_jumps_kernel
   use cmprsk, only : cmprsk_success
   implicit none

   integer, parameter :: n = 12
   integer, parameter :: np = 3
   integer, parameter :: ndf = 4
   real(dp), parameter :: tol = 2.0e-12_dp
   real(dp) :: beta(np)
   real(dp) :: h(np, np)
   real(dp) :: meat(np, np)
   real(dp) :: objective
   real(dp) :: residuals(np, ndf)
   real(dp) :: score(np)
   real(dp) :: time(n)
   real(dp) :: tf(ndf, 1)
   real(dp) :: variance_info(np, np)
   real(dp) :: wt(2, n)
   real(dp) :: x(n, 2)
   real(dp) :: x2(n, 1)
   real(dp) :: jump(ndf)
   integer :: event_code(n)
   integer :: group(n)
   integer :: i
   integer :: status
   logical :: censored(n)

   time = [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 8.0_dp]
   event_code = [1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2]
   group = [1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2]
   do i = 1, n
      x(i, 1) = (-1.0_dp)**i*0.2_dp*real(i, dp)
      x(i, 2) = real(mod(i, 3), dp)*0.35_dp - 0.2_dp
      x2(i, 1) = 0.1_dp*real(i, dp) - 0.4_dp
   end do
   tf(:, 1) = [0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
   beta = [0.15_dp, -0.25_dp, 0.1_dp]
   censored = event_code == 0
   call censoring_survival_left(time, censored, group, wt, status)
   if (status /= cmprsk_success) error stop 'censoring survival status'

   call crr_objective_score_info(time, event_code, x, x2, tf, wt, group, beta, objective, score, h)
   if (abs(objective - 8.8199617344466148_dp) > tol) error stop 'crr objective'
   if (maxval(abs(score - [0.99951776173152496_dp, -0.41884840286835345_dp, 0.69906045233095226_dp])) > tol) &
      error stop 'crr score'
   if (abs(h(1, 1) - 10.655279200658777_dp) > tol) error stop 'crr Hessian'

   call crr_variance_kernel(time, event_code, x, x2, tf, wt, group, beta, variance_info, meat)
   if (maxval(abs(variance_info - h)) > tol) error stop 'crr variance information'
   if (abs(meat(1, 1) - 7.2258912161691100_dp) > tol) error stop 'crr sandwich meat'

   call crr_score_residuals_kernel(time, event_code, x, x2, tf, wt, group, beta, residuals)
   if (abs(residuals(1, 3) + 2.0714294492677015_dp) > tol) error stop 'crr residual'
   call crr_baseline_jumps_kernel(time, event_code, x, x2, tf, wt, group, beta, jump)
   if (maxval(abs(jump - [0.081631102685078949_dp, 0.094360528000667229_dp, &
                           0.11182417309649641_dp, 0.14364440552120028_dp])) > tol) error stop 'crr jumps'
   print *, 'test_crr_kernels: PASS'
end program test_crr_kernels
