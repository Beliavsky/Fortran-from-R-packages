! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
program sde_demo
   use sde
   implicit none

   real(dp), parameter :: dt = 1.0_dp/252.0_dp
   real(dp), parameter :: theta(3) = [1.2_dp, 0.8_dp, 0.45_dp]
   real(dp), allocatable :: paths(:, :), observations(:)
   type(kernel_estimate) :: drift_estimate
   real(dp) :: exact_loglik, euler_loglik

   call seed_rng(20260725_i64)
   call simulate_ou_exact([0.5_dp], dt, 500, theta, paths)
   observations = paths(:, 1)

   exact_loglik = ou_log_likelihood(observations, dt, theta)
   euler_loglik = euler_log_likelihood(observations, dt, theta, ou_drift, ou_diffusion)
   call kernel_drift(observations, dt, drift_estimate, n_grid=25)

   write(*, '(a)') "sde-fortran demonstration"
   write(*, '(a, i0)') "observations:             ", size(observations)
   write(*, '(a, f12.6)') "first value:              ", observations(1)
   write(*, '(a, f12.6)') "last value:               ", observations(size(observations))
   write(*, '(a, f12.4)') "exact OU log likelihood:  ", exact_loglik
   write(*, '(a, f12.4)') "Euler log likelihood:     ", euler_loglik
   write(*, '(a, f12.6)') "kernel drift at midpoint: ", &
      drift_estimate%y((size(drift_estimate%y)+1)/2)

contains

   pure function ou_drift(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = local_theta(1)-local_theta(2)*x+0.0_dp*t
   end function ou_drift

   pure function ou_diffusion(t, x, local_theta) result(value)
      real(dp), intent(in) :: t, x, local_theta(:)
      real(dp) :: value
      value = local_theta(3)+0.0_dp*(t+x)
   end function ou_diffusion

end program sde_demo
