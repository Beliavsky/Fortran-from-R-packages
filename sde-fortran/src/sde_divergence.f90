! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_divergence
   use sde_kinds, only : dp
   use sde_interfaces, only : state_function, scalar_transform
   use sde_special, only : chi_square_cdf, nan_dp
   use sde_random, only : random_chi_square
   use sde_information, only : dc_log_likelihood
   use sde_optimization, only : optimization_result, nelder_mead_box
   implicit none
   private

   public :: sde_divergence_result
   public :: sde_divergence_test
   public :: fit_and_test_sde_divergence
   public :: phi_negative_log
   public :: phi_kl
   public :: phi_pearson
   public :: phi_hellinger

   type :: sde_divergence_result
      real(dp), allocatable :: theta1(:)
      real(dp), allocatable :: theta0(:)
      real(dp) :: log_likelihood1 = -huge(1.0_dp)
      real(dp) :: log_likelihood0 = -huge(1.0_dp)
      real(dp) :: likelihood_ratio = 0.0_dp
      real(dp) :: divergence = 0.0_dp
      real(dp) :: p_divergence = 0.0_dp
      real(dp) :: p_likelihood_ratio = 0.0_dp
      real(dp) :: c_phi = 0.0_dp
      real(dp) :: k_phi = 0.0_dp
      logical :: p_divergence_simulated = .false.
      logical :: hypotheses_switched = .false.
      type(optimization_result) :: optimizer
   end type sde_divergence_result

contains

   subroutine sde_divergence_test(x, dt, theta1, theta0, drift, diffusion, drift_x, diffusion_x, &
         diffusion_xx, result, phi, c_phi, k_phi, n_simulations)
      real(dp), intent(in) :: x(:), dt, theta1(:), theta0(:)
      procedure(state_function) :: drift, diffusion, drift_x, diffusion_x, diffusion_xx
      type(sde_divergence_result), intent(out) :: result
      procedure(scalar_transform), optional :: phi
      real(dp), intent(in), optional :: c_phi, k_phi
      integer, intent(in), optional :: n_simulations
      real(dp) :: c_value, k_value, difference, coefficient, ratio, statistic
      real(dp) :: chi_value, simulated_statistic
      integer :: i, nsim, exceedances

      if (size(x) < 2) error stop "sde_divergence_test: at least two observations are required"
      if (dt <= 0.0_dp) error stop "sde_divergence_test: dt must be positive"
      if (size(theta1) /= size(theta0)) error stop "sde_divergence_test: parameter dimensions differ"

      c_value = -1.0_dp
      k_value = 1.0_dp
      if (present(c_phi)) c_value = c_phi
      if (present(k_phi)) k_value = k_phi
      nsim = 500000
      if (present(n_simulations)) nsim = max(1, n_simulations)

      allocate(result%theta1(size(theta1)), result%theta0(size(theta0)))
      result%theta1 = theta1
      result%theta0 = theta0
      result%c_phi = c_value
      result%k_phi = k_value
      result%log_likelihood1 = dc_log_likelihood(x, dt, theta1, drift, diffusion, drift_x, &
         diffusion_x, diffusion_xx)
      result%log_likelihood0 = dc_log_likelihood(x, dt, theta0, drift, diffusion, drift_x, &
         diffusion_x, diffusion_xx)

      difference = result%log_likelihood1-result%log_likelihood0
      coefficient = 1.0_dp
      result%hypotheses_switched = difference > 0.0_dp
      if (result%hypotheses_switched) then
         difference = -difference
         coefficient = -1.0_dp
      end if

      ratio = exp(max(log(tiny(1.0_dp)), difference))
      result%likelihood_ratio = -2.0_dp*difference
      if (present(phi)) then
         result%divergence = phi(ratio)
      else
         result%divergence = phi_negative_log(ratio)
      end if
      result%p_likelihood_ratio = chi_square_cdf(result%likelihood_ratio, real(size(theta0), dp), &
         lower_tail=.false.)

      result%p_divergence = nan_dp()
      result%p_divergence_simulated = .false.
      if (abs(k_value) > epsilon(1.0_dp)) then
         exceedances = 0
         do i = 1, nsim
            chi_value = random_chi_square(real(size(theta0), dp))
            if (abs(c_value) <= epsilon(1.0_dp)) then
               simulated_statistic = 0.5_dp*k_value*chi_value*chi_value
            else
               simulated_statistic = coefficient*0.5_dp*c_value*chi_value+ &
                  0.5_dp*(c_value+k_value)*chi_value*chi_value
            end if
            if (simulated_statistic > result%divergence) exceedances = exceedances+1
         end do
         result%p_divergence = real(exceedances, dp)/real(nsim, dp)
         result%p_divergence_simulated = .true.
      else if (abs(c_value) > epsilon(1.0_dp)) then
         statistic = coefficient*result%divergence/c_value
         result%p_divergence = chi_square_cdf(statistic, real(size(theta0), dp), lower_tail=.false.)
      end if
   end subroutine sde_divergence_test

   subroutine fit_and_test_sde_divergence(x, dt, theta0, initial, drift, diffusion, drift_x, &
         diffusion_x, diffusion_xx, result, phi, c_phi, k_phi, n_simulations, lower, upper, &
         max_iterations)
      real(dp), intent(in) :: x(:), dt, theta0(:), initial(:)
      procedure(state_function) :: drift, diffusion, drift_x, diffusion_x, diffusion_xx
      type(sde_divergence_result), intent(out) :: result
      procedure(scalar_transform), optional :: phi
      real(dp), intent(in), optional :: c_phi, k_phi, lower(:), upper(:)
      integer, intent(in), optional :: n_simulations, max_iterations
      type(optimization_result) :: opt

      if (size(initial) /= size(theta0)) then
         error stop "fit_and_test_sde_divergence: initial and null parameter dimensions differ"
      end if
      call nelder_mead_box(euler_objective, initial, opt, lower=lower, upper=upper, &
         max_iterations=max_iterations)
      call sde_divergence_test(x, dt, opt%x, theta0, drift, diffusion, drift_x, diffusion_x, &
         diffusion_xx, result, phi=phi, c_phi=c_phi, k_phi=k_phi, n_simulations=n_simulations)
      result%optimizer = opt

   contains

      function euler_objective(theta) result(value)
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
         real(dp) :: sigma_value, variance, residual
         integer :: j

         value = 0.0_dp
         do j = 1, size(x)-1
            sigma_value = diffusion(x(j), theta)
            variance = dt*sigma_value*sigma_value
            if (variance <= tiny(1.0_dp)) then
               value = huge(1.0_dp)/16.0_dp
               return
            end if
            residual = x(j+1)-x(j)-dt*drift(x(j), theta)
            value = value+0.5_dp*(log(2.0_dp*acos(-1.0_dp)*variance)+residual*residual/variance)
         end do
      end function euler_objective

   end subroutine fit_and_test_sde_divergence

   pure function phi_negative_log(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (x <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = -log(x)
      end if
   end function phi_negative_log

   pure function phi_kl(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (x <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = x*log(x)-x+1.0_dp
      end if
   end function phi_kl

   pure function phi_pearson(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = 0.5_dp*(x-1.0_dp)*(x-1.0_dp)
   end function phi_pearson

   pure function phi_hellinger(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      if (x < 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = 2.0_dp*(sqrt(x)-1.0_dp)*(sqrt(x)-1.0_dp)
      end if
   end function phi_hellinger

end module sde_divergence
