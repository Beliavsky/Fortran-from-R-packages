! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_gmm
   use sde_kinds, only : dp
   use sde_interfaces, only : moment_function
   use sde_optimization, only : optimization_result, nelder_mead_box
   use sde_linalg, only : invert_matrix, quadratic_form, numerical_hessian
   implicit none
   private

   public :: gmm_result
   public :: fit_gmm
   public :: gmm_moment_mean
   public :: gmm_hac_covariance

   type :: gmm_result
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: moment_mean(:)
      real(dp), allocatable :: weighting(:, :)
      real(dp), allocatable :: hessian(:, :)
      real(dp) :: objective = huge(1.0_dp)
      integer :: stage2_iterations = 0
      logical :: converged = .false.
      type(optimization_result) :: first_stage
      type(optimization_result) :: final_stage
   end type gmm_result

contains

   subroutine gmm_moment_mean(x, dt, theta, n_moments, moments, mean_moments)
      real(dp), intent(in) :: x(:), dt, theta(:)
      integer, intent(in) :: n_moments
      procedure(moment_function) :: moments
      real(dp), intent(out) :: mean_moments(:)
      real(dp), allocatable :: values(:)
      integer :: i, n_transitions

      n_transitions = size(x)-1
      if (n_transitions < 1 .or. dt <= 0.0_dp .or. size(mean_moments) /= n_moments) then
         error stop "gmm_moment_mean: invalid dimensions"
      end if
      allocate(values(n_moments))
      mean_moments = 0.0_dp
      do i = 1, n_transitions
         call moments(x(i+1), x(i), theta, dt, values)
         mean_moments = mean_moments+values
      end do
      mean_moments = mean_moments/real(n_transitions, dp)
   end subroutine gmm_moment_mean

   subroutine gmm_hac_covariance(x, dt, theta, n_moments, moments, covariance, max_lag)
      real(dp), intent(in) :: x(:), dt, theta(:)
      integer, intent(in) :: n_moments
      procedure(moment_function) :: moments
      real(dp), intent(out) :: covariance(:, :)
      integer, intent(in), optional :: max_lag
      real(dp), allocatable :: u(:, :), values(:), gamma(:, :)
      real(dp) :: weight
      integer :: n_transitions, lag_limit, i, lag

      n_transitions = size(x)-1
      if (n_transitions < 2 .or. size(covariance, 1) /= n_moments .or. &
          size(covariance, 2) /= n_moments) then
         error stop "gmm_hac_covariance: invalid dimensions"
      end if
      lag_limit = n_transitions-1
      if (present(max_lag)) lag_limit = min(max_lag, n_transitions-1)
      lag_limit = max(0, lag_limit)
      allocate(u(n_moments, n_transitions), values(n_moments), gamma(n_moments, n_moments))
      do i = 1, n_transitions
         call moments(x(i+1), x(i), theta, dt, values)
         u(:, i) = values
      end do
      covariance = matmul(u, transpose(u))/real(size(x), dp)
      do lag = 1, lag_limit
         gamma = matmul(u(:, lag+1:n_transitions), transpose(u(:, 1:n_transitions-lag)))/real(size(x), dp)
         weight = 1.0_dp-real(lag, dp)/real(lag_limit+1, dp)
         covariance = covariance+weight*(gamma+transpose(gamma))
      end do
   end subroutine gmm_hac_covariance

   subroutine fit_gmm(x, dt, n_moments, moments, initial, result, lower, upper, max_stage2_iterations, &
         parameter_tolerance, objective_tolerance, max_lag, optimizer_iterations)
      real(dp), intent(in) :: x(:), dt, initial(:)
      integer, intent(in) :: n_moments
      procedure(moment_function) :: moments
      type(gmm_result), intent(out) :: result
      real(dp), intent(in), optional :: lower(:), upper(:)
      integer, intent(in), optional :: max_stage2_iterations, max_lag, optimizer_iterations
      real(dp), intent(in), optional :: parameter_tolerance, objective_tolerance
      real(dp), allocatable :: theta_old(:), theta_new(:), covariance(:, :), weighting(:, :), mean_values(:)
      real(dp) :: tol_parameter, tol_objective, ridge
      integer :: stage_limit, stage, lag_limit, opt_limit, inverse_status, n_parameters
      type(optimization_result) :: opt

      if (size(x) < 3 .or. dt <= 0.0_dp .or. n_moments <= 0 .or. size(initial) == 0) then
         error stop "fit_gmm: invalid input"
      end if
      n_parameters = size(initial)
      stage_limit = 30
      tol_parameter = 1.0e-3_dp
      tol_objective = 1.0e-3_dp
      lag_limit = size(x)-2
      opt_limit = 1500
      if (present(max_stage2_iterations)) stage_limit = max_stage2_iterations
      if (present(parameter_tolerance)) tol_parameter = parameter_tolerance
      if (present(objective_tolerance)) tol_objective = objective_tolerance
      if (present(max_lag)) lag_limit = max_lag
      if (present(optimizer_iterations)) opt_limit = optimizer_iterations
      allocate(theta_old(n_parameters), theta_new(n_parameters), covariance(n_moments, n_moments), &
         weighting(n_moments, n_moments), mean_values(n_moments))
      weighting = 0.0_dp
      weighting = identity_matrix(n_moments)

      call nelder_mead_box(unweighted_objective, initial, opt, lower=lower, upper=upper, &
         max_iterations=opt_limit)
      result%first_stage = opt
      theta_old = opt%x
      result%converged = .false.

      do stage = 1, stage_limit
         call gmm_hac_covariance(x, dt, theta_old, n_moments, moments, covariance, lag_limit)
         ridge = 0.0_dp
         do
            call invert_matrix(covariance, weighting, inverse_status, ridge)
            if (inverse_status == 0) exit
            if (ridge <= 0.0_dp) then
               ridge = 1.0e-10_dp*max(1.0_dp, maxval(abs(covariance)))
            else
               ridge = 10.0_dp*ridge
            end if
            if (ridge > 1.0e-2_dp*max(1.0_dp, maxval(abs(covariance)))) then
               error stop "fit_gmm: HAC covariance could not be regularized"
            end if
         end do
         call nelder_mead_box(weighted_objective, theta_old, opt, lower=lower, upper=upper, &
            max_iterations=opt_limit)
         theta_new = opt%x
         result%final_stage = opt
         if (sum(abs(theta_new-theta_old)) < tol_parameter .or. opt%value < tol_objective) then
            result%converged = .true.
            theta_old = theta_new
            exit
         end if
         theta_old = theta_new
      end do

      result%stage2_iterations = min(stage, stage_limit)
      allocate(result%estimate(n_parameters), result%moment_mean(n_moments), &
         result%weighting(n_moments, n_moments), result%hessian(n_parameters, n_parameters))
      result%estimate = theta_old
      result%weighting = weighting
      call gmm_moment_mean(x, dt, result%estimate, n_moments, moments, result%moment_mean)
      result%objective = quadratic_form(result%moment_mean, result%weighting)
      call numerical_hessian(weighted_objective, result%estimate, result%hessian)

   contains

      function unweighted_objective(theta) result(value)
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
         call gmm_moment_mean(x, dt, theta, n_moments, moments, mean_values)
         value = sum(mean_values**2)
      end function unweighted_objective

      function weighted_objective(theta) result(value)
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
         call gmm_moment_mean(x, dt, theta, n_moments, moments, mean_values)
         value = quadratic_form(mean_values, weighting)
      end function weighted_objective

   end subroutine fit_gmm

   pure function identity_matrix(n) result(matrix)
      integer, intent(in) :: n
      real(dp) :: matrix(n, n)
      integer :: i
      matrix = 0.0_dp
      do i = 1, n
         matrix(i, i) = 1.0_dp
      end do
   end function identity_matrix

end module sde_gmm
