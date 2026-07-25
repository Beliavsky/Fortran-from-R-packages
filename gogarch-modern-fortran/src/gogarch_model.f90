! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch_model
   use gogarch_kinds, only : dp
   use gogarch_linalg, only : determinant_matrix, symmetric_invsqrt
   use gogarch_core, only : initialize_gogarch, build_covariance_path
   use gogarch_univariate, only : fit_univariate, forecast_univariate
   use gogarch_distributions, only : random_innovation
   use gogarch_types, only : gogarch_fit, univariate_spec
   implicit none
   private
   public :: build_gogarch_fit, forecast_gogarch, standardized_residuals
   public :: simulate_fitted_gogarch, reconstruction_error, factor_coefficients
   public :: factor_coefficients_full

contains

   subroutine build_gogarch_fit(data, rotation, method, fit, max_garch_iterations, optimizer_iterations, parameters, &
      weights, factor_spec)
      real(dp), intent(in) :: data(:,:), rotation(:,:)
      character(len=*), intent(in) :: method
      type(gogarch_fit), intent(out) :: fit
      integer, intent(in), optional :: max_garch_iterations, optimizer_iterations
      real(dp), intent(in), optional :: parameters(:), weights(:)
      type(univariate_spec), intent(in), optional :: factor_spec
      real(dp) :: det_mixing, logdet
      logical :: ok, det_ok
      integer :: n, m, j, garch_maxit
      type(univariate_spec) :: spec
      n = size(data,1)
      m = size(data,2)
      garch_maxit = 500
      if (present(max_garch_iterations)) garch_maxit = max_garch_iterations
      spec = univariate_spec()
      if (present(factor_spec)) spec = factor_spec
      fit%n = n
      fit%m = m
      fit%method = method
      fit%factor_spec = spec
      allocate(fit%data(n,m),fit%sample_covariance(m,m),fit%covariance_sqrt(m,m),fit%covariance_invsqrt(m,m))
      allocate(fit%eigenvectors(m,m),fit%eigenvalues(m),fit%rotation(m,m),fit%mixing(m,m),fit%factors(n,m))
      allocate(fit%factor_models(m),fit%factor_variance(n,m),fit%covariance(m,m,n))
      fit%data = data
      fit%rotation = rotation
      call initialize_gogarch(data,fit%sample_covariance,fit%eigenvectors,fit%eigenvalues,fit%covariance_sqrt, &
         fit%covariance_invsqrt,ok)
      if (.not. ok) then
         fit%status = 2
         fit%log_likelihood = -huge(1.0_dp)
         fit%mixing = 0.0_dp
         fit%factors = 0.0_dp
         fit%factor_variance = 0.0_dp
         fit%covariance = 0.0_dp
         return
      end if
      fit%mixing = matmul(fit%covariance_sqrt,fit%rotation)
      fit%factors = matmul(matmul(data,fit%covariance_invsqrt),fit%rotation)
      fit%log_likelihood = 0.0_dp
      fit%status = 0
      do j = 1, m
         fit%factor_models(j) = fit_univariate(fit%factors(:,j),spec,garch_maxit)
         if (.not. allocated(fit%factor_models(j)%variance)) then
            fit%status = 3
            fit%factor_variance(:,j) = 1.0_dp
         else
            fit%factor_variance(:,j) = fit%factor_models(j)%variance
            fit%log_likelihood = fit%log_likelihood+fit%factor_models(j)%log_likelihood
            if (fit%factor_models(j)%status > 1) fit%status = 3
         end if
      end do
      det_mixing = determinant_matrix(fit%mixing,det_ok)
      if (det_ok .and. abs(det_mixing) > tiny(1.0_dp)) then
         logdet = log(abs(det_mixing))
         fit%log_likelihood = fit%log_likelihood-real(n,dp)*logdet
      else
         fit%status = 4
      end if
      call build_covariance_path(fit%mixing,fit%factor_variance,fit%covariance)
      if (present(optimizer_iterations)) fit%optimizer_iterations = optimizer_iterations
      if (present(parameters)) then
         allocate(fit%objective_parameters(size(parameters)))
         fit%objective_parameters = parameters
      end if
      if (present(weights)) then
         allocate(fit%mm_weights(size(weights)))
         fit%mm_weights = weights
      end if
   end subroutine build_gogarch_fit

   subroutine forecast_gogarch(fit, n_ahead, mean_forecast, covariance_forecast, factor_variance_forecast)
      type(gogarch_fit), intent(in) :: fit
      integer, intent(in) :: n_ahead
      real(dp), intent(out) :: mean_forecast(n_ahead,fit%m)
      real(dp), intent(out) :: covariance_forecast(fit%m,fit%m,n_ahead)
      real(dp), intent(out), optional :: factor_variance_forecast(n_ahead,fit%m)
      real(dp) :: factor_mean(n_ahead,fit%m), factor_var(n_ahead,fit%m)
      real(dp) :: scaled(fit%m,fit%m)
      integer :: j, k
      do j = 1, fit%m
         call forecast_univariate(fit%factor_models(j),n_ahead,factor_mean(:,j),factor_var(:,j))
      end do
      mean_forecast = matmul(factor_mean,transpose(fit%mixing))
      do k = 1, n_ahead
         scaled = fit%mixing
         do j = 1, fit%m
            scaled(:,j) = scaled(:,j)*factor_var(k,j)
         end do
         covariance_forecast(:,:,k) = matmul(scaled,transpose(fit%mixing))
         covariance_forecast(:,:,k) = 0.5_dp*(covariance_forecast(:,:,k)+transpose(covariance_forecast(:,:,k)))
      end do
      if (present(factor_variance_forecast)) factor_variance_forecast = factor_var
   end subroutine forecast_gogarch

   subroutine standardized_residuals(fit, residuals, ok)
      type(gogarch_fit), intent(in) :: fit
      real(dp), intent(out) :: residuals(fit%n,fit%m)
      logical, intent(out), optional :: ok
      real(dp) :: invroot(fit%m,fit%m), factor_mean(fit%m), asset_mean(fit%m), centered(fit%m)
      logical :: step_ok, all_ok
      integer :: t, j
      all_ok = .true.
      do j = 1, fit%m
         factor_mean(j) = fit%factor_models(j)%mean
      end do
      asset_mean = matmul(fit%mixing,factor_mean)
      do t = 1, fit%n
         invroot = symmetric_invsqrt(fit%covariance(:,:,t),1.0e-12_dp,step_ok)
         if (step_ok) then
            centered = fit%data(t,:)-asset_mean
            residuals(t,:) = matmul(invroot,centered)
         else
            residuals(t,:) = 0.0_dp
            all_ok = .false.
         end if
      end do
      if (present(ok)) ok = all_ok
   end subroutine standardized_residuals

   subroutine simulate_fitted_gogarch(fit, n_sim, data, factor_data, factor_variance)
      type(gogarch_fit), intent(in) :: fit
      integer, intent(in) :: n_sim
      real(dp), intent(out) :: data(n_sim,fit%m)
      real(dp), intent(out), optional :: factor_data(n_sim,fit%m), factor_variance(n_sim,fit%m)
      real(dp), allocatable :: factors(:,:), variances(:,:), eps_all(:), power_all(:)
      real(dp) :: gamma_i
      integer :: i, j, lag, n, index, p, q
      allocate(factors(n_sim,fit%m),variances(n_sim,fit%m))
      do j = 1, fit%m
         n = size(fit%factor_models(j)%variance)
         p = fit%factor_models(j)%p
         q = fit%factor_models(j)%q
         allocate(eps_all(n+n_sim),power_all(n+n_sim))
         eps_all(1:n) = fit%factor_models(j)%residuals
         power_all(1:n) = fit%factor_models(j)%power_scale
         do i = 1, n_sim
            index = n+i
            power_all(index) = fit%factor_models(j)%omega
            do lag = 1, p
               gamma_i = 0.0_dp
               if (lag <= fit%factor_models(j)%o) gamma_i = fit%factor_models(j)%leverage(lag)
               power_all(index) = power_all(index)+fit%factor_models(j)%arch(lag)* &
                  (abs(eps_all(index-lag))-gamma_i*eps_all(index-lag))**fit%factor_models(j)%delta
            end do
            do lag = 1, q
               power_all(index) = power_all(index)+fit%factor_models(j)%garch(lag)*power_all(index-lag)
            end do
            power_all(index) = max(power_all(index),1.0e-14_dp)
            eps_all(index) = power_all(index)**(1.0_dp/fit%factor_models(j)%delta)* &
               random_innovation(fit%factor_models(j)%distribution,fit%factor_models(j)%shape, &
               fit%factor_models(j)%skew)
            factors(i,j) = fit%factor_models(j)%mean+eps_all(index)
            variances(i,j) = power_all(index)**(2.0_dp/fit%factor_models(j)%delta)
         end do
         deallocate(eps_all,power_all)
      end do
      data = matmul(factors,transpose(fit%mixing))
      if (present(factor_data)) factor_data = factors
      if (present(factor_variance)) factor_variance = variances
   end subroutine simulate_fitted_gogarch

   pure function reconstruction_error(fit) result(error)
      type(gogarch_fit), intent(in) :: fit
      real(dp) :: error
      error = maxval(abs(fit%data-matmul(fit%factors,transpose(fit%mixing))))
   end function reconstruction_error

   subroutine factor_coefficients(fit, coefficients)
      type(gogarch_fit), intent(in) :: fit
      real(dp), intent(out) :: coefficients(fit%m,4)
      integer :: j
      do j = 1, fit%m
         coefficients(j,:) = [fit%factor_models(j)%mean,fit%factor_models(j)%omega, &
            fit%factor_models(j)%alpha,fit%factor_models(j)%beta]
      end do
   end subroutine factor_coefficients

   subroutine factor_coefficients_full(fit, means, omegas, arch, leverage, garch, delta, shape, skew)
      type(gogarch_fit), intent(in) :: fit
      real(dp), intent(out) :: means(fit%m), omegas(fit%m)
      real(dp), intent(out) :: arch(fit%m,fit%factor_spec%p)
      real(dp), intent(out) :: leverage(fit%m,fit%factor_spec%o)
      real(dp), intent(out) :: garch(fit%m,fit%factor_spec%q)
      real(dp), intent(out) :: delta(fit%m), shape(fit%m), skew(fit%m)
      integer :: j
      do j = 1, fit%m
         means(j) = fit%factor_models(j)%mean
         omegas(j) = fit%factor_models(j)%omega
         arch(j,:) = fit%factor_models(j)%arch
         if (fit%factor_spec%o > 0) leverage(j,:) = fit%factor_models(j)%leverage
         if (fit%factor_spec%q > 0) garch(j,:) = fit%factor_models(j)%garch
         delta(j) = fit%factor_models(j)%delta
         shape(j) = fit%factor_models(j)%shape
         skew(j) = fit%factor_models(j)%skew
      end do
   end subroutine factor_coefficients_full

end module gogarch_model
