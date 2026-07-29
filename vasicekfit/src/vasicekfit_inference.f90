! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Dmitriy Mayorov
module vasicekfit_inference
   use vasicekfit_kinds, only : dp
   use vasicekfit_normal, only : normal_pdf, normal_quantile
   use vasicekfit_linalg, only : invert_matrix, hac_covariance, symmetrize_matrix
   use vasicekfit_model, only : vasicek_fit_result, coefficients
   implicit none
   private

   type, public :: covariance_result
      logical :: ok = .false.
      character(len=:), allocatable :: message
      real(dp), allocatable :: covariance(:,:)
   end type covariance_result

   type, public :: confidence_interval_result
      logical :: ok = .false.
      character(len=:), allocatable :: message
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
   end type confidence_interval_result

   public :: vasicek_covariance, vasicek_confidence_intervals

contains

   function vasicek_covariance(fit, use_hac, hac_lag) result(output)
      type(vasicek_fit_result), intent(in) :: fit
      logical, intent(in), optional :: use_hac
      integer, intent(in), optional :: hac_lag
      type(covariance_result) :: output
      real(dp), allocatable :: base_covariance(:,:), joint_covariance(:,:), jacobian(:,:)
      real(dp), allocatable :: a(:,:), a_inv(:,:), weighted_x(:,:), influence(:,:)
      real(dp) :: sigma2, s, beta0, phi_beta0, var_sigma2
      integer :: n, k, m, j
      logical :: hac, ok

      output%message = ''
      if (.not. fit%ok) then
         output%message = 'cannot calculate covariance for an invalid fit'
         return
      end if
      hac = .false.
      if (present(use_hac)) hac = use_hac
      n = fit%n_observations
      m = fit%n_predictors
      k = m + 1
      sigma2 = fit%sigma2

      allocate(joint_covariance(k+1,k+1))
      joint_covariance = 0.0_dp
      if (.not. hac) then
         joint_covariance(1:k,1:k) = fit%beta_covariance
         if (fit%bias_correct) then
            var_sigma2 = 2.0_dp * sigma2**2 / real(n - m - 1,dp)
         else
            var_sigma2 = 2.0_dp * sigma2**2 * real(n - m - 1,dp) / real(n*n,dp)
         end if
         joint_covariance(k+1,k+1) = var_sigma2
      else
         a = matmul(transpose(fit%design), fit%design) / real(n,dp)
         call invert_matrix(a, a_inv, ok)
         if (.not. ok) then
            output%message = 'could not invert the regression moment matrix'
            return
         end if
         allocate(weighted_x(n,k), influence(n,k+1))
         weighted_x = fit%design
         do j = 1, k
            weighted_x(:,j) = weighted_x(:,j) * fit%residuals
         end do
         influence(:,1:k) = matmul(weighted_x, transpose(a_inv))
         influence(:,k+1) = fit%residuals**2 - sigma2
         if (present(hac_lag)) then
            joint_covariance = hac_covariance(influence, hac_lag)
         else
            joint_covariance = hac_covariance(influence)
         end if
      end if

      s = sqrt(1.0_dp + sigma2)
      beta0 = fit%beta(1)
      phi_beta0 = normal_pdf(beta0 / s)
      allocate(jacobian(k+1,k+1))
      jacobian = 0.0_dp
      jacobian(1,1) = phi_beta0 / s
      jacobian(1,k+1) = -phi_beta0 * beta0 / (2.0_dp * s**3)
      jacobian(2,k+1) = 1.0_dp / (1.0_dp + sigma2)**2
      do j = 1, m
         jacobian(2+j,1+j) = 1.0_dp / s
         jacobian(2+j,k+1) = -fit%beta(1+j) / (2.0_dp * s**3)
      end do

      base_covariance = matmul(jacobian, matmul(joint_covariance, transpose(jacobian)))
      call symmetrize_matrix(base_covariance)
      output%covariance = base_covariance
      output%ok = .true.
   end function vasicek_covariance

   function vasicek_confidence_intervals(fit, level, use_hac, hac_lag) result(output)
      type(vasicek_fit_result), intent(in) :: fit
      real(dp), intent(in), optional :: level
      logical, intent(in), optional :: use_hac
      integer, intent(in), optional :: hac_lag
      type(confidence_interval_result) :: output
      type(covariance_result) :: covariance
      real(dp), allocatable :: estimates(:), standard_errors(:), lower_bounds(:), upper_bounds(:)
      real(dp) :: confidence_level, critical
      integer :: i

      output%message = ''
      confidence_level = 0.95_dp
      if (present(level)) confidence_level = level
      if (confidence_level <= 0.0_dp .or. confidence_level >= 1.0_dp) then
         output%message = 'confidence level must lie strictly inside (0, 1)'
         return
      end if
      if (present(hac_lag)) then
         covariance = vasicek_covariance(fit, use_hac, hac_lag)
      else
         covariance = vasicek_covariance(fit, use_hac)
      end if
      if (.not. covariance%ok) then
         output%message = covariance%message
         return
      end if
      estimates = coefficients(fit)
      allocate(standard_errors(size(estimates)), lower_bounds(size(estimates)), upper_bounds(size(estimates)))
      do i = 1, size(estimates)
         standard_errors(i) = sqrt(max(0.0_dp, covariance%covariance(i,i)))
      end do
      critical = normal_quantile(0.5_dp + 0.5_dp * confidence_level)
      lower_bounds = estimates - critical * standard_errors
      upper_bounds = estimates + critical * standard_errors
      call move_alloc(estimates, output%estimate)
      call move_alloc(standard_errors, output%standard_error)
      call move_alloc(lower_bounds, output%lower)
      call move_alloc(upper_bounds, output%upper)
      output%ok = .true.
   end function vasicek_confidence_intervals

end module vasicekfit_inference
