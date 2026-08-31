! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! Multilevel R-squared calculations from mitml multilevelR2.R.
module mitml_r2
   use r_kinds, only : dp
   use mitml_numeric, only : mean_real, sample_variance, trace_matrix
   use mitml_types, only : MITML_ERR_ARGUMENT, MITML_ERR_DIMENSION, MITML_OK, multilevel_r2_result
   implicit none
   private

   public :: multilevel_r2
   public :: intraclass_correlation

contains

   subroutine multilevel_r2(x, beta, mu_z, cov_z, random_intercept_var, random_intercept_slope_cov, &
      random_slope_cov, residual_var, result, null_random_intercept_var, null_residual_var)
      real(dp), intent(in) :: x(:, :) !! Fixed-effect design excluding intercept, observations by coefficients.
      real(dp), intent(in) :: beta(:) !! Fixed-effect coefficients corresponding to columns of x.
      real(dp), intent(in) :: mu_z(:) !! Means of random-slope covariates; length q, possibly zero.
      real(dp), intent(in) :: cov_z(:, :) !! Covariance matrix of random-slope covariates, shape q by q.
      real(dp), intent(in) :: random_intercept_var !! Alternative-model random-intercept variance t0.1.
      real(dp), intent(in) :: random_intercept_slope_cov(:) !! Covariances between random intercept and q slopes.
      real(dp), intent(in) :: random_slope_cov(:, :) !! Alternative-model random-slope covariance matrix t11.1.
      real(dp), intent(in) :: residual_var !! Alternative-model residual variance s1.
      type(multilevel_r2_result), intent(out) :: result !! RB1, RB2, SB, and MVP statistics.
      real(dp), intent(in), optional :: null_random_intercept_var !! Null-model random-intercept variance t0.0.
      real(dp), intent(in), optional :: null_residual_var !! Null-model residual variance s0.
      real(dp), allocatable :: fitted(:)
      real(dp) :: total_variance
      real(dp) :: vyhat
      integer :: q

      call clear_r2_result(result)
      q = size(mu_z)
      if (size(x, 2) /= size(beta)) then
         call set_r2_error(result, MITML_ERR_DIMENSION, "x and beta dimensions do not match")
         return
      end if
      if (size(random_intercept_slope_cov) /= q) then
         call set_r2_error(result, MITML_ERR_DIMENSION, "random intercept-slope covariance length mismatch")
         return
      end if
      if (size(cov_z, 1) /= q .or. size(cov_z, 2) /= q) then
         call set_r2_error(result, MITML_ERR_DIMENSION, "cov_z must be q by q")
         return
      end if
      if (size(random_slope_cov, 1) /= q .or. size(random_slope_cov, 2) /= q) then
         call set_r2_error(result, MITML_ERR_DIMENSION, "random_slope_cov must be q by q")
         return
      end if
      if (random_intercept_var < 0.0_dp .or. residual_var < 0.0_dp) then
         call set_r2_error(result, MITML_ERR_ARGUMENT, "variance components must be nonnegative")
         return
      end if
      if (present(null_random_intercept_var) .neqv. present(null_residual_var)) then
         call set_r2_error(result, MITML_ERR_ARGUMENT, "both null-model variances must be supplied together")
         return
      end if

      allocate(fitted(size(x, 1)))
      if (size(beta) == 0) then
         fitted = 0.0_dp
         vyhat = 0.0_dp
      else
         fitted = matmul(x, beta)
         vyhat = sample_variance(fitted)
      end if
      total_variance = vyhat + random_intercept_var + residual_var
      if (q > 0) then
         total_variance = total_variance + 2.0_dp * dot_product(mu_z, random_intercept_slope_cov)
         total_variance = total_variance + dot_product(mu_z, matmul(random_slope_cov, mu_z))
         total_variance = total_variance + trace_matrix(matmul(random_slope_cov, cov_z))
      end if
      if (total_variance <= 0.0_dp) then
         call set_r2_error(result, MITML_ERR_ARGUMENT, "total model variance must be positive")
         return
      end if
      result%mvp = vyhat / total_variance

      if (present(null_random_intercept_var)) then
         if (null_random_intercept_var <= 0.0_dp .or. null_residual_var <= 0.0_dp) then
            call set_r2_error(result, MITML_ERR_ARGUMENT, "null-model variance components must be positive")
            return
         end if
         result%rb1 = 1.0_dp - residual_var / null_residual_var
         result%rb2 = 1.0_dp - random_intercept_var / null_random_intercept_var
         result%sb = 1.0_dp - (residual_var + random_intercept_var) / &
            (null_residual_var + null_random_intercept_var)
         result%has_reduction_measures = .true.
      end if
      result%status = MITML_OK
      result%message = "ok"
   end subroutine multilevel_r2

   pure elemental real(dp) function intraclass_correlation(random_intercept_var, residual_var) result(value)
      real(dp), intent(in) :: random_intercept_var !! Nonnegative random-intercept variance.
      real(dp), intent(in) :: residual_var !! Nonnegative residual variance.

      value = random_intercept_var / (random_intercept_var + residual_var)
   end function intraclass_correlation

   subroutine clear_r2_result(result)
      type(multilevel_r2_result), intent(out) :: result !! R-squared result reset before calculation.

      result%rb1 = 0.0_dp
      result%rb2 = 0.0_dp
      result%sb = 0.0_dp
      result%mvp = 0.0_dp
      result%has_reduction_measures = .false.
      result%status = MITML_OK
      result%message = ""
   end subroutine clear_r2_result

   subroutine set_r2_error(result, status, message)
      type(multilevel_r2_result), intent(inout) :: result !! Result object receiving an error status and message.
      integer, intent(in) :: status !! MITML status code for the failure.
      character(len=*), intent(in) :: message !! Human-readable failure explanation.

      result%status = status
      result%message = message
   end subroutine set_r2_error

end module mitml_r2
