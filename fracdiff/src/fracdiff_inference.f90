! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_inference
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use fracdiff_kinds, only : dp
   use fracdiff_status, only : fd_ok, fd_invalid_input, fd_singular_hessian
   use fracdiff_linalg, only : invert_matrix
   use fracdiff_optimize, only : evaluate_fixed_likelihood
   implicit none
   private

   public :: numerical_likelihood_hessian, covariance_from_hessian
   public :: inverse_normal_cdf

contains

   subroutine numerical_likelihood_hessian(x, m_terms, d, ar, ma, h, hessian, status)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: m_terms
      real(dp), intent(in) :: d, ar(:), ma(:), h
      real(dp), intent(out) :: hessian(:,:)
      integer, intent(out) :: status

      real(dp), allocatable :: theta(:), theta_work(:), ar_work(:), ma_work(:)
      real(dp) :: f0, fpp, fpm, fmp, fmm, variance, estimated_mean
      real(dp) :: hi, hj
      integer :: p, q, npar, i, j, local_status

      p = size(ar)
      q = size(ma)
      npar = p + q + 1
      status = fd_ok
      hessian = 0.0_dp
      if (h <= 0.0_dp .or. size(hessian,1) /= npar .or. size(hessian,2) /= npar) then
         status = fd_invalid_input
         return
      end if

      allocate(theta(npar), theta_work(npar), ar_work(p), ma_work(q))
      theta(1) = d
      if (p > 0) theta(2:p+1) = ar
      if (q > 0) theta(p+2:npar) = ma
      call likelihood_from_theta(theta, f0, local_status)
      if (local_status /= fd_ok) then
         status = local_status
         return
      end if

      do i = 1, npar
         hi = h
         theta_work = theta
         theta_work(i) = theta(i) + hi
         call likelihood_from_theta(theta_work, fpp, local_status)
         if (local_status /= fd_ok) then
            status = local_status
            return
         end if
         theta_work(i) = theta(i) - hi
         call likelihood_from_theta(theta_work, fmm, local_status)
         if (local_status /= fd_ok) then
            status = local_status
            return
         end if
         hessian(i,i) = (fpp - 2.0_dp*f0 + fmm)/(hi*hi)

         do j = i + 1, npar
            hj = h
            theta_work = theta
            theta_work(i) = theta(i) + hi
            theta_work(j) = theta(j) + hj
            call likelihood_from_theta(theta_work, fpp, local_status)
            if (local_status /= fd_ok) then
               status = local_status
               return
            end if
            theta_work = theta
            theta_work(i) = theta(i) + hi
            theta_work(j) = theta(j) - hj
            call likelihood_from_theta(theta_work, fpm, local_status)
            if (local_status /= fd_ok) then
               status = local_status
               return
            end if
            theta_work = theta
            theta_work(i) = theta(i) - hi
            theta_work(j) = theta(j) + hj
            call likelihood_from_theta(theta_work, fmp, local_status)
            if (local_status /= fd_ok) then
               status = local_status
               return
            end if
            theta_work = theta
            theta_work(i) = theta(i) - hi
            theta_work(j) = theta(j) - hj
            call likelihood_from_theta(theta_work, fmm, local_status)
            if (local_status /= fd_ok) then
               status = local_status
               return
            end if
            hessian(i,j) = (fpp - fpm - fmp + fmm)/(4.0_dp*hi*hj)
            hessian(j,i) = hessian(i,j)
         end do
      end do

   contains

      subroutine likelihood_from_theta(parameters, value, eval_status)
         real(dp), intent(in) :: parameters(:)
         real(dp), intent(out) :: value
         integer, intent(out) :: eval_status

         if (p > 0) ar_work = parameters(2:p+1)
         if (q > 0) ma_work = parameters(p+2:npar)
         if (parameters(1) <= -0.5_dp .or. parameters(1) >= 0.5_dp) then
            eval_status = fd_invalid_input
            value = -huge(1.0_dp)
            return
         end if
         call evaluate_fixed_likelihood(x, parameters(1), m_terms, ar_work, ma_work, &
                                        value, variance, estimated_mean, eval_status)
      end subroutine likelihood_from_theta

   end subroutine numerical_likelihood_hessian

   subroutine covariance_from_hessian(hessian, covariance, standard_error, correlation, status)
      real(dp), intent(in) :: hessian(:,:)
      real(dp), intent(out) :: covariance(:,:), standard_error(:), correlation(:,:)
      integer, intent(out) :: status

      real(dp), allocatable :: information(:,:)
      real(dp) :: nan_value
      integer :: n, i, j, inverse_status

      n = size(hessian,1)
      status = fd_ok
      if (size(hessian,2) /= n .or. size(covariance,1) /= n .or. size(covariance,2) /= n .or. &
          size(standard_error) /= n .or. size(correlation,1) /= n .or. size(correlation,2) /= n) then
         status = fd_invalid_input
         return
      end if
      allocate(information(n,n))
      information = -0.5_dp*(hessian + transpose(hessian))
      call invert_matrix(information, covariance, inverse_status)
      if (inverse_status /= 0) then
         status = fd_singular_hessian
         nan_value = ieee_value(1.0_dp, ieee_quiet_nan)
         covariance = nan_value
         standard_error = nan_value
         correlation = nan_value
         return
      end if

      do i = 1, n
         if (covariance(i,i) <= 0.0_dp) then
            status = fd_singular_hessian
            standard_error(i) = 0.0_dp
         else
            standard_error(i) = sqrt(covariance(i,i))
         end if
      end do
      correlation = 0.0_dp
      do j = 1, n
         do i = 1, n
            if (standard_error(i) > 0.0_dp .and. standard_error(j) > 0.0_dp) then
               correlation(i,j) = covariance(i,j)/(standard_error(i)*standard_error(j))
            end if
         end do
      end do
   end subroutine covariance_from_hessian

   pure function inverse_normal_cdf(probability) result(value)
      real(dp), intent(in) :: probability
      real(dp) :: value
      real(dp) :: q, r
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
         4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: dcoef(4) = [ &
         7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
         2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: p_low = 0.02425_dp
      real(dp), parameter :: p_high = 1.0_dp - p_low

      if (probability <= 0.0_dp) then
         value = -huge(1.0_dp)
      else if (probability >= 1.0_dp) then
         value = huge(1.0_dp)
      else if (probability < p_low) then
         q = sqrt(-2.0_dp*log(probability))
         value = (((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
                 ((((dcoef(1)*q + dcoef(2))*q + dcoef(3))*q + dcoef(4))*q + 1.0_dp)
      else if (probability <= p_high) then
         q = probability - 0.5_dp
         r = q*q
         value = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6))*q / &
                 (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + 1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp - probability))
         value = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
                  ((((dcoef(1)*q + dcoef(2))*q + dcoef(3))*q + dcoef(4))*q + 1.0_dp)
      end if
   end function inverse_normal_cdf

end module fracdiff_inference
