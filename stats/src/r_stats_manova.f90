! SPDX-License-Identifier: MIT
! SPDX-FileComment: Multivariate linear models and MANOVA tests corresponding to R stats.
module r_stats_manova
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use r_distributions, only: r_pf
   use r_kinds, only: dp
   use r_linalg, only: general_real_eigenvalues, inverse_matrix, least_squares
   use r_optional, only: optval
   use r_stats_types, only: manova_result_t, mlm_fit_t
   implicit none
   private

   public :: manova_lm, mlm_fit

contains

   function mlm_fit(y, xpred, intercept) result(fit)
      !! Fits a full-rank multivariate linear model by shared dense least squares.
      real(dp), intent(in) :: y(:, :) !! Response matrix with shape `(n,m)`.
      real(dp), intent(in) :: xpred(:, :) !! Predictor matrix `(n,p)` without an intercept.
      logical, intent(in), optional :: intercept !! Include an intercept; defaults to true.
      type(mlm_fit_t) :: fit
      real(dp), allocatable :: design(:, :)
      integer :: info, k, n, p

      n = size(y, 1)
      p = size(xpred, 2)
      fit%has_intercept = optval(intercept, .true.)
      k = p + merge(1, 0, fit%has_intercept)
      if (size(xpred, 1) /= n .or. size(y, 2) < 2 .or. n <= k .or. &
          .not. all(ieee_is_finite(y)) .or. .not. all(ieee_is_finite(xpred))) then
         fit%status = 1
         return
      end if
      allocate (design(n, k))
      if (fit%has_intercept) then
         design(:, 1) = 1.0_dp
         if (p > 0) design(:, 2:) = xpred
      else
         design = xpred
      end if
      allocate (fit%coefficients(k, size(y, 2)))
      call least_squares(design, y, fit%coefficients, info)
      if (info /= 0) then
         fit%status = info
         return
      end if
      fit%y = y
      fit%fitted = matmul(design, fit%coefficients)
      fit%residuals = y - fit%fitted
      fit%residual_sscp = matmul(transpose(fit%residuals), fit%residuals)
      fit%rank = k
      fit%df_residual = n - k
   end function mlm_fit

   function manova_lm(reduced, full) result(result)
      !! Computes R-compatible MANOVA criteria for two nested multivariate fits.
      type(mlm_fit_t), intent(in) :: reduced !! Smaller model fitted to the same responses.
      type(mlm_fit_t), intent(in) :: full !! Larger model fitted to the same responses.
      type(manova_result_t) :: result
      real(dp), allocatable :: error_inverse(:, :), hypothesis(:, :), imag(:), matrix(:, :)
      real(dp), allocatable :: real_eigenvalues(:)
      real(dp) :: m, n_value, p, q, s, temporary1, temporary2, temporary3
      integer :: info

      if (.not. allocated(reduced%y) .or. .not. allocated(full%y)) then
         result%status = 1
         return
      end if
      if (any(shape(reduced%y) /= shape(full%y)) .or. reduced%rank >= full%rank .or. &
          full%df_residual <= 0 .or. any(abs(reduced%y - full%y) > 0.0_dp)) then
         result%status = 1
         return
      end if
      hypothesis = reduced%residual_sscp - full%residual_sscp
      call inverse_matrix(full%residual_sscp, error_inverse, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      matrix = matmul(error_inverse, hypothesis)
      call general_real_eigenvalues(matrix, real_eigenvalues, imag, info)
      if (info /= 0 .or. any(abs(imag) > 1.0e-8_dp)) then
         result%status = merge(info, 1, info /= 0)
         return
      end if
      result%eigenvalues = max(0.0_dp, real_eigenvalues)
      result%hypothesis_df = full%rank - reduced%rank
      p = real(size(full%y, 2), dp)
      q = real(result%hypothesis_df, dp)
      s = min(p, q)
      m = 0.5_dp*(abs(p - q) - 1.0_dp)
      n_value = 0.5_dp*(real(full%df_residual, dp) - p - 1.0_dp)

      result%statistic(1) = sum(result%eigenvalues/(1.0_dp + result%eigenvalues))
      temporary1 = 2.0_dp*m + s + 1.0_dp
      temporary2 = 2.0_dp*n_value + s + 1.0_dp
      result%approximate_f(1) = temporary2/temporary1*result%statistic(1)/ &
                                (s - result%statistic(1))
      result%numerator_df(1) = s*temporary1
      result%denominator_df(1) = s*temporary2

      result%statistic(2) = product(1.0_dp/(1.0_dp + result%eigenvalues))
      temporary1 = real(full%df_residual, dp) - 0.5_dp*(p - q + 1.0_dp)
      temporary2 = (p*q - 2.0_dp)/4.0_dp
      temporary3 = p**2 + q**2 - 5.0_dp
      if (temporary3 > 0.0_dp) then
         temporary3 = sqrt(((p*q)**2 - 4.0_dp)/temporary3)
      else
         temporary3 = 1.0_dp
      end if
      result%approximate_f(2) = (result%statistic(2)**(-1.0_dp/temporary3) - 1.0_dp)* &
                                (temporary1*temporary3 - 2.0_dp*temporary2)/(p*q)
      result%numerator_df(2) = p*q
      result%denominator_df(2) = temporary1*temporary3 - 2.0_dp*temporary2

      result%statistic(3) = sum(result%eigenvalues)
      temporary1 = 2.0_dp*m + s + 1.0_dp
      temporary2 = 2.0_dp*(s*n_value + 1.0_dp)
      result%approximate_f(3) = temporary2*result%statistic(3)/(s*s*temporary1)
      result%numerator_df(3) = s*temporary1
      result%denominator_df(3) = temporary2

      result%statistic(4) = maxval(result%eigenvalues)
      temporary1 = max(p, q)
      temporary2 = real(full%df_residual, dp) - temporary1 + q
      result%approximate_f(4) = temporary2*result%statistic(4)/temporary1
      result%numerator_df(4) = temporary1
      result%denominator_df(4) = temporary2
      do info = 1, 4
         result%p_value(info) = r_pf(result%approximate_f(info), result%numerator_df(info), &
                                     result%denominator_df(info), lower_tail=.false.)
      end do
   end function manova_lm

end module r_stats_manova
