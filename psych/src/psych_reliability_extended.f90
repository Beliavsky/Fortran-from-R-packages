module psych_reliability_extended
   use r_kinds, only : dp
   use psych_correlations, only : squared_multiple_correlations
   implicit none
   private
   public :: icc_complete_matrix, omega_statistics
contains
   pure subroutine icc_complete_matrix(x, icc, mean_squares)
      real(dp), intent(in) :: x(:, :) !! Complete subject-by-rater/item measurement matrix.
      real(dp), intent(out) :: icc(6) !! ICC1, ICC2, ICC3 and corresponding average-rater ICC1k, ICC2k, ICC3k.
      real(dp), intent(out) :: mean_squares(4) !! Subject, rater, residual, and pooled-within mean squares.
      real(dp), allocatable :: row_mean(:), col_mean(:)
      real(dp) :: grand, ss_subject, ss_rater, ss_total, ss_error, msb, msj, mse, msw
      integer :: i, j, n, k
      n = size(x, 1)
      k = size(x, 2)
      allocate(row_mean(n), col_mean(k))
      row_mean = sum(x, dim=2) / real(max(k, 1), dp)
      col_mean = sum(x, dim=1) / real(max(n, 1), dp)
      grand = sum(x) / real(max(n * k, 1), dp)
      ss_subject = real(k, dp) * sum((row_mean - grand) ** 2)
      ss_rater = real(n, dp) * sum((col_mean - grand) ** 2)
      ss_total = sum((x - grand) ** 2)
      ss_error = max(ss_total - ss_subject - ss_rater, 0.0_dp)
      if (n > 1) then
         msb = ss_subject / real(n - 1, dp)
      else
         msb = 0.0_dp
      end if
      if (k > 1) then
         msj = ss_rater / real(k - 1, dp)
      else
         msj = 0.0_dp
      end if
      if (n > 1 .and. k > 1) then
         mse = ss_error / real((n - 1) * (k - 1), dp)
      else
         mse = 0.0_dp
      end if
      if (n * (k - 1) > 0) then
         msw = (ss_rater + ss_error) / real(n * (k - 1), dp)
      else
         msw = 0.0_dp
      end if
      mean_squares = [msb, msj, mse, msw]
      icc(1) = safe_ratio(msb - msw, msb + real(k - 1, dp) * msw)
      icc(2) = safe_ratio(msb - mse, msb + real(k - 1, dp) * mse + real(k, dp) * (msj - mse) / real(max(n, 1), dp))
      icc(3) = safe_ratio(msb - mse, msb + real(k - 1, dp) * mse)
      icc(4) = safe_ratio(msb - msw, msb)
      icc(5) = safe_ratio(msb - mse, msb + (msj - mse) / real(max(n, 1), dp))
      icc(6) = safe_ratio(msb - mse, msb)
      i = n
      j = k
      if (i + j < 0) icc = 0.0_dp
   end subroutine icc_complete_matrix

   subroutine omega_statistics(r, loadings, omega_h, omega_total, omega_limit, alpha, lambda6)
      real(dp), intent(in) :: r(:, :) !! Item covariance or correlation matrix.
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loadings with the general factor in column one.
      real(dp), intent(out) :: omega_h !! General-factor omega hierarchical estimate.
      real(dp), intent(out) :: omega_total !! Reliability due to all supplied common factors.
      real(dp), intent(out) :: omega_limit !! General-factor variance relative to common-score variance.
      real(dp), intent(out) :: alpha !! Cronbach alpha computed from r.
      real(dp), intent(out) :: lambda6 !! Guttman lambda-6 estimate from squared multiple correlations.
      real(dp), allocatable :: smc(:)
      real(dp) :: total_variance, item_variance, general_sum_sq, uniqueness, common_variance, sum_smc
      integer :: i, p, info
      p = size(r, 1)
      total_variance = sum(r)
      item_variance = 0.0_dp
      uniqueness = 0.0_dp
      do i = 1, p
         item_variance = item_variance + r(i, i)
         uniqueness = uniqueness + max(r(i, i) - sum(loadings(i, :) ** 2), 0.0_dp)
      end do
      general_sum_sq = sum(loadings(:, 1)) ** 2
      common_variance = total_variance - uniqueness
      omega_h = safe_ratio(general_sum_sq, total_variance)
      omega_total = safe_ratio(common_variance, total_variance)
      omega_limit = safe_ratio(general_sum_sq, common_variance)
      if (p > 1) then
         alpha = real(p, dp) / real(p - 1, dp) * safe_ratio(total_variance - item_variance, total_variance)
      else
         alpha = 0.0_dp
      end if
      allocate(smc(p))
      call squared_multiple_correlations(r, smc, info)
      if (info /= 0) smc = 0.0_dp
      sum_smc = sum(smc)
      lambda6 = safe_ratio(total_variance + sum_smc - item_variance, total_variance)
   end subroutine omega_statistics

   pure real(dp) function safe_ratio(numerator, denominator) result(value)
      real(dp), intent(in) :: numerator !! Numerator of a guarded scalar ratio.
      real(dp), intent(in) :: denominator !! Denominator; near-zero values yield zero.
      if (abs(denominator) > tiny(1.0_dp)) then
         value = numerator / denominator
      else
         value = 0.0_dp
      end if
   end function safe_ratio
end module psych_reliability_extended
