! SPDX-License-Identifier: GPL-2.0-or-later
module pbkrtest_tests
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq, r_pf
   use r_special, only : r_regularized_gamma_q
   use pbkrtest_types, only : bootstrap_result_t, lrt_result_t, pbkr_invalid_argument, &
      pbkr_success
   implicit none
   private
   public :: bootstrap_p_values
   public :: likelihood_ratio_test

contains

   pure subroutine likelihood_ratio_test(loglik_full, loglik_small, nparam_full, nparam_small, result, status)
      real(dp), intent(in) :: loglik_full !! Maximized log-likelihood of the larger model on the comparison scale.
      real(dp), intent(in) :: loglik_small !! Maximized log-likelihood of the nested smaller model on the same scale.
      integer, intent(in) :: nparam_full !! Number of estimated parameters in the larger model.
      integer, intent(in) :: nparam_small !! Number of estimated parameters in the smaller model.
      type(lrt_result_t), intent(out) :: result !! Likelihood-ratio statistic, parameter-count difference, and chi-square p-value.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.

      result%statistic = 2.0_dp * (loglik_full - loglik_small)
      result%df = nparam_full - nparam_small
      if (result%df <= 0) then
         status = pbkr_invalid_argument
         result%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      result%p_value = r_pchisq(result%statistic, real(result%df, dp), lower_tail=.false.)
      status = pbkr_success
   end subroutine likelihood_ratio_test

   pure subroutine bootstrap_p_values(lrt_statistic, ndf, reference, result, status)
      real(dp), intent(in) :: lrt_statistic !! Observed likelihood-ratio statistic.
      integer, intent(in) :: ndf !! Chi-square numerator degrees of freedom for the model comparison.
      real(dp), intent(in) :: reference(:) !! Parametric-bootstrap reference statistics; only positive values enter moment fits.
      type(bootstrap_result_t), intent(out) :: result !! Bootstrap, Bartlett, gamma, and moment-matched F calibration summaries.
      integer, intent(out) :: status !! `pbkr_success` on success or a package error code.
      real(dp), allocatable :: positive(:)
      real(dp) :: ee2
      real(dp) :: nan_value
      real(dp) :: se_raw
      integer :: i
      integer :: k

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      result%nsim = size(reference)
      if (ndf <= 0 .or. result%nsim <= 0) then
         status = pbkr_invalid_argument
         return
      end if
      result%npos = count(reference > 0.0_dp)
      allocate(positive(result%npos))
      k = 0
      do i = 1, result%nsim
         if (reference(i) > 0.0_dp) then
            k = k + 1
            positive(k) = reference(i)
         end if
      end do
      result%n_extreme = count(positive > lrt_statistic)
      result%p_chisq = r_pchisq(lrt_statistic, real(ndf, dp), lower_tail=.false.)
      result%p_bootstrap_all = real(1 + result%n_extreme, dp) / real(1 + result%nsim, dp)
      result%f_statistic = lrt_statistic / real(ndf, dp)

      if (result%npos == 0) then
         result%mean_positive = nan_value
         result%variance_positive = nan_value
         result%p_bootstrap = nan_value
         result%standard_error = nan_value
         result%ci_low = nan_value
         result%ci_high = nan_value
         result%bartlett_statistic = nan_value
         result%p_bartlett = nan_value
         result%gamma_scale = nan_value
         result%gamma_shape = nan_value
         result%p_gamma = nan_value
         result%f_ddf = nan_value
         result%p_f = nan_value
         status = pbkr_success
         return
      end if

      result%mean_positive = sum(positive) / real(result%npos, dp)
      result%p_bootstrap = real(1 + result%n_extreme, dp) / real(1 + result%npos, dp)
      se_raw = sqrt(result%p_bootstrap * (1.0_dp - result%p_bootstrap) / real(result%npos, dp))
      result%standard_error = round_four(se_raw)
      result%ci_low = round_four(result%p_bootstrap - 1.96_dp * result%standard_error)
      result%ci_high = round_four(result%p_bootstrap + 1.96_dp * result%standard_error)
      result%bartlett_statistic = real(ndf, dp) * lrt_statistic / result%mean_positive
      result%p_bartlett = r_pchisq(result%bartlett_statistic, real(ndf, dp), lower_tail=.false.)

      if (result%npos > 1) then
         result%variance_positive = sum((positive - result%mean_positive) ** 2) / &
            real(result%npos - 1, dp)
      else
         result%variance_positive = nan_value
      end if
      if (result%npos > 1 .and. result%variance_positive > 0.0_dp) then
         result%gamma_scale = result%variance_positive / result%mean_positive
         result%gamma_shape = result%mean_positive ** 2 / result%variance_positive
         result%p_gamma = r_regularized_gamma_q(result%gamma_shape, &
            lrt_statistic / result%gamma_scale)
      else
         result%gamma_scale = nan_value
         result%gamma_shape = nan_value
         result%p_gamma = nan_value
      end if

      ee2 = result%mean_positive / real(ndf, dp)
      if (abs(ee2 - 1.0_dp) <= epsilon(1.0_dp)) then
         result%f_ddf = nan_value
         result%p_f = nan_value
      else
         result%f_ddf = 2.0_dp * ee2 / (ee2 - 1.0_dp)
         if (result%f_ddf > 2.0_dp) then
            result%p_f = r_pf(result%f_statistic, real(ndf, dp), result%f_ddf, lower_tail=.false.)
         else
            result%p_f = nan_value
         end if
      end if
      status = pbkr_success
   end subroutine bootstrap_p_values

   pure elemental real(dp) function round_four(x) result(value)
      real(dp), intent(in) :: x !! Value rounded to four decimal places using nearest-integer arithmetic.

      value = anint(10000.0_dp * x) / 10000.0_dp
   end function round_four

end module pbkrtest_tests
