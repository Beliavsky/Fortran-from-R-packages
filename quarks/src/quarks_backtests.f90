module quarks_backtests
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use quarks_kinds, only : dp
   use quarks_stats, only : chi_square_upper_tail, binomial_cdf, &
      safe_log_probability, nan_value
   use quarks_types, only : coverage_result, traffic_result, loss_result, &
      quarks_ok, quarks_invalid_input, quarks_no_violations
   implicit none
   private

   public :: cvgtest, trftest, lossfun

contains

   function cvgtest(loss, var, p, confidence_level, upstream_formula) result(result)
      real(dp), intent(in) :: loss(:), var(:), p
      real(dp), intent(in), optional :: confidence_level
      logical, intent(in), optional :: upstream_formula
      type(coverage_result) :: result
      logical, allocatable :: hit(:)
      logical :: compatibility
      real(dp) :: puc, p01, p11, log_tp, log_t1, log_t2
      integer :: i, n, n0, n1

      result%p = p
      if (present(confidence_level)) result%confidence_level = confidence_level
      compatibility = .true.
      if (present(upstream_formula)) compatibility = upstream_formula
      n = size(loss)
      if (n <= 1 .or. size(var) /= n .or. p <= 0.0_dp .or. p >= 1.0_dp .or. &
          any(.not. ieee_is_finite(loss)) .or. any(.not. ieee_is_finite(var)) .or. &
          result%confidence_level <= 0.0_dp .or. &
          result%confidence_level >= 1.0_dp) then
         result%status = quarks_invalid_input
         result%message = 'invalid coverage-test inputs'
         result%p_uc = nan_value()
         result%p_ind = nan_value()
         result%p_cc = nan_value()
         return
      end if
      allocate(hit(n))
      hit = loss > var
      n1 = count(hit)
      n0 = n - n1
      result%violations = n1
      if (n1 == 0) then
         result%status = quarks_no_violations
         result%message = 'no VaR violations; tests are not applicable'
         result%p_uc = nan_value()
         result%p_ind = nan_value()
         result%p_cc = nan_value()
         return
      end if
      do i = 1, n - 1
         if (.not. hit(i) .and. .not. hit(i + 1)) result%n00 = result%n00 + 1
         if (.not. hit(i) .and. hit(i + 1)) result%n01 = result%n01 + 1
         if (hit(i) .and. .not. hit(i + 1)) result%n10 = result%n10 + 1
         if (hit(i) .and. hit(i + 1)) result%n11 = result%n11 + 1
      end do
      puc = real(n1, dp) / real(n, dp)
      p01 = real(result%n01, dp) / real(max(1, result%n00 + result%n01), dp)
      p11 = real(result%n11, dp) / real(max(1, result%n10 + result%n11), dp)
      log_tp = safe_log_probability(p, n0) + safe_log_probability(1.0_dp - p, n1)
      log_t1 = safe_log_probability(1.0_dp - puc, n0) + &
         safe_log_probability(puc, n1)
      if (result%n11 == 0) then
         log_t2 = safe_log_probability(1.0_dp - p01, result%n00) + &
            safe_log_probability(p01, result%n01)
      else if (compatibility) then
         log_t2 = safe_log_probability(1.0_dp - p01, result%n00) + &
            safe_log_probability(p01, result%n10) + &
            safe_log_probability(1.0_dp - p11, result%n10) + &
            safe_log_probability(p11, result%n11)
      else
         log_t2 = safe_log_probability(1.0_dp - p01, result%n00) + &
            safe_log_probability(p01, result%n01) + &
            safe_log_probability(1.0_dp - p11, result%n10) + &
            safe_log_probability(p11, result%n11)
      end if
      result%lr_uc = max(0.0_dp, -2.0_dp * (log_tp - log_t1))
      result%lr_ind = max(0.0_dp, -2.0_dp * (log_t1 - log_t2))
      result%lr_cc = max(0.0_dp, -2.0_dp * (log_tp - log_t2))
      result%p_uc = chi_square_upper_tail(result%lr_uc, 1)
      result%p_ind = chi_square_upper_tail(result%lr_ind, 1)
      result%p_cc = chi_square_upper_tail(result%lr_cc, 2)
   end function cvgtest

   function trftest(loss, var, p) result(result)
      real(dp), intent(in) :: loss(:), var(:), p
      type(traffic_result) :: result
      if (size(loss) <= 1 .or. size(var) /= size(loss) .or. &
          p <= 0.0_dp .or. p >= 1.0_dp .or. &
          any(.not. ieee_is_finite(loss)) .or. any(.not. ieee_is_finite(var))) then
         result%status = quarks_invalid_input
         result%message = 'invalid traffic-light-test inputs'
         result%cumulative_probability = nan_value()
         return
      end if
      result%observations = size(loss)
      result%violations = count(var < loss)
      result%breach_probability = 1.0_dp - p
      result%cumulative_probability = binomial_cdf(result%violations, &
         result%observations, result%breach_probability)
   end function trftest

   function lossfun(loss, es, beta) result(result)
      real(dp), intent(in) :: loss(:), es(:)
      real(dp), intent(in), optional :: beta
      type(loss_result) :: result
      real(dp) :: penalty, rlf, flf1, flf2, flf23, flf3
      logical, allocatable :: exceed(:), below(:), nonnegative_below(:), negative_below(:)
      penalty = 1.0e-4_dp
      if (present(beta)) penalty = beta
      if (size(loss) <= 1 .or. size(es) /= size(loss) .or. &
          any(.not. ieee_is_finite(loss)) .or. any(.not. ieee_is_finite(es))) then
         result%status = quarks_invalid_input
         result%message = 'invalid ES loss-function inputs'
         result%lossfun1 = nan_value()
         result%lossfun2 = nan_value()
         result%lossfun3 = nan_value()
         result%lossfun4 = nan_value()
         return
      end if
      allocate(exceed(size(loss)), below(size(loss)), &
         nonnegative_below(size(loss)), negative_below(size(loss)))
      exceed = loss > es
      below = .not. exceed
      nonnegative_below = below .and. loss >= 0.0_dp
      negative_below = below .and. loss < 0.0_dp
      rlf = sum((loss - es)**2, mask=exceed)
      flf1 = penalty * sum(es, mask=below)
      flf2 = penalty * sum(abs(loss - es), mask=below)
      flf23 = penalty * sum(abs(loss - es), mask=nonnegative_below)
      flf3 = penalty * sum(es, mask=negative_below)
      result%lossfun1 = 10000.0_dp * rlf
      result%lossfun2 = 10000.0_dp * (rlf + flf1)
      result%lossfun3 = 10000.0_dp * (rlf + flf2)
      result%lossfun4 = 10000.0_dp * (rlf + flf23 + flf3)
   end function lossfun

end module quarks_backtests
