! SPDX-License-Identifier: BSD-3-Clause
! Numerical translation of posterior 1.7.0, R/pit.R, R/discrete-summaries.R,
! and R/uniformity_test.R. Randomized PIT uses caller-supplied uniforms.
module posterior_predictive
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   use r_kinds, only: dp, r_pi
   use r_ordering, only: r_order, r_sort_values_in_place
   use r_special, only: r_regularized_beta
   use posterior_diagnostics, only: is_constant
   use posterior_pareto, only: gpdfit, pgeneralized_pareto, ps_tail_length
   implicit none
   private
   public :: entropy, dissent, modal_category, pit, pareto_pit, uniformity_test
contains
   pure function nan_value() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_value

   pure function entropy(x) result(value)
      real(dp), intent(in) :: x(:) !! Numeric category labels; missing values propagate, and observed categories define normalization.
      real(dp) :: value, probability
      integer :: n, categories, i, multiplicity
      logical :: first(size(x))
      value = nan_value()
      n = size(x)
      if (any(ieee_is_nan(x))) return
      first = .true.
      do i = 1, n
         first(i) = .not. any(x(:i - 1) == x(i))
      end do
      categories = count(first)
      value = 0.0_dp
      if (categories <= 1) return
      do i = 1, n
         if (.not. first(i)) cycle
         multiplicity = count(x == x(i))
         probability = real(multiplicity, dp) / real(n, dp)
         value = value - probability * log(probability)
      end do
      value = value / log(real(categories, dp))
   end function entropy

   pure function dissent(x) result(value)
      real(dp), intent(in) :: x(:) !! Finite numeric ordinal scores; range defines category distance, empty sample returns zero.
      real(dp) :: value, distance, center
      value = nan_value()
      if (any(.not. ieee_is_finite(x))) return
      value = 0.0_dp
      if (size(x) == 0) return
      distance = maxval(x) - minval(x)
      if (distance == 0.0_dp) return
      center = sum(x) / real(size(x), dp)
      value = -sum(log(1.0_dp - abs(x - center) / distance)) / real(size(x), dp) / log(2.0_dp)
   end function dissent

   pure function modal_category(x) result(value)
      real(dp), intent(in) :: x(:) !! Numeric category labels; ties select the first category encountered in the sample.
      real(dp) :: value
      integer :: i, frequency, best
      value = nan_value()
      if (size(x) == 0 .or. any(ieee_is_nan(x))) return
      best = 0
      do i = 1, size(x)
         frequency = count(x == x(i))
         if (frequency > best) then
            value = x(i)
            best = frequency
         end if
      end do
   end function modal_category

   pure subroutine normalize_weights(n, weights, log_weights, normalized, status)
      integer, intent(in) :: n !! Number of predictive draws to weight.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative linear weights or log weights; absent means equal weights.
      logical, intent(in), optional :: log_weights !! Interpret supplied weights on log scale, default false.
      real(dp), intent(out) :: normalized(n) !! Linear nonnegative probabilities summing to one on success.
      integer, intent(out) :: status !! Zero success, one invalid weights or empty input.
      logical :: logarithm
      real(dp) :: maximum
      normalized = nan_value()
      status = 1
      if (n < 1) return
      logarithm = .false.
      if (present(log_weights)) logarithm = log_weights
      if (.not. present(weights)) then
         normalized = 1.0_dp / real(n, dp)
         status = 0
         return
      end if
      if (size(weights) /= n .or. any(ieee_is_nan(weights))) return
      if (logarithm) then
         maximum = maxval(weights)
         if (.not. ieee_is_finite(maximum)) return
         normalized = exp(weights - maximum)
      else
         if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
         maximum = maxval(weights)
         if (maximum <= 0.0_dp) return
         normalized = weights / maximum
      end if
      normalized = normalized / sum(normalized)
      status = 0
   end subroutine normalize_weights

   pure function pit(x, y, uniform, weights, log_weights) result(value)
      real(dp), intent(in) :: x(:) !! Predictive draws for one observed scalar, allowing infinities but no NaNs.
      real(dp), intent(in) :: y !! Observed scalar, compared strictly below draws and within any tie mass.
      real(dp), intent(in), optional :: uniform !! Caller-supplied U[0,1] tie randomizer; default 0.5 gives mid-PIT.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative importance weights matching predictive draws.
      logical, intent(in), optional :: log_weights !! Whether weights are logarithms, default false.
      real(dp) :: value, u, normalized(size(x))
      integer :: status
      value = nan_value()
      if (any(ieee_is_nan(x)) .or. ieee_is_nan(y)) return
      u = 0.5_dp
      if (present(uniform)) u = uniform
      if (.not. ieee_is_finite(u) .or. u < 0.0_dp .or. u > 1.0_dp) return
      call normalize_weights(size(x), weights, log_weights, normalized, status)
      if (status /= 0) return
      value = min(1.0_dp, sum(normalized, mask = x < y) + u * sum(normalized, mask = x == y))
   end function pit

   pure function pareto_pit(x, y, uniform, weights, log_weights, ndraws_tail) result(value)
      real(dp), intent(in) :: x(:) !! Predictive draws for a single observed variable.
      real(dp), intent(in) :: y !! Observed scalar at which to evaluate the smoothed predictive CDF.
      real(dp), intent(in), optional :: uniform !! Supplied U[0,1] tie randomizer, default 0.5.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative importance weights matching draws.
      logical, intent(in), optional :: log_weights !! Interpret weights as log weights; default false.
      integer, intent(in), optional :: ndraws_tail !! Tail fitting count, default ps_tail_length(N,1), capped at floor(N/2).
      real(dp) :: value, normalized(size(x)), sorted(size(x)), fitted_weights(size(x))
      real(dp) :: cutoff, k, sigma, proportion, cdf, bound
      integer, allocatable :: order(:)
      integer :: n, nt, status, side
      logical :: replaced
      value = pit(x, y, uniform, weights, log_weights)
      if (ieee_is_nan(value)) return
      n = size(x)
      call normalize_weights(n, weights, log_weights, normalized, status)
      nt = ps_tail_length(n, 1.0_dp)
      if (present(ndraws_tail)) nt = ndraws_tail
      if (nt >= 5 .and. all(ieee_is_finite(x)) .and. .not. is_constant(x)) then
         nt = min(nt, n / 2)
         replaced = .false.
         if (nt >= 2 .and. nt < n) then
            do side = 1, 2
               if (replaced) exit
               if (side == 1) then
                  sorted = x
               else
                  sorted = -x
               end if
               call r_order(sorted, order)
               sorted = sorted(order)
               fitted_weights = normalized(order)
               if (is_constant(sorted(n - nt + 1:))) cycle
               cutoff = sorted(n - nt)
               if (cutoff == sorted(n - nt + 1)) cutoff = cutoff - epsilon(1.0_dp)
               if (present(weights)) then
                  call gpdfit(sorted(n - nt + 1:) - cutoff, k, sigma, sort_x = .false., &
                     weights = fitted_weights(n - nt + 1:), status = status)
               else
                  call gpdfit(sorted(n - nt + 1:) - cutoff, k, sigma, sort_x = .false., status = status)
               end if
               if (status /= 0 .or. .not. ieee_is_finite(k)) cycle
               proportion = sum(fitted_weights(n - nt + 1:))
               if (side == 1 .and. y > cutoff) then
                  cdf = pgeneralized_pareto(y, mu = cutoff, sigma = sigma, k = k)
                  value = 1.0_dp - proportion * (1.0_dp - cdf)
                  replaced = .true.
               else if (side == 2 .and. -y > cutoff) then
                  cdf = pgeneralized_pareto(-y, mu = cutoff, sigma = sigma, k = k)
                  value = proportion * (1.0_dp - cdf)
               end if
            end do
         end if
      end if
      bound = 1.0_dp / real(n, dp) / 10000.0_dp
      value = min(1.0_dp - bound, max(bound, value))
   end function pareto_pit

   pure function binomial_cdf(k, n, p) result(value)
      integer, intent(in) :: k !! Integer count upper bound.
      integer, intent(in) :: n !! Number of binomial trials, nonnegative.
      real(dp), intent(in) :: p !! Success probability in [0,1].
      real(dp) :: value
      if (k < 0) then
         value = 0.0_dp
      else if (k >= n) then
         value = 1.0_dp
      else if (p == 0.0_dp) then
         value = 1.0_dp
      else if (p == 1.0_dp) then
         value = 0.0_dp
      else
         value = r_regularized_beta(1.0_dp - p, real(n - k, dp), real(k + 1, dp))
      end if
   end function binomial_cdf

   pure function negative_cauchy_quantile(p) result(value)
      real(dp), intent(in) :: p !! Cauchy cumulative probability in [0,1].
      real(dp) :: value
      if (p == 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else if (p == 1.0_dp) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else if (p < 0.5_dp) then
         value = 1.0_dp / tan(r_pi * p)
      else if (p > 0.5_dp) then
         value = -1.0_dp / tan(r_pi * (1.0_dp - p))
      else
         value = 0.0_dp
      end if
   end function negative_cauchy_quantile

   pure subroutine uniformity_test(values, method, pvalue, pointwise, status)
      real(dp), intent(in) :: values(:) !! PIT values in [0,1], finite and nonempty.
      character(*), intent(in) :: method !! One of 'POT', 'PRIT', or 'PIET', case sensitive.
      real(dp), intent(out) :: pvalue !! Combined Cauchy-test p-value, NaN when the test cannot be evaluated.
      real(dp), intent(out) :: pointwise(size(values)) !! Shapley contributions in sorted-PIT order, matching upstream.
      integer, intent(out), optional :: status !! Zero success, one invalid inputs, two no p-values below 0.5 for truncation.
      real(dp) :: sorted(size(values)), p(size(values)), cauchy(size(values)), cdf, combined, harmonic
      integer :: i, n, count_le, included
      logical :: truncate
      pvalue = nan_value()
      pointwise = nan_value()
      if (present(status)) status = 1
      n = size(values)
      if (n == 0 .or. any(.not. ieee_is_finite(values))) return
      if (any(values < 0.0_dp) .or. any(values > 1.0_dp)) return
      if (method /= 'POT' .and. method /= 'PRIT' .and. method /= 'PIET') return
      sorted = values
      call r_sort_values_in_place(sorted)
      do i = 1, n
         select case (method)
         case ('POT')
            cdf = r_regularized_beta(sorted(i), real(i, dp), real(n - i + 1, dp))
            p(i) = 2.0_dp * min(cdf, 1.0_dp - cdf)
         case ('PRIT')
            count_le = count(sorted <= sorted(i))
            p(i) = 2.0_dp * min(1.0_dp - binomial_cdf(count_le - 1, n, sorted(i)), &
               binomial_cdf(count_le, n, sorted(i)))
         case ('PIET')
            p(i) = 2.0_dp * min(sorted(i), 1.0_dp - sorted(i))
         end select
      end do
      truncate = method /= 'PIET'
      included = n
      if (truncate) included = count(p < 0.5_dp)
      if (included == 0) then
         if (present(status)) status = 2
         return
      end if
      combined = 0.0_dp
      do i = 1, n
         if (truncate .and. p(i) >= 0.5_dp) cycle
         combined = combined + negative_cauchy_quantile(p(i))
      end do
      combined = combined / real(included, dp)
      if (ieee_is_nan(combined)) then
         pvalue = nan_value()
      else if (combined > 1.0_dp) then
         pvalue = atan(1.0_dp / combined) / r_pi
      else if (combined < -1.0_dp) then
         pvalue = 1.0_dp - atan(-1.0_dp / combined) / r_pi
      else
         pvalue = 0.5_dp - atan(combined) / r_pi
      end if
      cauchy = tan((0.5_dp - p) * r_pi)
      if (n == 1) then
         pointwise = 0.0_dp
      else
         harmonic = sum([(1.0_dp / real(i, dp), i = 1, n)])
         do i = 1, n
            pointwise(i) = cauchy(i) / real(n, dp) + (harmonic - 1.0_dp) / real(n, dp) * &
               (cauchy(i) - (sum(cauchy(:i - 1)) + sum(cauchy(i + 1:))) / real(n - 1, dp))
         end do
      end if
      if (present(status)) status = 0
   end subroutine uniformity_test
end module posterior_predictive
