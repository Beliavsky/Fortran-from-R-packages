! SPDX-License-Identifier: BSD-3-Clause
! Numerical translation of posterior 1.7.0, R/convergence.R and R/nested_rhat.R.
module posterior_diagnostics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf
   use r_kinds, only: dp
   use r_distributions, only: r_qnorm
   use r_quantiles, only: r_quantile_type7
   use r_ordering, only: r_order, r_sort_values_in_place
   use r_time_series, only: r_autocovariance
   use r_special, only: r_regularized_beta
   implicit none
   private
   public :: dp, rhat, rhat_basic, rhat_nested, ess_basic, ess_bulk, ess_tail
   public :: ess_mean, ess_sd, ess_quantile, ess_median
   public :: mcse_mean, mcse_sd, mcse_quantile, mcse_median
   public :: autocovariance, autocorrelation, r_scale, u_scale, z_scale
   public :: is_constant, quantile2, summarise_draws
contains
   pure function nan_value() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_value

   pure function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:) !! Sample observations; at least two are required.
      real(dp) :: value, center
      value = nan_value()
      if (size(x) < 2) return
      center = sum(x) / real(size(x), dp)
      value = sum((x - center)**2) / real(size(x) - 1, dp)
   end function sample_variance

   pure function is_constant(x, tol) result(value)
      real(dp), intent(in) :: x(:) !! Numeric values whose overall range is checked.
      real(dp), intent(in), optional :: tol !! Finite nonnegative range tolerance, default machine epsilon; invalid tolerance returns false.
      logical :: value
      real(dp) :: tolerance, lowest, highest
      tolerance = epsilon(1.0_dp)
      if (present(tol)) tolerance = tol
      value = .false.
      if (.not. ieee_is_finite(tolerance)) return
      if (tolerance < 0.0_dp) return
      if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) return
      lowest = minval(x)
      highest = maxval(x)
      if (highest > 0.0_dp .and. lowest < 0.0_dp) then
         ! A mathematical range above the largest finite real cannot be below tol.
         if (highest > huge(1.0_dp) + lowest) return
      end if
      value = highest - lowest < tolerance
   end function is_constant

   pure function invalid_draws(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains for a single scalar variable.
      logical :: value
      value = .true.
      if (size(x) == 0 .or. any(.not. ieee_is_finite(x))) return
      value = is_constant(reshape(x, [size(x)]))
   end function invalid_draws

   pure function split_for_diagnostics(x) result(y)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains; the central odd iteration is discarded.
      real(dp), allocatable :: y(:, :)
      integer :: n, m, h
      n = size(x, 1)
      m = size(x, 2)
      if (n <= 1) then
         y = x
      else
         h = n / 2
         if (h == 1) then
            ! R's single-row extraction drops its dimension before cbind.
            allocate(y(m, 2))
            y(:, 1) = x(1, :)
            y(:, 2) = x(n, :)
            return
         end if
         allocate(y(h, 2 * m))
         y(:, :m) = x(:h, :)
         y(:, m + 1:) = x(n - h + 1:, :)
      end if
   end function split_for_diagnostics

   pure function r_scale(x) result(y)
      real(dp), intent(in) :: x(:) !! Numeric values, with average ranks for ties and NaNs preserved.
      real(dp) :: y(size(x))
      real(dp), allocatable :: valid(:)
      integer, allocatable :: ids(:), order(:)
      integer :: i, j, n
      y = nan_value()
      n = count(.not. ieee_is_nan(x))
      valid = pack(x, .not. ieee_is_nan(x))
      ids = pack([(i, i = 1, size(x))], .not. ieee_is_nan(x))
      call r_order(valid, order)
      i = 1
      do while (i <= n)
         j = i
         do while (j < n)
            if (valid(order(j + 1)) /= valid(order(i))) exit
            j = j + 1
         end do
         y(ids(order(i:j))) = real(i + j, dp) / 2.0_dp
         i = j + 1
      end do
   end function r_scale

   pure function u_scale(x, offset) result(y)
      real(dp), intent(in) :: x(:) !! Numeric observations; NaNs remain missing and still count in the denominator.
      real(dp), intent(in), optional :: offset !! Blom fractional offset, default 3/8.
      real(dp) :: y(size(x)), c
      c = 0.375_dp
      if (present(offset)) c = offset
      y = (r_scale(x) - c) / (real(size(x), dp) - 2.0_dp * c + 1.0_dp)
   end function u_scale

   pure function z_scale(x, offset) result(y)
      real(dp), intent(in) :: x(:) !! Numeric observations ranked jointly before inverse-normal transformation.
      real(dp), intent(in), optional :: offset !! Fractional rank offset, default 3/8.
      real(dp) :: y(size(x)), u(size(x))
      integer :: i
      u = u_scale(x, offset)
      y = nan_value()
      do i = 1, size(x)
         if (.not. ieee_is_nan(u(i))) y(i) = r_qnorm(u(i))
      end do
   end function z_scale

   pure function quantile2(x, probs, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Numeric observations, with infinities supported at exact order statistics.
      real(dp), intent(in) :: probs(:) !! Type-7 quantile probabilities in [0,1].
      logical, intent(in), optional :: na_rm !! Remove NaNs before quantiles; default false.
      real(dp) :: y(size(probs))
      real(dp), allocatable :: values(:)
      real(dp) :: h, fraction
      integer :: j, k, n
      logical :: remove
      remove = .false.
      if (present(na_rm)) remove = na_rm
      y = nan_value()
      if (.not. remove .and. any(ieee_is_nan(x))) return
      values = pack(x, .not. ieee_is_nan(x))
      n = size(values)
      if (n == 0) return
      if (all(ieee_is_finite(values))) then
         do k = 1, size(probs)
            y(k) = r_quantile_type7(values, probs(k))
         end do
         return
      end if
      call r_sort_values_in_place(values)
      do k = 1, size(probs)
         if (.not. ieee_is_finite(probs(k))) cycle
         if (probs(k) < 0.0_dp .or. probs(k) > 1.0_dp) cycle
         h = real(n - 1, dp) * probs(k)
         j = min(n, floor(h) + 1)
         fraction = h - real(j - 1, dp)
         if (fraction == 0.0_dp .or. j == n) then
            y(k) = values(j)
         else if (values(j) == values(j + 1)) then
            y(k) = values(j)
         else
            y(k) = (1.0_dp - fraction) * values(j) + fraction * values(j + 1)
         end if
      end do
   end function quantile2

   pure function autocovariance(x) result(y)
      real(dp), intent(in) :: x(:) !! Ordered observations; all lags use biased denominator N.
      real(dp) :: y(size(x))
      real(dp), allocatable :: covariance(:)
      integer :: status
      y = nan_value()
      if (size(x) < 2) return
      call r_autocovariance(x, covariance, lag_max = size(x) - 1, status = status)
      if (status == 0) y = covariance
   end function autocovariance

   pure function autocorrelation(x) result(y)
      real(dp), intent(in) :: x(:) !! Ordered observations; constant or invalid series yield NaNs.
      real(dp) :: y(size(x)), variance
      y = autocovariance(x)
      if (size(x) == 0) return
      variance = y(1)
      if (variance > 0.0_dp) then
         y = y / variance
      else
         y = nan_value()
      end if
   end function autocorrelation

   pure function basic_rhat(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Unsplit iterations by chains; at least two chains and iterations.
      real(dp) :: value, within, between, means(size(x, 2)), variances(size(x, 2))
      integer :: j, n, m
      value = nan_value()
      n = size(x, 1)
      m = size(x, 2)
      if (n < 2 .or. m < 2 .or. invalid_draws(x)) return
      do j = 1, m
         means(j) = sum(x(:, j)) / real(n, dp)
         variances(j) = sample_variance(x(:, j))
      end do
      within = sum(variances) / real(m, dp)
      between = real(n, dp) * sample_variance(means)
      if (within == 0.0_dp) then
         if (between > 0.0_dp) value = ieee_value(0.0_dp, ieee_positive_inf)
         return
      end if
      value = sqrt((between / within + real(n - 1, dp)) / real(n, dp))
   end function basic_rhat

   pure function rhat_basic(x, split) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains for one scalar variable.
      logical, intent(in), optional :: split !! Split each chain into halves; default true, dropping a central odd draw.
      real(dp) :: value
      logical :: divide
      divide = .true.
      if (present(split)) divide = split
      if (divide) then
         value = basic_rhat(split_for_diagnostics(x))
      else
         value = basic_rhat(x)
      end if
   end function rhat_basic

   pure function rhat(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains; all chains are ranked jointly.
      real(dp) :: value, bulk, tail, median(1)
      real(dp), allocatable :: split(:, :), folded(:, :)
      value = nan_value()
      if (size(x) == 0) return
      split = split_for_diagnostics(x)
      bulk = basic_rhat(reshape(z_scale(reshape(split, [size(split)])), shape(split)))
      median = quantile2(reshape(x, [size(x)]), [0.5_dp])
      folded = split_for_diagnostics(abs(x - median(1)))
      tail = basic_rhat(reshape(z_scale(reshape(folded, [size(folded)])), shape(folded)))
      if (ieee_is_nan(bulk) .or. ieee_is_nan(tail)) return
      value = max(bulk, tail)
   end function rhat

   pure function basic_ess(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Unsplit iterations by chains, requiring at least three iterations.
      real(dp) :: value, mean_var, var_plus, even, odd, tau, total
      real(dp) :: ac(size(x, 1)), rho(size(x, 1)), means(size(x, 2))
      integer :: n, m, j, t, max_t
      value = nan_value()
      n = size(x, 1)
      m = size(x, 2)
      if (n < 3 .or. invalid_draws(x)) return
      ac = 0.0_dp
      do j = 1, m
         ac = ac + autocovariance(x(:, j))
         means(j) = sum(x(:, j)) / real(n, dp)
      end do
      ac = ac / real(m, dp)
      mean_var = ac(1) * real(n, dp) / real(n - 1, dp)
      var_plus = ac(1)
      if (m > 1) var_plus = var_plus + sample_variance(means)
      if (var_plus <= 0.0_dp) return
      rho = 0.0_dp
      t = 0
      even = 1.0_dp
      rho(1) = even
      odd = 1.0_dp - (mean_var - ac(2)) / var_plus
      rho(2) = odd
      do while (t < n - 5)
         if (ieee_is_nan(even + odd)) exit
         if (even + odd <= 0.0_dp) exit
         t = t + 2
         even = 1.0_dp - (mean_var - ac(t + 1)) / var_plus
         odd = 1.0_dp - (mean_var - ac(t + 2)) / var_plus
         if (even + odd >= 0.0_dp) then
            rho(t + 1) = even
            rho(t + 2) = odd
         end if
      end do
      max_t = t
      if (even > 0.0_dp) rho(max_t + 1) = even
      t = 0
      do while (t <= max_t - 4)
         t = t + 2
         if (rho(t + 1) + rho(t + 2) > rho(t - 1) + rho(t)) then
            rho(t + 1) = (rho(t - 1) + rho(t)) / 2.0_dp
            rho(t + 2) = rho(t + 1)
         end if
      end do
      total = real(n * m, dp)
      ! R's 1:0 index includes element one; preserve the short-chain/antithetic convention.
      tau = -1.0_dp + 2.0_dp * sum(rho(1:max(1, max_t))) + rho(max_t + 1)
      tau = max(tau, 1.0_dp / log10(total))
      value = total / tau
   end function basic_ess

   pure function ess_basic(x, split) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains for one variable.
      logical, intent(in), optional :: split !! Split chains before ESS; default true.
      real(dp) :: value
      logical :: divide
      divide = .true.
      if (present(split)) divide = split
      if (divide) then
         value = basic_ess(split_for_diagnostics(x))
      else
         value = basic_ess(x)
      end if
   end function ess_basic

   pure function ess_bulk(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains; pooled rank-normalization precedes ESS.
      real(dp) :: value
      real(dp), allocatable :: split(:, :)
      value = nan_value()
      if (size(x) == 0) return
      split = split_for_diagnostics(x)
      value = basic_ess(reshape(z_scale(reshape(split, [size(split)])), shape(split)))
   end function ess_bulk

   pure function ess_mean(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains used for the mean's ESS.
      real(dp) :: value
      value = ess_basic(x)
   end function ess_mean

   pure function ess_sd(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains used for squared-centered-draw ESS.
      real(dp) :: value
      value = nan_value()
      if (size(x) == 0) return
      value = ess_mean((x - sum(x) / real(size(x), dp))**2)
   end function ess_sd

   pure function ess_quantile(x, probs) result(values)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains for a scalar variable.
      real(dp), intent(in) :: probs(:) !! Quantile probabilities in [0,1]; 1 uses the top mid-rank threshold.
      real(dp) :: values(size(probs)), p, threshold(1)
      integer :: j
      values = nan_value()
      if (invalid_draws(x)) return
      do j = 1, size(probs)
         p = probs(j)
         if (.not. ieee_is_finite(p)) cycle
         if (p < 0.0_dp .or. p > 1.0_dp) cycle
         if (p == 1.0_dp) p = 1.0_dp - 0.5_dp / real(size(x), dp)
         threshold = quantile2(reshape(x, [size(x)]), [p])
         values(j) = basic_ess(split_for_diagnostics(merge(1.0_dp, 0.0_dp, x <= threshold(1))))
      end do
   end function ess_quantile

   pure function ess_tail(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains for the minimum 5% and 95% quantile ESS.
      real(dp) :: value, tails(2)
      tails = ess_quantile(x, [0.05_dp, 0.95_dp])
      value = nan_value()
      if (any(ieee_is_nan(tails))) return
      value = minval(tails)
   end function ess_tail

   pure function ess_median(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains for median ESS.
      real(dp) :: value, median(1)
      median = ess_quantile(x, [0.5_dp])
      value = median(1)
   end function ess_median

   pure function mcse_mean(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains; MCSE is pooled sample SD divided by sqrt(ESS).
      real(dp) :: value
      value = sqrt(sample_variance(reshape(x, [size(x)])) / ess_mean(x))
   end function mcse_mean

   pure function mcse_sd(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains used in the moments-of-moments SD error estimate.
      real(dp) :: value, centered(size(x, 1), size(x, 2)), evar, varvar
      value = nan_value()
      if (invalid_draws(x)) return
      centered = x - sum(x) / real(size(x), dp)
      evar = sum(centered**2) / real(size(x), dp)
      varvar = (sum(centered**4) / real(size(x), dp) - evar**2) / ess_mean(centered**2)
      if (.not. ieee_is_finite(varvar)) return
      value = sqrt(max(0.0_dp, varvar / evar / 4.0_dp))
   end function mcse_sd

   pure function beta_quantile(p, a, b) result(value)
      real(dp), intent(in) :: p !! Interior beta cumulative probability.
      real(dp), intent(in) :: a !! Positive first beta shape parameter.
      real(dp), intent(in) :: b !! Positive second beta shape parameter.
      real(dp) :: value, lo, hi, mid
      integer :: i
      value = nan_value()
      if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) return
      if (a <= 0.0_dp .or. b <= 0.0_dp) return
      lo = 0.0_dp
      hi = 1.0_dp
      do i = 1, 100
         mid = (lo + hi) / 2.0_dp
         if (r_regularized_beta(mid, a, b) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      value = (lo + hi) / 2.0_dp
   end function beta_quantile

   pure function mcse_quantile(x, probs) result(values)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains for order-statistic uncertainty intervals.
      real(dp), intent(in) :: probs(:) !! Quantile probabilities in [0,1].
      real(dp) :: values(size(probs)), ess(size(probs)), lower, upper, p
      real(dp), allocatable :: sorted(:)
      integer :: j, n, lo, hi
      values = nan_value()
      if (invalid_draws(x)) return
      sorted = reshape(x, [size(x)])
      call r_sort_values_in_place(sorted)
      n = size(sorted)
      ess = ess_quantile(x, probs)
      do j = 1, size(probs)
         if (.not. ieee_is_finite(ess(j))) cycle
         p = probs(j)
         lower = beta_quantile(0.1586553_dp, ess(j) * p + 1.0_dp, ess(j) * (1.0_dp - p) + 1.0_dp)
         upper = beta_quantile(0.8413447_dp, ess(j) * p + 1.0_dp, ess(j) * (1.0_dp - p) + 1.0_dp)
         if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper)) cycle
         lo = max(1, floor(lower * real(n, dp)))
         hi = min(n, ceiling(upper * real(n, dp)))
         values(j) = (sorted(hi) - sorted(lo)) / 2.0_dp
      end do
   end function mcse_quantile

   pure function mcse_median(x) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains for the median's Monte Carlo standard error.
      real(dp) :: value, median(1)
      median = mcse_quantile(x, [0.5_dp])
      value = median(1)
   end function mcse_median

   pure function rhat_nested(x, superchain_ids) result(value)
      real(dp), intent(in) :: x(:, :) !! Iterations by chains, permitting one iteration per chain.
      integer, intent(in) :: superchain_ids(:) !! One group ID per chain; groups must have equal positive sizes.
      real(dp) :: value, means(size(x, 2)), variances(size(x, 2)), group_mean(size(x, 2))
      real(dp) :: within, group_within, group_between
      integer :: groups(size(x, 2)), ngroup, n, m, j, k, per_group
      logical :: selected(size(x, 2))
      value = nan_value()
      n = size(x, 1)
      m = size(x, 2)
      if (invalid_draws(x) .or. size(superchain_ids) /= m) return
      ngroup = 0
      do j = 1, m
         if (any(groups(:ngroup) == superchain_ids(j))) cycle
         ngroup = ngroup + 1
         groups(ngroup) = superchain_ids(j)
      end do
      if (ngroup < 2) return
      per_group = count(superchain_ids == groups(1))
      do j = 1, m
         means(j) = sum(x(:, j)) / real(n, dp)
         variances(j) = 0.0_dp
         if (n > 1) variances(j) = sample_variance(x(:, j))
      end do
      within = 0.0_dp
      do k = 1, ngroup
         selected = superchain_ids == groups(k)
         if (count(selected) /= per_group) return
         group_mean(k) = sum(means, mask = selected) / real(per_group, dp)
         group_within = sum(variances, mask = selected) / real(per_group, dp)
         group_between = 0.0_dp
         if (per_group > 1) group_between = sample_variance(pack(means, selected))
         within = within + group_within + group_between
      end do
      within = within / real(ngroup, dp)
      if (within == 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else
         value = sqrt(1.0_dp + sample_variance(group_mean(:ngroup)) / within)
      end if
   end function rhat_nested

   pure function summarise_draws(x) result(values)
      real(dp), intent(in) :: x(:, :, :) !! Iterations by chains by variables; returns nine default summaries per variable.
      real(dp) :: values(9, size(x, 3)), quantiles(3), center, flat(size(x, 1) * size(x, 2))
      integer :: j
      values = nan_value()
      if (size(x, 1) * size(x, 2) == 0) return
      do j = 1, size(x, 3)
         flat = reshape(x(:, :, j), [size(flat)])
         center = sum(flat) / real(size(flat), dp)
         quantiles = quantile2(flat, [0.5_dp, 0.05_dp, 0.95_dp])
         values(1, j) = center
         values(2, j) = quantiles(1)
         values(3, j) = sqrt(sample_variance(flat))
         quantiles(1:1) = quantile2(abs(flat - quantiles(1)), [0.5_dp])
         values(4, j) = 1.4826_dp * quantiles(1)
         values(5:6, j) = quantiles(2:3)
         values(7, j) = rhat(x(:, :, j))
         values(8, j) = ess_bulk(x(:, :, j))
         values(9, j) = ess_tail(x(:, :, j))
      end do
   end function summarise_draws
end module posterior_diagnostics
