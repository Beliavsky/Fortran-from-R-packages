! SPDX-License-Identifier: BSD-3-Clause
! Translation of posterior 1.7.0, R/gpd.R and R/pareto_smooth.R.
module posterior_pareto
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   use r_kinds, only: dp
   use r_ordering, only: r_order
   use r_transforms, only: r_expm1, r_log1p
   use posterior_diagnostics, only: ess_tail, is_constant
   implicit none
   private
   type, public :: pareto_diagnostic
      real(dp) :: khat, min_ss, khat_threshold, convergence_rate
   end type pareto_diagnostic
   public :: gpdfit, qgeneralized_pareto, pgeneralized_pareto
   public :: ps_tail, ps_tail_length, ps_min_ss, ps_khat_threshold, ps_convergence_rate
   public :: pareto_smooth, pareto_khat, pareto_diags
   public :: pareto_min_ss, pareto_khat_threshold, pareto_convergence_rate
contains
   pure function nan_value() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_value

   pure elemental function qgeneralized_pareto(p, mu, sigma, k, lower_tail, log_p) result(value)
      real(dp), intent(in) :: p !! Cumulative probability, or log probability if log_p is true.
      real(dp), intent(in), optional :: mu !! Finite location parameter; default zero.
      real(dp), intent(in), optional :: sigma !! Positive finite scale parameter; default one.
      real(dp), intent(in), optional :: k !! Finite shape parameter; default zero, the exponential distribution.
      logical, intent(in), optional :: lower_tail !! Interpret p as lower-tail probability; default true.
      logical, intent(in), optional :: log_p !! Interpret p on the natural-log scale; default false.
      real(dp) :: value, location, scale, shape, probability, log_survival
      logical :: lower, logarithm
      location = 0.0_dp
      scale = 1.0_dp
      shape = 0.0_dp
      lower = .true.
      logarithm = .false.
      if (present(mu)) location = mu
      if (present(sigma)) scale = sigma
      if (present(k)) shape = k
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) logarithm = log_p
      value = nan_value()
      if (.not. ieee_is_finite(location) .or. .not. ieee_is_finite(scale)) return
      if (.not. ieee_is_finite(shape) .or. scale <= 0.0_dp .or. ieee_is_nan(p)) return
      probability = p
      if (logarithm) probability = exp(p)
      if (probability < 0.0_dp .or. probability > 1.0_dp) return
      if (.not. lower) probability = 1.0_dp - probability
      if (probability == 0.0_dp) then
         value = location
         return
      end if
      if (probability == 1.0_dp) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
         if (shape < 0.0_dp) value = location - scale / shape
         return
      end if
      log_survival = r_log1p(-probability)
      if (shape == 0.0_dp) then
         value = location - scale * log_survival
      else
         value = location + scale * r_expm1(-shape * log_survival) / shape
      end if
   end function qgeneralized_pareto

   pure elemental function pgeneralized_pareto(q, mu, sigma, k, lower_tail, log_p) result(value)
      real(dp), intent(in) :: q !! Quantile whose cumulative probability is requested; infinite bounds are supported.
      real(dp), intent(in), optional :: mu !! Finite location parameter; default zero.
      real(dp), intent(in), optional :: sigma !! Positive finite scale; default one.
      real(dp), intent(in), optional :: k !! Finite generalized-Pareto shape; default zero.
      logical, intent(in), optional :: lower_tail !! Return lower-tail probability; default true.
      logical, intent(in), optional :: log_p !! Return log probability; default false.
      real(dp) :: value, location, scale, shape, z
      logical :: lower, logarithm
      location = 0.0_dp
      scale = 1.0_dp
      shape = 0.0_dp
      lower = .true.
      logarithm = .false.
      if (present(mu)) location = mu
      if (present(sigma)) scale = sigma
      if (present(k)) shape = k
      if (present(lower_tail)) lower = lower_tail
      if (present(log_p)) logarithm = log_p
      value = nan_value()
      if (.not. ieee_is_finite(location) .or. .not. ieee_is_finite(scale)) return
      if (.not. ieee_is_finite(shape) .or. scale <= 0.0_dp .or. ieee_is_nan(q)) return
      z = (q - location) / scale
      if (z <= 0.0_dp) then
         value = 0.0_dp
      else if (.not. ieee_is_finite(z)) then
         value = 1.0_dp
      else if (abs(shape) < 1.0e-15_dp) then
         value = -r_expm1(-z)
      else if (1.0_dp + shape * z <= 0.0_dp) then
         value = 1.0_dp
      else
         value = -r_expm1(r_log1p(shape * z) / (-shape))
      end if
      value = min(1.0_dp, max(0.0_dp, value))
      if (.not. lower) value = 1.0_dp - value
      if (logarithm) then
         if (value == 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_negative_inf)
         else
            value = log(value)
         end if
      end if
   end function pgeneralized_pareto

   pure subroutine gpdfit(x, k, sigma, wip, min_grid_pts, sort_x, weights, status)
      real(dp), intent(in) :: x(:) !! Nonnegative excess observations, at least two, with a nonconstant lower quartile.
      real(dp), intent(out) :: k !! Estimated shape, shrunk toward 0.5 when wip is true; NaN on invalid data.
      real(dp), intent(out) :: sigma !! Estimated positive scale, before weak-prior shape shrinkage.
      logical, intent(in), optional :: wip !! Apply weakly informative Gaussian shape prior; default true.
      integer, intent(in), optional :: min_grid_pts !! Positive minimum grid count; default 30, plus floor(sqrt(N)).
      logical, intent(in), optional :: sort_x !! Sort sample and corresponding weights; default true.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative finite weights with positive total, matching observations, normalized internally to N.
      integer, intent(out), optional :: status !! Zero success, one invalid input, two unsuitable/degenerate fit.
      real(dp), allocatable :: sorted(:), w(:), theta(:), shape(:), likelihood(:), grid_weights(:)
      integer, allocatable :: order(:)
      real(dp) :: xstar, theta_hat, normalizer
      integer :: n, m, minimum, j, quartile
      logical :: prior, sorting
      k = nan_value()
      sigma = nan_value()
      if (present(status)) status = 1
      prior = .true.
      sorting = .true.
      minimum = 30
      if (present(wip)) prior = wip
      if (present(sort_x)) sorting = sort_x
      if (present(min_grid_pts)) minimum = min_grid_pts
      n = size(x)
      if (n < 2 .or. minimum < 1) return
      if (any(.not. ieee_is_finite(x)) .or. any(x < 0.0_dp)) return
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n) return
         if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
         if (maxval(weights) <= 0.0_dp) return
         w = weights / maxval(weights)
         w = w / sum(w) * real(n, dp)
      end if
      sorted = x
      if (sorting) then
         call r_order(x, order)
         sorted = x(order)
         w = w(order)
      else if (any(x(2:) < x(:n - 1))) then
         return
      end if
      quartile = max(1, floor(real(n, dp) / 4.0_dp + 0.5_dp))
      xstar = sorted(quartile)
      if (present(status)) status = 2
      if (xstar <= sorted(1)) return
      m = minimum + floor(sqrt(real(n, dp)))
      allocate(theta(m), shape(m), likelihood(m), grid_weights(m))
      do j = 1, m
         theta(j) = 1.0_dp / sorted(n) + &
            (1.0_dp - sqrt(real(m, dp) / (real(j, dp) - 0.5_dp))) / 3.0_dp / xstar
         shape(j) = sum(w * r_log1p(-theta(j) * sorted)) / real(n, dp)
         if (theta(j) == 0.0_dp .or. shape(j) == 0.0_dp) return
         likelihood(j) = real(n, dp) * (log(-theta(j) / shape(j)) - shape(j) - 1.0_dp)
      end do
      if (any(.not. ieee_is_finite(likelihood))) return
      normalizer = maxval(likelihood)
      grid_weights = exp(likelihood - normalizer)
      grid_weights = grid_weights / sum(grid_weights)
      theta_hat = sum(theta * grid_weights)
      if (theta_hat == 0.0_dp) return
      k = sum(w * r_log1p(-theta_hat * sorted)) / real(n, dp)
      sigma = -k / theta_hat
      if (prior) k = (k * real(n, dp) + 5.0_dp) / real(n + 10, dp)
      if (ieee_is_nan(k)) then
         k = ieee_value(0.0_dp, ieee_positive_inf)
         sigma = nan_value()
         return
      end if
      if (present(status)) status = 0
   end subroutine gpdfit

   pure elemental function ps_min_ss(k) result(value)
      real(dp), intent(in) :: k !! Estimated Pareto shape; k >= 1 requires infinite sample size.
      real(dp) :: value, exponent
      value = nan_value()
      if (ieee_is_nan(k)) return
      value = ieee_value(0.0_dp, ieee_positive_inf)
      if (k >= 1.0_dp) return
      exponent = 1.0_dp / (1.0_dp - max(0.0_dp, k))
      if (exponent < log10(huge(1.0_dp))) value = 10.0_dp**exponent
   end function ps_min_ss

   pure elemental function ps_khat_threshold(ndraws) result(value)
      integer, intent(in) :: ndraws !! Positive number of draws used by the diagnostic.
      real(dp) :: value
      value = nan_value()
      if (ndraws < 1) return
      if (ndraws == 1) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else
         value = 1.0_dp - 1.0_dp / log10(real(ndraws, dp))
      end if
   end function ps_khat_threshold

   pure elemental function ps_convergence_rate(k, ndraws) result(value)
      real(dp), intent(in) :: k !! Pareto shape; upstream sets rate zero at k=0 or k=1.
      integer, intent(in) :: ndraws !! Number of draws, strictly greater than one.
      real(dp) :: value, n, power
      value = nan_value()
      if (ndraws <= 1 .or. ieee_is_nan(k)) return
      value = 0.0_dp
      n = real(ndraws, dp)
      if (k < 0.0_dp) then
         value = 1.0_dp
      else if (k == 0.5_dp) then
         value = 1.0_dp - 1.0_dp / log(n)
      else if (k > 0.0_dp .and. k < 1.0_dp) then
         power = n**(2.0_dp * k)
         value = max(0.0_dp, (2.0_dp * (k - 1.0_dp) * power * n + &
            (1.0_dp - 2.0_dp * k) * power + n**2) / ((n - 1.0_dp) * (n - power)))
      end if
   end function ps_convergence_rate

   pure elemental function ps_tail_length(ndraws, r_eff) result(value)
      integer, intent(in) :: ndraws !! Nonnegative total number of draws.
      real(dp), intent(in) :: r_eff !! Positive relative effective sample size, ESS divided by draws.
      integer :: value
      real(dp) :: n
      value = -1
      if (ndraws < 0 .or. .not. ieee_is_finite(r_eff) .or. r_eff <= 0.0_dp) return
      n = real(ndraws, dp)
      if (n * r_eff > 225.0_dp) then
         value = floor(3.0_dp * sqrt(n / r_eff))
      else
         value = ndraws / 5
      end if
   end function ps_tail_length

   pure subroutine ps_tail(x, ndraws_tail, y, k, tail, smooth_draws, are_log_weights, status)
      real(dp), intent(in) :: x(:) !! Finite draws or unnormalized finite log weights.
      integer, intent(in) :: ndraws_tail !! Tail observations to fit, at least five and strictly smaller than total draws.
      real(dp), intent(out) :: y(size(x)) !! Smoothed draws in original order; log weights are shifted by max(x).
      real(dp), intent(out) :: k !! Tail generalized-Pareto shape, NaN if no suitable fit exists.
      character(*), intent(in), optional :: tail !! Tail to process, 'right' (default) or 'left'.
      logical, intent(in), optional :: smooth_draws !! Replace order statistics by fitted quantiles; default true.
      logical, intent(in), optional :: are_log_weights !! Exponentiate tail excesses after a stable log shift; default false.
      integer, intent(out), optional :: status !! Zero success, one invalid input, two unsuitable tail.
      real(dp), allocatable :: sorted(:), excess(:), fitted(:)
      integer, allocatable :: order(:)
      real(dp) :: cutoff, max_tail, sigma, probability
      integer :: n, nt, j, fit_status
      logical :: left, smooth, log_weights
      y = x
      k = nan_value()
      if (present(status)) status = 1
      n = size(x)
      nt = ndraws_tail
      if (nt < 5 .or. nt >= n .or. any(.not. ieee_is_finite(x))) return
      left = .false.
      smooth = .true.
      log_weights = .false.
      if (present(tail)) then
         if (tail /= 'left' .and. tail /= 'right') return
         left = tail == 'left'
      end if
      if (present(smooth_draws)) smooth = smooth_draws
      if (present(are_log_weights)) log_weights = are_log_weights
      if (log_weights) y = y - maxval(y)
      if (left) y = -y
      call r_order(y, order)
      sorted = y(order)
      if (present(status)) status = 2
      if (is_constant(sorted(n - nt + 1:))) then
         if (left) y = -y
         return
      end if
      cutoff = sorted(n - nt)
      if (cutoff == sorted(n - nt + 1)) cutoff = cutoff - epsilon(1.0_dp)
      max_tail = sorted(n)
      excess = sorted(n - nt + 1:)
      if (log_weights) then
         excess = exp(excess)
         cutoff = exp(cutoff)
         ! Scalar and vector EXP may differ by one ulp in optimized builds.
         ! Ordered tail values imply nonnegative excesses, including exact ties.
         excess = max(0.0_dp, excess - cutoff)
      else
         excess = excess - cutoff
      end if
      call gpdfit(excess, k, sigma, sort_x = .false., status = fit_status)
      if (fit_status == 0 .and. ieee_is_finite(k) .and. smooth) then
         allocate(fitted(nt))
         do j = 1, nt
            probability = (real(j, dp) - 0.5_dp) / real(nt, dp)
            fitted(j) = qgeneralized_pareto(probability, mu = cutoff, sigma = sigma, k = k)
         end do
         if (log_weights) fitted = log(fitted)
         y(order(n - nt + 1:)) = min(fitted, max_tail)
      end if
      if (left) y = -y
      if (present(status)) status = fit_status
   end subroutine ps_tail

   pure subroutine pareto_smooth(x, y, k, tail, r_eff, ndraws_tail, are_log_weights, smooth_draws, status)
      real(dp), intent(in) :: x(:) !! Finite scalar draws; provide r_eff to preserve external multi-chain ESS.
      real(dp), intent(out) :: y(size(x)) !! Smoothed draws, or the original draws when fitting is unavailable.
      real(dp), intent(out) :: k !! Largest estimated tail shape; NaN on invalid input.
      character(*), intent(in), optional :: tail !! 'both' (default), 'right', or 'left'; log weights force right.
      real(dp), intent(in), optional :: r_eff !! Positive relative efficiency; default tail ESS treating x as one chain.
      integer, intent(in), optional :: ndraws_tail !! Tail length; default automatic, then constrained as in posterior.
      logical, intent(in), optional :: are_log_weights !! Treat x as log weights; default false, forcing the right tail.
      logical, intent(in), optional :: smooth_draws !! Replace tails by fitted values; default true, false for diagnostics only.
      integer, intent(out), optional :: status !! Zero success, one invalid input, two unsuitable tail.
      character(5) :: side
      logical :: logarithm, smooth
      real(dp) :: efficiency, left_k, right_k, work(size(x))
      integer :: nt, n, first_status, second_status
      y = x
      k = nan_value()
      if (present(status)) status = 1
      n = size(x)
      if (n < 6 .or. any(.not. ieee_is_finite(x))) return
      if (is_constant(x)) return
      logarithm = .false.
      smooth = .true.
      side = 'both'
      if (present(are_log_weights)) logarithm = are_log_weights
      if (present(smooth_draws)) smooth = smooth_draws
      if (logarithm) then
         side = 'right'
      else if (present(tail)) then
         if (tail /= 'both' .and. tail /= 'right' .and. tail /= 'left') return
         side = tail
      end if
      efficiency = 1.0_dp
      if (present(r_eff)) then
         efficiency = r_eff
      else if (.not. present(ndraws_tail)) then
         efficiency = ess_tail(reshape(x, [n, 1])) / real(n, dp)
      end if
      nt = ps_tail_length(n, efficiency)
      if (present(ndraws_tail)) nt = ndraws_tail
      if (nt < 0) return
      nt = max(5, nt)
      if (nt >= n) return
      if (side == 'both') then
         nt = min(nt, n / 2)
         call ps_tail(x, nt, work, left_k, 'left', smooth, logarithm, first_status)
         call ps_tail(work, nt, y, right_k, 'right', smooth, logarithm, second_status)
         if (.not. ieee_is_nan(left_k) .and. .not. ieee_is_nan(right_k)) then
            k = max(left_k, right_k)
         else if (.not. ieee_is_nan(left_k)) then
            k = left_k
         else if (.not. ieee_is_nan(right_k)) then
            k = right_k
         else
            ! Preserve max(c(NA,NA), na.rm=TRUE) from the upstream both-tail path.
            k = ieee_value(0.0_dp, ieee_negative_inf)
         end if
         if (present(status)) status = min(first_status, second_status)
      else
         call ps_tail(x, nt, y, k, trim(side), smooth, logarithm, first_status)
         if (present(status)) status = first_status
      end if
   end subroutine pareto_smooth

   pure function pareto_khat(x, tail, r_eff, ndraws_tail, are_log_weights) result(k)
      real(dp), intent(in) :: x(:) !! Draws for one variable, flattened in chain-major order.
      character(*), intent(in), optional :: tail !! 'both' (default), 'right', or 'left'.
      real(dp), intent(in), optional :: r_eff !! Positive relative efficiency supplied by a chain-aware caller.
      integer, intent(in), optional :: ndraws_tail !! Explicit tail length, at least five and less than total draws.
      logical, intent(in), optional :: are_log_weights !! Interpret x as log weights and diagnose its right tail.
      real(dp) :: k, y(size(x))
      call pareto_smooth(x, y, k, tail, r_eff, ndraws_tail, are_log_weights, smooth_draws = .false.)
   end function pareto_khat

   pure function pareto_diags(x, tail, r_eff, ndraws_tail, are_log_weights) result(diags)
      real(dp), intent(in) :: x(:) !! Scalar draws to diagnose.
      character(*), intent(in), optional :: tail !! 'both' (default), 'right', or 'left'.
      real(dp), intent(in), optional :: r_eff !! Positive relative effective sample size.
      integer, intent(in), optional :: ndraws_tail !! Explicit number of tail observations.
      logical, intent(in), optional :: are_log_weights !! Whether draws are log weights, forcing the right tail.
      type(pareto_diagnostic) :: diags
      diags%khat = pareto_khat(x, tail, r_eff, ndraws_tail, are_log_weights)
      diags%min_ss = nan_value()
      diags%khat_threshold = nan_value()
      diags%convergence_rate = nan_value()
      if (ieee_is_nan(diags%khat)) return
      diags%min_ss = ps_min_ss(diags%khat)
      diags%khat_threshold = ps_khat_threshold(size(x))
      diags%convergence_rate = ps_convergence_rate(diags%khat, size(x))
   end function pareto_diags

   pure function pareto_min_ss(x) result(value)
      real(dp), intent(in) :: x(:) !! Scalar draws, using both-tail diagnostics and one-chain default ESS.
      real(dp) :: value
      value = ps_min_ss(pareto_khat(x))
   end function pareto_min_ss

   pure function pareto_khat_threshold(x) result(value)
      real(dp), intent(in) :: x(:) !! Draws whose total count determines the threshold.
      real(dp) :: value
      value = ps_khat_threshold(size(x))
   end function pareto_khat_threshold

   pure function pareto_convergence_rate(x) result(value)
      real(dp), intent(in) :: x(:) !! Scalar draws, using both-tail diagnostics and one-chain default ESS.
      real(dp) :: value
      value = ps_convergence_rate(pareto_khat(x), size(x))
   end function pareto_convergence_rate
end module posterior_pareto
