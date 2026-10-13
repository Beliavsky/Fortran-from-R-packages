! SPDX-License-Identifier: GPL-2.0-or-later
module binom_intervals
   use binom_types
   use binom_numerics
   use r_distributions, only: r_qnorm, r_pnorm, r_dnorm
   use r_transforms, only: r_log1p, r_expm1
   use binom_profile_module, only: binom_profile
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: dp, binom_interval, binom_methods, binom_confint, binom_confint_all
   public :: binom_wilson, binom_agresti_coull, binom_exact, binom_asymp, binom_prop_test
   public :: binom_cloglog, binom_logit, binom_probit, binom_bayes, binom_lrt, binom_profile
contains
   pure elemental function binom_confint(x, n, method, conf_level) result(ci)
      real(dp), intent(in) :: x !! Success count in [0,n], including fractional counts for power calculations.
      real(dp), intent(in) :: n !! Positive trial count.
      character(len=*), intent(in), optional :: method !! Exact binom_methods name or ac/asymp alias; default wilson.
      real(dp), intent(in), optional :: conf_level !! Confidence level strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      real(dp) :: level, alpha, p, z, z2, center, se, denom, q, correction
      character(len=32) :: selected
      level = 0.95_dp
      if (present(conf_level)) level = conf_level
      ci%mean = binom_nan()
      ci%lower = binom_nan()
      ci%upper = binom_nan()
      ci%status = binom_domain
      if (.not. valid_binomial(x, n, level)) return
      p = x / n
      ci%mean = p
      ci%status = binom_ok
      selected = 'wilson'
      if (present(method)) selected = method
      alpha = 1.0_dp - level
      z = r_qnorm(0.5_dp * alpha, lower_tail=.false.)
      z2 = z * z
      select case (trim(selected))
      case ('agresti-coull', 'ac')
         center = (x + 0.5_dp * z2) / (n + z2)
         se = sqrt(center * (1.0_dp - center) / (n + z2))
         ci%lower = center - z * se
         ci%upper = center + z * se
      case ('asymptotic', 'asymp')
         se = sqrt(p * (1.0_dp - p) / n)
         ci%lower = p - z * se
         ci%upper = p + z * se
      case ('wilson')
         center = p + 0.5_dp * z2 / n
         se = z * sqrt((p * (1.0_dp - p) + 0.25_dp * z2 / n) / n)
         denom = 1.0_dp + z2 / n
         ci%lower = (center - se) / denom
         ci%upper = (center + se) / denom
         if (x == 0.0_dp) ci%lower = 0.0_dp
         if (x == n) ci%upper = 1.0_dp
      case ('prop.test')
         correction = min(0.5_dp, abs(x - 0.5_dp * n)) / n
         q = p + correction
         ci%upper = 1.0_dp
         if (q < 1.0_dp) ci%upper = (q + z2 / (2.0_dp * n) + &
            z * sqrt(q * (1.0_dp - q) / n + z2 / (4.0_dp * n * n))) / (1.0_dp + z2 / n)
         q = p - correction
         ci%lower = 0.0_dp
         if (q > 0.0_dp) ci%lower = (q + z2 / (2.0_dp * n) - &
            z * sqrt(q * (1.0_dp - q) / n + z2 / (4.0_dp * n * n))) / (1.0_dp + z2 / n)
         ci%lower = max(0.0_dp, ci%lower)
         ci%upper = min(1.0_dp, ci%upper)
      case ('exact')
         ci%lower = 0.0_dp
         ci%upper = 1.0_dp
         if (x > 0.0_dp) ci%lower = beta_quantile(0.5_dp * alpha, x, n - x + 1.0_dp)
         if (x < n) ci%upper = 1.0_dp - beta_quantile(0.5_dp * alpha, n - x, x + 1.0_dp)
      case ('logit', 'probit', 'cloglog')
         if (x == 0.0_dp .or. x == n) then
            ci%lower = 0.0_dp
            ci%upper = 1.0_dp
            if (x == 0.0_dp) ci%upper = -r_expm1(log(0.5_dp * alpha) / n)
            if (x == n) ci%lower = exp(log(0.5_dp * alpha) / n)
         else if (trim(selected) == 'logit') then
            center = logit(p)
            se = sqrt(1.0_dp / (n * p * (1.0_dp - p)))
            ci%lower = logistic(center - z * se)
            ci%upper = logistic(center + z * se)
         else if (trim(selected) == 'probit') then
            center = r_qnorm(p)
            se = sqrt(p * (1.0_dp - p) / n) / r_dnorm(center)
            ci%lower = r_pnorm(center - z * se)
            ci%upper = r_pnorm(center + z * se)
         else
            center = log(-log(p))
            se = sqrt((1.0_dp - p) / (n * p)) / abs(log(p))
            ci%lower = cloglog_inverse(center + z * se)
            ci%upper = cloglog_inverse(center - z * se)
         end if
      case ('bayes')
         ci = binom_bayes(x, n, level)
      case ('lrt')
         ci = binom_lrt(x, n, level)
      case ('profile')
         ci = binom_profile(x, n, level)
      case default
         ci%status = binom_bad_method
      end select
      if (ci%status == binom_ok) then
         if (.not. all(ieee_is_finite([ci%mean, ci%lower, ci%upper]))) ci%status = binom_convergence
      end if
   end function binom_confint

   pure function binom_confint_all(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Success count in the closed interval zero through n.
      real(dp), intent(in) :: n !! Positive trial count.
      real(dp), intent(in), optional :: conf_level !! Confidence level strictly between zero and one; default 0.95.
      type(binom_interval) :: ci(11)
      integer :: i
      do i = 1, size(binom_methods)
         ci(i) = binom_confint(x, n, binom_methods(i), conf_level)
      end do
   end function binom_confint_all

   pure elemental function binom_bayes(x, n, conf_level, interval_type, prior_shape1, prior_shape2, tol, maxit) result(ci)
      real(dp), intent(in) :: x !! Success count in the closed interval zero through n.
      real(dp), intent(in) :: n !! Positive trial count.
      real(dp), intent(in), optional :: conf_level !! Posterior probability strictly between zero and one; default 0.95.
      character(len=*), intent(in), optional :: interval_type !! highest (default) or central posterior interval.
      real(dp), intent(in), optional :: prior_shape1 !! Positive prior first beta shape; default 0.5.
      real(dp), intent(in), optional :: prior_shape2 !! Positive prior second beta shape; default 0.5.
      real(dp), intent(in), optional :: tol !! Positive relative tail-allocation tolerance; default 1e-12.
      integer, intent(in), optional :: maxit !! Positive maximum bisection iterations for highest-density intervals; default 200.
      type(binom_interval) :: ci
      real(dp) :: level, a, b, alpha, lo, hi, t, density_diff, tolerance
      integer :: iterations, i
      logical :: highest
      level = 0.95_dp
      if (present(conf_level)) level = conf_level
      a = 0.5_dp
      if (present(prior_shape1)) a = prior_shape1
      b = 0.5_dp
      if (present(prior_shape2)) b = prior_shape2
      tolerance = 1.0e-12_dp
      if (present(tol)) tolerance = tol
      iterations = 200
      if (present(maxit)) iterations = maxit
      ci%mean = binom_nan()
      ci%lower = binom_nan()
      ci%upper = binom_nan()
      ci%sig = binom_nan()
      ci%status = binom_domain
      if (.not. valid_binomial(x, n, level)) return
      if (.not. all(ieee_is_finite([a, b, tolerance]))) return
      if (a <= 0.0_dp .or. b <= 0.0_dp .or. tolerance <= 0.0_dp .or. iterations < 1) return
      highest = .true.
      if (present(interval_type)) then
         if (interval_type /= 'highest' .and. interval_type /= 'central') then
            ci%status = binom_bad_method
            return
         end if
         highest = interval_type == 'highest'
      end if
      a = a + x
      b = b + n - x
      ci%shape1 = a
      ci%shape2 = b
      ci%mean = a / (a + b)
      alpha = 1.0_dp - level
      ci%status = binom_ok
      if (x == 0.0_dp) then
         ci%lower = 0.0_dp
         ci%upper = 1.0_dp - beta_quantile(alpha, b, a)
      else if (x == n) then
         ci%lower = beta_quantile(alpha, a, b)
         ci%upper = 1.0_dp
      else
         t = 0.5_dp * alpha
         if (highest) then
            ! R's interior highest-density algorithm assumes a unimodal posterior.
            if (a <= 1.0_dp .or. b <= 1.0_dp) then
               ci%status = binom_domain
               return
            end if
            lo = 0.0_dp
            hi = alpha
            do i = 1, iterations
               t = lo + 0.5_dp * (hi - lo)
               if (t == lo .or. t == hi) then
                  ci%status = binom_convergence
                  exit
               end if
               ci%lower = beta_quantile(t, a, b)
               ci%upper = 1.0_dp - beta_quantile(alpha - t, b, a)
               density_diff = (a - 1.0_dp) * (log(ci%lower) - log(ci%upper)) + &
                  (b - 1.0_dp) * (r_log1p(-ci%lower) - r_log1p(-ci%upper))
               if (density_diff < 0.0_dp) then
                  lo = t
               else
                  hi = t
               end if
               if (hi - lo <= tolerance * alpha) exit
            end do
            if (i > iterations) ci%status = binom_convergence
         end if
         ci%lower = beta_quantile(t, a, b)
         ci%upper = 1.0_dp - beta_quantile(alpha - t, b, a)
      end if
      ci%sig = beta_probability(ci%lower, a, b) + beta_probability(ci%upper, a, b, .true.)
      if (.not. all(ieee_is_finite([ci%mean, ci%lower, ci%upper, ci%sig]))) ci%status = binom_convergence
   end function binom_bayes

   pure elemental function binom_lrt(x, n, conf_level, bayes, conf_adj, prior_shape1, prior_shape2, tol, maxit) result(ci)
      real(dp), intent(in) :: x !! Success count in the closed interval zero through n.
      real(dp), intent(in) :: n !! Positive trial count.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      logical, intent(in), optional :: bayes !! Adjust boundary observations with beta prior pseudocounts; default false.
      logical, intent(in), optional :: conf_adj !! Adjust tails at boundary or tolerance-clipped observations; default false.
      real(dp), intent(in), optional :: prior_shape1 !! Positive first prior shape used when bayes is true; default 0.5.
      real(dp), intent(in), optional :: prior_shape2 !! Positive second prior shape used when bayes is true; default 0.5.
      real(dp), intent(in), optional :: tol !! Positive bracket/root tolerance below 0.5; omitted uses R's two separate defaults.
      integer, intent(in), optional :: maxit !! Brent iterations per endpoint, 1 through 1000000; default 1000.
      type(binom_interval) :: ci
      real(dp) :: level, alpha, xx, nn, a, b, z, low, high, p, bracket_tol, root_tol
      logical :: adjust, adjust_level
      integer :: iterations, info
      level = 0.95_dp
      if (present(conf_level)) level = conf_level
      adjust = .false.
      if (present(bayes)) adjust = bayes
      adjust_level = .false.
      if (present(conf_adj)) adjust_level = conf_adj
      a = 0.5_dp
      if (present(prior_shape1)) a = prior_shape1
      b = 0.5_dp
      if (present(prior_shape2)) b = prior_shape2
      bracket_tol = sqrt(epsilon(1.0_dp))
      root_tol = sqrt(bracket_tol)
      if (present(tol)) then
         bracket_tol = tol
         root_tol = tol
      end if
      iterations = 1000
      if (present(maxit)) iterations = maxit
      ci%mean = binom_nan()
      ci%lower = binom_nan()
      ci%upper = binom_nan()
      ci%status = binom_domain
      if (.not. valid_binomial(x, n, level)) return
      if (.not. all(ieee_is_finite([a, b, bracket_tol]))) return
      if (a <= 0.0_dp .or. b <= 0.0_dp .or. bracket_tol <= 0.0_dp .or. bracket_tol >= 0.5_dp) return
      if (iterations < 1) return
      if (iterations > 1000000) then
         ci%status = binom_resource
         return
      end if
      ci%mean = x / n
      ci%status = binom_ok
      xx = x
      nn = n
      alpha = 1.0_dp - level
      z = r_qnorm(0.5_dp * alpha, lower_tail=.false.)
      if ((x == 0.0_dp .or. x == n) .and. adjust_level) z = r_qnorm(0.25_dp * alpha, lower_tail=.false.)
      if ((x == 0.0_dp .or. x == n) .and. adjust) then
         ! Reject unrepresentable posterior counts before performing overflowing additions.
         if (n > huge(1.0_dp) - a) then
            ci%status = binom_resource
            return
         end if
         nn = n + a
         if (nn > huge(1.0_dp) - b) then
            ci%status = binom_resource
            return
         end if
         nn = nn + b
         xx = x + a
      end if
      p = xx / nn
      ci%lower = 0.0_dp
      ci%upper = 1.0_dp
      if (x > 0.0_dp .and. bracket_tol < p) then
         high = p
         if (p == 1.0_dp) high = 1.0_dp - bracket_tol
         if (lrt_deviate(bracket_tol, xx, nn, bracket_tol) + z > 0.0_dp) then
            if (adjust_level) z = r_qnorm(alpha, lower_tail=.false.)
         else
            call lrt_brent(bracket_tol, high, xx, nn, -z, bracket_tol, root_tol, iterations, ci%lower, info)
            if (info /= 0) ci%status = binom_convergence
         end if
      end if
      if (x < n .and. p < 1.0_dp - bracket_tol) then
         low = p
         if (lrt_deviate(1.0_dp - bracket_tol, xx, nn, bracket_tol) - z < 0.0_dp) then
            if (adjust_level) z = r_qnorm(alpha, lower_tail=.false.)
            if (lrt_deviate(bracket_tol, xx, nn, bracket_tol) + z <= 0.0_dp) then
               call lrt_brent(bracket_tol, p, xx, nn, -z, bracket_tol, root_tol, iterations, ci%lower, info)
               if (info /= 0) ci%status = binom_convergence
            end if
         else
            call lrt_brent(low, 1.0_dp - bracket_tol, xx, nn, z, bracket_tol, root_tol, iterations, ci%upper, info)
            if (info /= 0) ci%status = binom_convergence
         end if
      end if
   end function binom_lrt

   pure elemental real(dp) function lrt_deviate(candidate, x, n, tol) result(value)
      real(dp), intent(in) :: candidate !! Candidate probability in [0,1]; endpoints occur only at the MLE.
      real(dp), intent(in) :: x !! Adjusted success count in [0,n].
      real(dp), intent(in) :: n !! Positive adjusted trial count.
      real(dp), intent(in) :: tol !! Width around the MLE within which upstream sets signed deviance to zero.
      real(dp) :: difference
      difference = candidate - x / n
      value = 0.0_dp
      if (abs(difference) < tol) return
      value = sign(sqrt(binomial_deviance(candidate, x, n)), difference)
   end function lrt_deviate

   pure subroutine lrt_brent(ax, bx, x, n, bound, bracket_tol, tol, maxit, root, status)
      real(dp), intent(in) :: ax !! Lower probability bracket endpoint.
      real(dp), intent(in) :: bx !! Upper probability bracket endpoint.
      real(dp), intent(in) :: x !! Adjusted number of successes.
      real(dp), intent(in) :: n !! Positive adjusted trial count.
      real(dp), intent(in) :: bound !! Signed normal deviate targeted by the likelihood root.
      real(dp), intent(in) :: bracket_tol !! MLE-neighborhood zero threshold used by the upstream deviance function.
      real(dp), intent(in) :: tol !! Positive absolute Brent root tolerance.
      integer, intent(in) :: maxit !! Maximum number of Brent iterations.
      real(dp), intent(out) :: root !! Best root approximation; NaN if the initial bracket is invalid.
      integer, intent(out) :: status !! Zero on convergence, two on invalid bracket or iteration exhaustion.
      real(dp) :: a, b, c, fa, fb, fc, previous, active_tol, step, p, q, cb, t1, t2
      integer :: i
      a = ax
      b = bx
      fa = lrt_deviate(a, x, n, bracket_tol) - bound
      fb = lrt_deviate(b, x, n, bracket_tol) - bound
      root = binom_nan()
      status = binom_convergence
      if (.not. all(ieee_is_finite([fa, fb])) .or. ax > bx) return
      if ((fa > 0.0_dp .and. fb > 0.0_dp) .or. (fa < 0.0_dp .and. fb < 0.0_dp)) return
      c = a
      fc = fa
      ! Brent's zeroin iteration, translated from upstream binom/src/binom.c.
      do i = 1, maxit + 1
         previous = b - a
         if (abs(fc) < abs(fb)) then
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
         end if
         active_tol = 2.0_dp * epsilon(1.0_dp) * abs(b) + 0.5_dp * tol
         step = 0.5_dp * (c - b)
         if (abs(step) <= active_tol .or. fb == 0.0_dp) then
            root = b
            status = binom_ok
            return
         end if
         if (abs(previous) >= active_tol .and. abs(fa) > abs(fb)) then
            cb = c - b
            if (a == c) then
               t1 = fb / fa
               p = cb * t1
               q = 1.0_dp - t1
            else
               q = fa / fc
               t1 = fb / fc
               t2 = fb / fa
               p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0_dp))
               q = (q - 1.0_dp) * (t1 - 1.0_dp) * (t2 - 1.0_dp)
            end if
            if (p > 0.0_dp) then
               q = -q
            else
               p = -p
            end if
            if (p < 0.75_dp * cb * q - 0.5_dp * abs(active_tol * q) .and. &
               p < 0.5_dp * abs(previous * q)) step = p / q
         end if
         if (abs(step) < active_tol) then
            if (step > 0.0_dp) then
               step = active_tol
            else
               step = -active_tol
            end if
         end if
         a = b
         fa = fb
         b = b + step
         fb = lrt_deviate(b, x, n, bracket_tol) - bound
         if (.not. ieee_is_finite(fb)) return
         if ((fb > 0.0_dp .and. fc > 0.0_dp) .or. (fb < 0.0_dp .and. fc < 0.0_dp)) then
            c = a
            fc = fa
         end if
      end do
      root = b
   end subroutine lrt_brent

   pure elemental function binom_wilson(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Number of successes between zero and n.
      real(dp), intent(in) :: n !! Positive number of trials.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      ci = binom_confint(x, n, 'wilson', conf_level)
   end function binom_wilson

   pure elemental function binom_agresti_coull(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Number of successes between zero and n.
      real(dp), intent(in) :: n !! Positive number of trials.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      ci = binom_confint(x, n, 'agresti-coull', conf_level)
   end function binom_agresti_coull

   pure elemental function binom_exact(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Number of successes between zero and n.
      real(dp), intent(in) :: n !! Positive number of trials.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      ci = binom_confint(x, n, 'exact', conf_level)
   end function binom_exact

   pure elemental function binom_asymp(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Number of successes between zero and n.
      real(dp), intent(in) :: n !! Positive number of trials.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      ci = binom_confint(x, n, 'asymptotic', conf_level)
   end function binom_asymp

   pure elemental function binom_prop_test(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Number of successes between zero and n.
      real(dp), intent(in) :: n !! Positive number of trials.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      ci = binom_confint(x, n, 'prop.test', conf_level)
   end function binom_prop_test

   pure elemental function binom_cloglog(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Number of successes between zero and n.
      real(dp), intent(in) :: n !! Positive number of trials.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      ci = binom_confint(x, n, 'cloglog', conf_level)
   end function binom_cloglog

   pure elemental function binom_logit(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Number of successes between zero and n.
      real(dp), intent(in) :: n !! Positive number of trials.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      ci = binom_confint(x, n, 'logit', conf_level)
   end function binom_logit

   pure elemental function binom_probit(x, n, conf_level) result(ci)
      real(dp), intent(in) :: x !! Number of successes between zero and n.
      real(dp), intent(in) :: n !! Positive number of trials.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      type(binom_interval) :: ci
      ci = binom_confint(x, n, 'probit', conf_level)
   end function binom_probit
   pure elemental real(dp) function cloglog_inverse(eta) result(value)
      real(dp), intent(in) :: eta !! Complementary-log-log coordinate, including extreme finite tails.
      if (eta > log(huge(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(-exp(eta))
      end if
   end function cloglog_inverse
end module binom_intervals
