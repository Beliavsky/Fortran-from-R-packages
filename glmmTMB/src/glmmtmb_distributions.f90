! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
module glmmtmb_distributions
   use glmmtmb_kinds, only: dp
   use glmmtmb_math, only: invlogit, lambert_w0, log_bell_number, logaddexp, log1p_safe
   use tmb_distributions, only: dnorm, pnorm
   use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_quiet_nan, ieee_value
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: dbell, dbetabinom_robust, dcauchy, dcompois2, dgenpois, dnbinom_robust
   public :: dpois_glmmtmb, dskewnorm, dtweedie_compound, compois_variance, dlkj
   public :: dinvwishart_vec, dwishart_vec, mvlgamma
contains
   pure elemental real(dp) function logspace_gamma(x) result(ans)
      real(dp), intent(in) :: x !! Logarithm of a positive gamma-function argument.
      if (x < -150.0_dp) then
         ans = -x
      else
         ans = log_gamma(exp(x))
      end if
   end function logspace_gamma

   pure elemental real(dp) function dpois_glmmtmb(y, mu, give_log) result(ans)
      real(dp), intent(in) :: y !! Nonnegative count at which the Poisson mass is evaluated.
      real(dp), intent(in) :: mu !! Nonnegative Poisson mean.
      logical, intent(in) :: give_log !! Return the log probability when true.
      real(dp) :: logans
      if (y < 0.0_dp .or. mu < 0.0_dp) then
         logans = -ieee_value(logans, ieee_positive_inf)
      else if (mu == 0.0_dp) then
         if (y == 0.0_dp) then
            logans = 0.0_dp
         else
            logans = -ieee_value(logans, ieee_positive_inf)
         end if
      else
         logans = y * log(mu) - mu - log_gamma(y + 1.0_dp)
      end if
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dpois_glmmtmb

   pure elemental real(dp) function dnbinom_robust(y, log_mu, log_var_minus_mu, give_log) result(ans)
      real(dp), intent(in) :: y !! Nonnegative count at which the negative-binomial mass is evaluated.
      real(dp), intent(in) :: log_mu !! Natural logarithm of the positive distribution mean.
      real(dp), intent(in) :: log_var_minus_mu !! Natural logarithm of variance minus mean.
      logical, intent(in) :: give_log !! Return the log probability when true.
      real(dp) :: log_size, log_den, log_p, log_q, size_par, logans
      if (y < 0.0_dp) then
         logans = -ieee_value(logans, ieee_positive_inf)
      else
         log_size = 2.0_dp * log_mu - log_var_minus_mu
         if (log_size > log(huge(1.0_dp)) - 2.0_dp) then
            logans = dpois_glmmtmb(y, exp(log_mu), .true.)
         else
            size_par = exp(log_size)
            log_den = logaddexp(log_mu, log_var_minus_mu)
            log_p = log_mu - log_den
            log_q = log_var_minus_mu - log_den
            logans = log_gamma(y + size_par) - log_gamma(size_par) - log_gamma(y + 1.0_dp)
            logans = logans + size_par * log_p + y * log_q
         end if
      end if
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dnbinom_robust

   pure elemental real(dp) function dbetabinom_robust(y, loga, logb, size_n, give_log) result(ans)
      real(dp), intent(in) :: y !! Number of successes, normally an integer between zero and size_n.
      real(dp), intent(in) :: loga !! Natural logarithm of the first beta shape parameter.
      real(dp), intent(in) :: logb !! Natural logarithm of the second beta shape parameter.
      real(dp), intent(in) :: size_n !! Number of beta-binomial trials.
      logical, intent(in) :: give_log !! Return the log probability when true.
      real(dp) :: a, b, logy, lognmy, logans
      if (y < 0.0_dp .or. y > size_n) then
         logans = -ieee_value(logans, ieee_positive_inf)
      else
         a = exp(loga)
         b = exp(logb)
         if (y == 0.0_dp) then
            logy = -ieee_value(logy, ieee_positive_inf)
         else
            logy = log(y)
         end if
         if (size_n == y) then
            lognmy = -ieee_value(lognmy, ieee_positive_inf)
         else
            lognmy = log(size_n - y)
         end if
         logans = log_gamma(size_n + 1.0_dp) - log_gamma(y + 1.0_dp) - log_gamma(size_n - y + 1.0_dp)
         logans = logans + logspace_gamma(logaddexp(logy, loga))
         logans = logans + logspace_gamma(logaddexp(lognmy, logb))
         logans = logans - log_gamma(size_n + a + b) + log_gamma(a + b)
         logans = logans - logspace_gamma(loga) - logspace_gamma(logb)
      end if
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dbetabinom_robust

   pure elemental real(dp) function dgenpois(y, theta, lambda, give_log) result(ans)
      real(dp), intent(in) :: y !! Nonnegative count at which the generalized-Poisson mass is evaluated.
      real(dp), intent(in) :: theta !! Positive generalized-Poisson theta parameter.
      real(dp), intent(in) :: lambda !! Generalized-Poisson dependence parameter subject to positive support.
      logical, intent(in) :: give_log !! Return the log probability when true.
      real(dp) :: base, logans
      base = theta + lambda * y
      if (y < 0.0_dp .or. theta <= 0.0_dp .or. base <= 0.0_dp) then
         logans = -ieee_value(logans, ieee_positive_inf)
      else
         logans = log(theta) + (y - 1.0_dp) * log(base) - theta - lambda * y - log_gamma(y + 1.0_dp)
      end if
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dgenpois

   pure elemental real(dp) function dskewnorm(y, mu, sigma, alpha, give_log) result(ans)
      real(dp), intent(in) :: y !! Observation at which the skew-normal density is evaluated.
      real(dp), intent(in) :: mu !! Mean on the data scale under glmmTMB's parameterization.
      real(dp), intent(in) :: sigma !! Positive data-scale standard deviation.
      real(dp), intent(in) :: alpha !! Skew-normal slant parameter.
      logical, intent(in) :: give_log !! Return the log density when true.
      real(dp) :: delta, omega, xi, z, p, logans
      if (sigma <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      delta = alpha / sqrt(1.0_dp + alpha * alpha)
      omega = sigma / sqrt(1.0_dp - 2.0_dp * delta * delta / pi)
      xi = mu - omega * delta * sqrt(2.0_dp / pi)
      z = (y - xi) / omega
      p = pnorm(alpha * z, 0.0_dp, 1.0_dp)
      if (p <= 0.0_dp) then
         logans = -ieee_value(logans, ieee_positive_inf)
      else
         logans = log(2.0_dp) - log(omega) + dnorm(z, 0.0_dp, 1.0_dp, .true.) + log(p)
      end if
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dskewnorm

   pure elemental real(dp) function dcauchy(x, location, scale, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point for the Cauchy density.
      real(dp), intent(in) :: location !! Cauchy location parameter.
      real(dp), intent(in) :: scale !! Strictly positive Cauchy scale parameter.
      logical, intent(in) :: give_log !! Return the log density when true.
      real(dp) :: resid, logans
      if (scale <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      resid = (x - location) / scale
      logans = -log(pi) - log(scale) - log1p_safe(resid * resid)
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dcauchy

   pure elemental real(dp) function dbell(x, theta, give_log) result(ans)
      real(dp), intent(in) :: x !! Nonnegative integer-valued Bell-distribution observation.
      real(dp), intent(in) :: theta !! Strictly positive Bell parameter satisfying mean=theta*exp(theta).
      logical, intent(in) :: give_log !! Return the log probability when true.
      integer :: n
      real(dp) :: logans
      n = nint(x)
      if (x < 0.0_dp .or. theta <= 0.0_dp .or. abs(x - real(n, dp)) > 32.0_dp * epsilon(1.0_dp)) then
         logans = -ieee_value(logans, ieee_positive_inf)
      else
         logans = x * log(theta) - exp(theta) + 1.0_dp + log_bell_number(n) - log_gamma(x + 1.0_dp)
      end if
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dbell

   pure subroutine compois_stats_from_loglambda(loglambda, nu, logz, mean, variance, status)
      real(dp), intent(in) :: loglambda !! Natural logarithm of the Conway-Maxwell-Poisson lambda parameter.
      real(dp), intent(in) :: nu !! Strictly positive Conway-Maxwell-Poisson dispersion exponent.
      real(dp), intent(out) :: logz !! Log normalizing constant for the requested lambda and nu.
      real(dp), intent(out) :: mean !! Mean implied by lambda and nu.
      real(dp), intent(out) :: variance !! Variance implied by lambda and nu.
      integer, intent(out) :: status !! Zero on success, nonzero for invalid parameters or excessive required support.
      integer :: k, kmax
      real(dp) :: mode_est, lw, maxlw, sw, sk, sk2, w
      if (nu <= 0.0_dp) then
         status = 1
         logz = ieee_value(logz, ieee_quiet_nan)
         mean = logz
         variance = logz
         return
      end if
      mode_est = exp(min(loglambda / nu, log(100001.0_dp)))
      if (mode_est > 100000.0_dp) then
         status = 2
         logz = ieee_value(logz, ieee_quiet_nan)
         mean = logz
         variance = logz
         return
      end if
      kmax = max(80, ceiling(mode_est + 20.0_dp * sqrt(mode_est + 1.0_dp) + 60.0_dp))
      maxlw = 0.0_dp
      do k = 0, kmax
         lw = real(k, dp) * loglambda - nu * log_gamma(real(k + 1, dp))
         maxlw = max(maxlw, lw)
      end do
      sw = 0.0_dp
      sk = 0.0_dp
      sk2 = 0.0_dp
      do k = 0, kmax
         lw = real(k, dp) * loglambda - nu * log_gamma(real(k + 1, dp))
         w = exp(lw - maxlw)
         sw = sw + w
         sk = sk + real(k, dp) * w
         sk2 = sk2 + real(k * k, dp) * w
      end do
      logz = maxlw + log(sw)
      mean = sk / sw
      variance = max(0.0_dp, sk2 / sw - mean * mean)
      status = 0
   end subroutine compois_stats_from_loglambda

   pure subroutine compois_loglambda_from_mean(target_mean, nu, loglambda, logz, variance, status)
      real(dp), intent(in) :: target_mean !! Nonnegative desired Conway-Maxwell-Poisson mean.
      real(dp), intent(in) :: nu !! Strictly positive Conway-Maxwell-Poisson dispersion exponent.
      real(dp), intent(out) :: loglambda !! Solved logarithm of lambda matching target_mean.
      real(dp), intent(out) :: logz !! Log normalizing constant at the solved lambda.
      real(dp), intent(out) :: variance !! Variance at the solved lambda.
      integer, intent(out) :: status !! Zero on convergence, nonzero for invalid or numerically unsupported inputs.
      real(dp) :: low, high, mid, m, v, z
      integer :: iter, st
      if (target_mean < 0.0_dp .or. nu <= 0.0_dp) then
         status = 1
         loglambda = ieee_value(loglambda, ieee_quiet_nan)
         logz = loglambda
         variance = loglambda
         return
      else if (target_mean == 0.0_dp) then
         loglambda = -ieee_value(loglambda, ieee_positive_inf)
         logz = 0.0_dp
         variance = 0.0_dp
         status = 0
         return
      end if
      mid = nu * log(max(target_mean, 1.0e-12_dp))
      low = mid - 8.0_dp
      high = mid + 8.0_dp
      do iter = 1, 20
         call compois_stats_from_loglambda(low, nu, z, m, v, st)
         if (st /= 0) then
            status = st
            loglambda = ieee_value(loglambda, ieee_quiet_nan)
            logz = loglambda
            variance = loglambda
            return
         end if
         if (m <= target_mean) exit
         low = low - 8.0_dp
      end do
      do iter = 1, 20
         call compois_stats_from_loglambda(high, nu, z, m, v, st)
         if (st /= 0) then
            status = st
            loglambda = ieee_value(loglambda, ieee_quiet_nan)
            logz = loglambda
            variance = loglambda
            return
         end if
         if (m >= target_mean) exit
         high = high + 8.0_dp
      end do
      do iter = 1, 80
         mid = 0.5_dp * (low + high)
         call compois_stats_from_loglambda(mid, nu, z, m, v, st)
         if (st /= 0) exit
         if (m < target_mean) then
            low = mid
         else
            high = mid
         end if
         if (abs(m - target_mean) <= 2.0e-12_dp * max(1.0_dp, target_mean)) exit
      end do
      if (st /= 0) then
         status = st
         loglambda = ieee_value(loglambda, ieee_quiet_nan)
         logz = loglambda
         variance = loglambda
         return
      end if
      loglambda = mid
      logz = z
      variance = v
      status = 0
   end subroutine compois_loglambda_from_mean

   pure real(dp) function dcompois2(y, mean, nu, give_log) result(ans)
      real(dp), intent(in) :: y !! Nonnegative integer count for the mean-parameterized COM-Poisson mass.
      real(dp), intent(in) :: mean !! Nonnegative exact mean of the COM-Poisson distribution.
      real(dp), intent(in) :: nu !! Strictly positive COM-Poisson dispersion exponent.
      logical, intent(in) :: give_log !! Return the log probability when true.
      real(dp) :: loglambda, logz, variance, logans
      integer :: n, status
      n = nint(y)
      if (y < 0.0_dp .or. abs(y - real(n, dp)) > 32.0_dp * epsilon(1.0_dp)) then
         logans = -ieee_value(logans, ieee_positive_inf)
      else if (mean == 0.0_dp) then
         if (n == 0) then
            logans = 0.0_dp
         else
            logans = -ieee_value(logans, ieee_positive_inf)
         end if
      else
         call compois_loglambda_from_mean(mean, nu, loglambda, logz, variance, status)
         if (status /= 0) then
            logans = ieee_value(logans, ieee_quiet_nan)
         else
            logans = real(n, dp) * loglambda - nu * log_gamma(real(n + 1, dp)) - logz
         end if
      end if
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dcompois2

   pure real(dp) function compois_variance(mean, nu) result(ans)
      real(dp), intent(in) :: mean !! Nonnegative exact mean of the COM-Poisson distribution.
      real(dp), intent(in) :: nu !! Strictly positive COM-Poisson dispersion exponent.
      real(dp) :: loglambda, logz
      integer :: status
      call compois_loglambda_from_mean(mean, nu, loglambda, logz, ans, status)
      if (status /= 0) ans = ieee_value(ans, ieee_quiet_nan)
   end function compois_variance

   pure real(dp) function dtweedie_compound(y, mu, phi, power, give_log) result(ans)
      real(dp), intent(in) :: y !! Nonnegative observation for the Tweedie compound-Poisson density.
      real(dp), intent(in) :: mu !! Strictly positive Tweedie mean.
      real(dp), intent(in) :: phi !! Strictly positive Tweedie dispersion parameter.
      real(dp), intent(in) :: power !! Tweedie power strictly between one and two.
      logical, intent(in) :: give_log !! Return the log density or mass when true.
      real(dp) :: lambda, alpha, scale, logterm, logsum, jump_est
      integer :: n, nmax
      if (mu <= 0.0_dp .or. phi <= 0.0_dp .or. power <= 1.0_dp .or. power >= 2.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      else if (y < 0.0_dp) then
         logsum = -ieee_value(logsum, ieee_positive_inf)
      else
         lambda = mu**(2.0_dp - power) / (phi * (2.0_dp - power))
         if (y == 0.0_dp) then
            logsum = -lambda
         else
            alpha = (2.0_dp - power) / (power - 1.0_dp)
            scale = phi * (power - 1.0_dp) * mu**(power - 1.0_dp)
            jump_est = max(lambda, y / max(alpha * scale, tiny(1.0_dp)))
            nmax = max(200, ceiling(jump_est + 25.0_dp * sqrt(jump_est + 1.0_dp) + 80.0_dp))
            nmax = min(nmax, 100000)
            logsum = -ieee_value(logsum, ieee_positive_inf)
            do n = 1, nmax
               logterm = -lambda + real(n, dp) * log(lambda) - log_gamma(real(n + 1, dp))
               logterm = logterm + (real(n, dp) * alpha - 1.0_dp) * log(y) - y / scale
               logterm = logterm - log_gamma(real(n, dp) * alpha) - real(n, dp) * alpha * log(scale)
               logsum = logaddexp(logsum, logterm)
               if (n > ceiling(jump_est) + 20 .and. logterm < logsum - 45.0_dp) exit
            end do
         end if
      end if
      if (give_log) then
         ans = logsum
      else
         ans = exp(logsum)
      end if
   end function dtweedie_compound


   pure real(dp) function mvlgamma(x, p) result(ans)
      real(dp), intent(in) :: x !! Scalar argument of the multivariate log-gamma function.
      integer, intent(in) :: p !! Positive matrix dimension of the multivariate gamma function.
      integer :: i
      if (p < 1 .or. x <= 0.5_dp * real(p - 1, dp)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      ans = 0.25_dp * real(p * (p - 1), dp) * log(pi)
      do i = 0, p - 1
         ans = ans + log_gamma(x - 0.5_dp * real(i, dp))
      end do
   end function mvlgamma

   pure subroutine inverse_unit_lower(l, linv, status)
      real(dp), intent(in) :: l(:, :) !! Unit-diagonal lower-triangular matrix to invert.
      real(dp), intent(out) :: linv(:, :) !! Unit-diagonal lower-triangular inverse of l.
      integer, intent(out) :: status !! Zero on success, nonzero if l or linv is not square and conformable.
      real(dp) :: value
      integer :: i, j, k, n
      n = size(l, 1)
      if (size(l, 2) /= n .or. size(linv, 1) /= n .or. size(linv, 2) /= n) then
         linv = 0.0_dp
         status = 1
         return
      end if
      linv = 0.0_dp
      do i = 1, n
         linv(i, i) = 1.0_dp
         do j = 1, i - 1
            value = 0.0_dp
            do k = j, i - 1
               value = value + l(i, k) * linv(k, j)
            end do
            linv(i, j) = -value
         end do
      end do
      status = 0
   end subroutine inverse_unit_lower

   pure real(dp) function dwishart_vec(x, df, scale, give_log) result(ans)
      real(dp), intent(in) :: x(:) !! glmmTMB unconstrained Wishart vector: log diagonals then strict lower-factor entries.
      real(dp), intent(in) :: df !! Wishart degrees of freedom, required to exceed matrix dimension minus one.
      real(dp), intent(in) :: scale(:) !! Scale matrix in the same unconstrained vector parameterization as x.
      logical, intent(in) :: give_log !! Return the log density when true.
      real(dp), allocatable :: lx(:, :), ls(:, :), invls(:, :), a(:, :), xx(:, :), ssinv(:, :)
      real(dp), allocatable :: log_diag_x(:), log_diag_s(:), log_diag_d(:)
      real(dp) :: root, log_det_x, log_det_s, logres
      integer :: i, j, k, n, status
      if (size(x) /= size(scale)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      root = sqrt(1.0_dp + 8.0_dp * real(size(x), dp))
      n = nint(0.5_dp * (-1.0_dp + root))
      if (n * (n + 1) / 2 /= size(x) .or. df <= real(n - 1, dp)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      allocate(lx(n, n), ls(n, n), invls(n, n), a(n, n), xx(n, n), ssinv(n, n))
      allocate(log_diag_x(n), log_diag_s(n), log_diag_d(n))
      lx = 0.0_dp
      ls = 0.0_dp
      do i = 1, n
         lx(i, i) = 1.0_dp
         ls(i, i) = 1.0_dp
      end do
      k = n
      do i = 1, n
         do j = 1, i - 1
            k = k + 1
            lx(i, j) = x(k)
            ls(i, j) = scale(k)
         end do
      end do
      do i = 1, n
         log_diag_x(i) = log(sum(lx(i, :) * lx(i, :)))
         log_diag_s(i) = log(sum(ls(i, :) * ls(i, :)))
      end do
      log_det_x = 2.0_dp * sum(x(1:n)) - sum(log_diag_x)
      log_det_s = 2.0_dp * sum(scale(1:n)) - sum(log_diag_s)
      call inverse_unit_lower(ls, invls, status)
      if (status /= 0) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      ssinv = matmul(transpose(invls), invls)
      xx = matmul(lx, transpose(lx))
      a = ssinv * xx
      log_diag_d = x(1:n) - scale(1:n) - 0.5_dp * (log_diag_x - log_diag_s)
      do j = 1, n
         do i = 1, n
            a(i, j) = a(i, j) * exp(log_diag_d(i) + log_diag_d(j))
         end do
      end do
      logres = -0.5_dp * (df * log_det_s + (-df + real(n + 1, dp)) * log_det_x + &
         df * real(n, dp) * log(2.0_dp) + 2.0_dp * mvlgamma(0.5_dp * df, n) + sum(a))
      if (give_log) then
         ans = logres
      else
         ans = exp(logres)
      end if
   end function dwishart_vec

   pure real(dp) function dinvwishart_vec(x, df, scale, give_log) result(ans)
      real(dp), intent(in) :: x(:) !! glmmTMB unconstrained inverse-Wishart vector: log diagonals then strict lower entries.
      real(dp), intent(in) :: df !! Inverse-Wishart degrees of freedom, required to exceed matrix dimension minus one.
      real(dp), intent(in) :: scale(:) !! Scale matrix in the same unconstrained vector parameterization as x.
      logical, intent(in) :: give_log !! Return the log density when true.
      real(dp), allocatable :: lx(:, :), ls(:, :), invlx(:, :), a(:, :), ss(:, :), xinv(:, :)
      real(dp), allocatable :: log_diag_x(:), log_diag_s(:), log_diag_d(:)
      real(dp) :: root, log_det_x, log_det_s, logres
      integer :: i, j, k, n, status
      if (size(x) /= size(scale)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      root = sqrt(1.0_dp + 8.0_dp * real(size(x), dp))
      n = nint(0.5_dp * (-1.0_dp + root))
      if (n * (n + 1) / 2 /= size(x) .or. df <= real(n - 1, dp)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      allocate(lx(n, n), ls(n, n), invlx(n, n), a(n, n), ss(n, n), xinv(n, n))
      allocate(log_diag_x(n), log_diag_s(n), log_diag_d(n))
      lx = 0.0_dp
      ls = 0.0_dp
      do i = 1, n
         lx(i, i) = 1.0_dp
         ls(i, i) = 1.0_dp
      end do
      k = n
      do i = 1, n
         do j = 1, i - 1
            k = k + 1
            lx(i, j) = x(k)
            ls(i, j) = scale(k)
         end do
      end do
      do i = 1, n
         log_diag_x(i) = log(sum(lx(i, :) * lx(i, :)))
         log_diag_s(i) = log(sum(ls(i, :) * ls(i, :)))
      end do
      log_det_x = 2.0_dp * sum(x(1:n)) - sum(log_diag_x)
      log_det_s = 2.0_dp * sum(scale(1:n)) - sum(log_diag_s)
      call inverse_unit_lower(lx, invlx, status)
      if (status /= 0) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      ss = matmul(ls, transpose(ls))
      xinv = matmul(transpose(invlx), invlx)
      a = ss * xinv
      log_diag_d = scale(1:n) - x(1:n) - 0.5_dp * (log_diag_s - log_diag_x)
      do j = 1, n
         do i = 1, n
            a(i, j) = a(i, j) * exp(log_diag_d(i) + log_diag_d(j))
         end do
      end do
      logres = -0.5_dp * (-df * log_det_s + (df + real(n + 1, dp)) * log_det_x + &
         df * real(n, dp) * log(2.0_dp) + 2.0_dp * mvlgamma(0.5_dp * df, n) + sum(a))
      if (give_log) then
         ans = logres
      else
         ans = exp(logres)
      end if
   end function dinvwishart_vec
   pure real(dp) function dlkj(x, eta, give_log) result(ans)
      real(dp), intent(in) :: x(:) !! Unconstrained strict-lower-triangle correlation parameters filled row-wise.
      real(dp), intent(in) :: eta !! Positive LKJ shape used by glmmTMB's transformed-parameter prior kernel.
      logical, intent(in) :: give_log !! Return the log kernel when true.
      integer :: i, j, k, n
      real(dp), allocatable :: rownorm2(:)
      real(dp) :: root, logres
      if (size(x) == 0) then
         logres = 0.0_dp
      else
         root = sqrt(1.0_dp + 8.0_dp * real(size(x), dp))
         n = nint(0.5_dp * (1.0_dp + root))
         if (n * (n - 1) / 2 /= size(x) .or. eta <= 0.0_dp) then
            logres = ieee_value(logres, ieee_quiet_nan)
         else
            allocate(rownorm2(n))
            rownorm2 = 1.0_dp
            k = 0
            do i = 1, n
               do j = 1, i - 1
                  k = k + 1
                  rownorm2(i) = rownorm2(i) + x(k) * x(k)
               end do
            end do
            logres = -(eta - 1.0_dp) * sum(log(rownorm2))
         end if
      end if
      if (give_log) then
         ans = logres
      else
         ans = exp(logres)
      end if
   end function dlkj
end module glmmtmb_distributions
