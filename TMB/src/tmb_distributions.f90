module tmb_distributions
   use tmb_kinds, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: dnorm, pnorm, qnorm, dexp, pexp, qexp, dweibull, pweibull, qweibull
   public :: dbinom, dbinom_robust, dbeta, dgamma, dlnorm, dlogis
   public :: df_density, dt_density, dsn, dmultinom, dshasho, pshasho, qshasho, norm2shasho
contains
   pure elemental real(dp) function log1p_local(x) result(ans)
      real(dp), intent(in) :: x !! Value greater than -1 for which log(1+x) is required.
      real(dp) :: term, total
      integer :: k
      if (abs(x) > 1.0e-4_dp) then
         ans = log(1.0_dp + x)
         return
      end if
      total = 0.0_dp
      term = x
      do k = 1, 12
         if (mod(k, 2) == 1) then
            total = total + term / real(k, dp)
         else
            total = total - term / real(k, dp)
         end if
         term = term * x
      end do
      ans = total
   end function log1p_local

   pure elemental real(dp) function expm1_local(x) result(ans)
      real(dp), intent(in) :: x !! Value for which exp(x)-1 is required.
      real(dp) :: term, total
      integer :: k
      if (abs(x) > 1.0e-4_dp) then
         ans = exp(x) - 1.0_dp
         return
      end if
      total = 0.0_dp
      term = 1.0_dp
      do k = 1, 12
         term = term * x / real(k, dp)
         total = total + term
      end do
      ans = total
   end function expm1_local

   pure elemental real(dp) function softplus(x) result(ans)
      real(dp), intent(in) :: x !! Real argument of log(1+exp(x)), evaluated without overflow.
      if (x > 0.0_dp) then
         ans = x + log1p_local(exp(-x))
      else
         ans = log1p_local(exp(x))
      end if
   end function softplus
   pure elemental real(dp) function dnorm(x, mean, sd, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point.
      real(dp), intent(in) :: mean !! Distribution mean.
      real(dp), intent(in) :: sd !! Positive standard deviation.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: z, logans
      if (sd <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      z = (x - mean) / sd
      logans = -0.5_dp * log(2.0_dp * pi) - log(sd) - 0.5_dp * z * z
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dnorm

   pure elemental real(dp) function pnorm(q, mean, sd) result(ans)
      real(dp), intent(in) :: q !! Quantile at which the normal CDF is evaluated.
      real(dp), intent(in) :: mean !! Distribution mean.
      real(dp), intent(in) :: sd !! Positive standard deviation.
      if (sd <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else
         ans = 0.5_dp * erfc(-(q - mean) / (sd * sqrt(2.0_dp)))
      end if
   end function pnorm

   pure elemental real(dp) function qnorm(p, mean, sd) result(ans)
      real(dp), intent(in) :: p !! Probability in the closed interval [0,1].
      real(dp), intent(in) :: mean !! Distribution mean.
      real(dp), intent(in) :: sd !! Positive standard deviation.
      real(dp) :: x
      if (sd <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      if (p == 0.0_dp) then
         ans = -ieee_value(ans, ieee_positive_inf)
         return
      else if (p == 1.0_dp) then
         ans = ieee_value(ans, ieee_positive_inf)
         return
      end if
      x = inv_norm01(p)
      ans = mean + sd * x
   end function qnorm

   pure elemental real(dp) function inv_norm01(p) result(x)
      real(dp), intent(in) :: p !! Probability strictly between zero and one.
      real(dp), parameter :: a(6) = [-3.969683028665376e1_dp, 2.209460984245205e2_dp, &
                         -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp]
      real(dp), parameter :: b(5) = [-5.447609879822406e1_dp, 1.615858368580409e2_dp, &
                                     -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp]
      real(dp), parameter :: c(6) = [-7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
                               -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp]
      real(dp), parameter :: d(4) = [7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
                                     2.445134137142996_dp, 3.754408661907416_dp]
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow
      real(dp) :: q, r
      if (p < plow) then
         q = sqrt(-2.0_dp * log(p))
         x = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
             ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q * q
         x = (((((a(1) * r + a(2)) * r + a(3)) * r + a(4)) * r + a(5)) * r + a(6)) * q / &
             (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + b(5)) * r + 1.0_dp)
      else
         q = sqrt(-2.0_dp * log(1.0_dp - p))
         x = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
             ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
      end if
   end function inv_norm01

   pure elemental real(dp) function pexp(x, rate) result(ans)
      real(dp), intent(in) :: x !! Evaluation point.
      real(dp), intent(in) :: rate !! Strictly positive exponential rate.
      if (rate <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (x < 0.0_dp) then
         ans = 0.0_dp
      else
         ans = -expm1_local(-rate * x)
      end if
   end function pexp

   pure elemental real(dp) function dexp(x, rate, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point.
      real(dp), intent(in) :: rate !! Strictly positive exponential rate.
      logical, intent(in) :: give_log !! Return log density when true.
      if (rate <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (x < 0.0_dp) then
         if (give_log) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else
            ans = 0.0_dp
         end if
      else if (give_log) then
         ans = log(rate) - rate * x
      else
         ans = rate * exp(-rate * x)
      end if
   end function dexp

   pure elemental real(dp) function qexp(p, rate) result(ans)
      real(dp), intent(in) :: p !! Probability in the closed interval [0,1].
      real(dp), intent(in) :: rate !! Strictly positive exponential rate.
      if (rate <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else
         ans = -log(1.0_dp - p) / rate
      end if
   end function qexp

   pure elemental real(dp) function pweibull(x, shape, scale) result(ans)
      real(dp), intent(in) :: x !! Evaluation point.
      real(dp), intent(in) :: shape !! Strictly positive Weibull shape.
      real(dp), intent(in) :: scale !! Strictly positive Weibull scale.
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (x < 0.0_dp) then
         ans = 0.0_dp
      else
         ans = -expm1_local(-(x / scale)**shape)
      end if
   end function pweibull

   pure elemental real(dp) function dweibull(x, shape, scale, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point.
      real(dp), intent(in) :: shape !! Strictly positive Weibull shape.
      real(dp), intent(in) :: scale !! Strictly positive Weibull scale.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: logans
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (x < 0.0_dp) then
         if (give_log) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else
            ans = 0.0_dp
         end if
      else if (x == 0.0_dp .and. shape /= 1.0_dp) then
         if (shape < 1.0_dp) then
            ans = ieee_value(ans, ieee_positive_inf)
         else if (give_log) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else
            ans = 0.0_dp
         end if
      else
         logans = log(shape) - log(scale) + (shape - 1.0_dp) * log(x / scale) - (x / scale)**shape
         if (give_log) then
            ans = logans
         else
            ans = exp(logans)
         end if
      end if
   end function dweibull

   pure elemental real(dp) function qweibull(p, shape, scale) result(ans)
      real(dp), intent(in) :: p !! Probability in the closed interval [0,1].
      real(dp), intent(in) :: shape !! Strictly positive Weibull shape.
      real(dp), intent(in) :: scale !! Strictly positive Weibull scale.
      if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else
         ans = scale * (-log(1.0_dp - p))**(1.0_dp / shape)
      end if
   end function qweibull

   pure elemental real(dp) function dbinom(k, size_n, prob, give_log) result(ans)
      real(dp), intent(in) :: k !! Number of successes, normally an integer in [0,size_n].
      real(dp), intent(in) :: size_n !! Number of Bernoulli trials.
      real(dp), intent(in) :: prob !! Success probability in [0,1].
      logical, intent(in) :: give_log !! Return log probability when true.
      real(dp) :: logans
      if (k < 0.0_dp .or. k > size_n .or. prob < 0.0_dp .or. prob > 1.0_dp) then
         if (give_log) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else
            ans = 0.0_dp
         end if
         return
      end if
      logans = log_gamma(size_n + 1.0_dp) - log_gamma(k + 1.0_dp) - log_gamma(size_n - k + 1.0_dp)
      if (k > 0.0_dp) logans = logans + k * log(prob)
      if (size_n > k) logans = logans + (size_n - k) * log(1.0_dp - prob)
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dbinom

   pure elemental real(dp) function dbinom_robust(k, size_n, logit_p, give_log) result(ans)
      real(dp), intent(in) :: k !! Number of successes, normally an integer in [0,size_n].
      real(dp), intent(in) :: size_n !! Number of Bernoulli trials.
      real(dp), intent(in) :: logit_p !! Log-odds of success probability.
      logical, intent(in) :: give_log !! Return log probability when true.
      real(dp) :: logans, logp, logq
      if (k < 0.0_dp .or. k > size_n) then
         if (give_log) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else
            ans = 0.0_dp
         end if
         return
      end if
      logp = -softplus(-logit_p)
      logq = -softplus(logit_p)
      logans = log_gamma(size_n + 1.0_dp) - log_gamma(k + 1.0_dp) - log_gamma(size_n - k + 1.0_dp)
      logans = logans + k * logp + (size_n - k) * logq
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dbinom_robust

   pure elemental real(dp) function dbeta(x, shape1, shape2, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point in [0,1].
      real(dp), intent(in) :: shape1 !! Strictly positive first beta shape.
      real(dp), intent(in) :: shape2 !! Strictly positive second beta shape.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: logans
      if (shape1 <= 0.0_dp .or. shape2 <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      logans = log_gamma(shape1 + shape2) - log_gamma(shape1) - log_gamma(shape2) + &
               (shape1 - 1.0_dp) * log(x) + (shape2 - 1.0_dp) * log(1.0_dp - x)
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dbeta

   pure elemental real(dp) function dgamma(x, shape, scale, give_log) result(ans)
      real(dp), intent(in) :: x !! Nonnegative evaluation point.
      real(dp), intent(in) :: shape !! Strictly positive gamma shape.
      real(dp), intent(in) :: scale !! Strictly positive gamma scale.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: logans
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (x < 0.0_dp) then
         if (give_log) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else
            ans = 0.0_dp
         end if
      else
         logans = (shape - 1.0_dp) * log(x) - x / scale - log_gamma(shape) - shape * log(scale)
         if (give_log) then
            ans = logans
         else
            ans = exp(logans)
         end if
      end if
   end function dgamma

   pure elemental real(dp) function dlnorm(x, meanlog, sdlog, give_log) result(ans)
      real(dp), intent(in) :: x !! Positive evaluation point.
      real(dp), intent(in) :: meanlog !! Mean of log(X).
      real(dp), intent(in) :: sdlog !! Positive standard deviation of log(X).
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: logans
      if (sdlog <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (x <= 0.0_dp) then
         if (give_log) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else
            ans = 0.0_dp
         end if
      else
         logans = dnorm(log(x), meanlog, sdlog, .true.) - log(x)
         if (give_log) then
            ans = logans
         else
            ans = exp(logans)
         end if
      end if
   end function dlnorm

   pure elemental real(dp) function dlogis(x, location, scale, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point.
      real(dp), intent(in) :: location !! Logistic location parameter.
      real(dp), intent(in) :: scale !! Strictly positive logistic scale.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: z, logans
      if (scale <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      z = (x - location) / scale
      logans = -log(scale) - softplus(z) - softplus(-z)
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dlogis

   pure elemental real(dp) function df_density(x, df1, df2, give_log) result(ans)
      real(dp), intent(in) :: x !! Nonnegative evaluation point for the F density.
      real(dp), intent(in) :: df1 !! Strictly positive numerator degrees of freedom.
      real(dp), intent(in) :: df2 !! Strictly positive denominator degrees of freedom.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: logans
      if (x < 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
         if (give_log) then
            ans = -ieee_value(ans, ieee_positive_inf)
         else
            ans = 0.0_dp
         end if
         return
      end if
      logans = log_gamma(0.5_dp * (df1 + df2)) - log_gamma(0.5_dp * df1) - log_gamma(0.5_dp * df2) + &
               0.5_dp * df1 * log(df1 / df2) + (0.5_dp * df1 - 1.0_dp) * log(x) - &
               0.5_dp * (df1 + df2) * log(1.0_dp + df1 * x / df2)
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function df_density

   pure elemental real(dp) function dt_density(x, dfree, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point for the Student-t density.
      real(dp), intent(in) :: dfree !! Strictly positive degrees of freedom.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: logans
      if (dfree <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      logans = log_gamma(0.5_dp * (dfree + 1.0_dp)) - 0.5_dp * log(dfree * pi) - &
               log_gamma(0.5_dp * dfree) - 0.5_dp * (dfree + 1.0_dp) * log(1.0_dp + x * x / dfree)
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dt_density

   pure elemental real(dp) function dsn(x, alpha, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point for the standard skew-normal density.
      real(dp), intent(in) :: alpha !! Skew-normal slant parameter.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: logans
      logans = log(2.0_dp) + dnorm(x, 0.0_dp, 1.0_dp, .true.) + log(pnorm(alpha * x, 0.0_dp, 1.0_dp))
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dsn

   pure real(dp) function dmultinom(x, prob, give_log) result(ans)
      real(dp), intent(in) :: x(:) !! Nonnegative class counts.
      real(dp), intent(in) :: prob(:) !! Class probabilities, same size as x and summing to one.
      logical, intent(in) :: give_log !! Return log probability when true.
      real(dp) :: logans
      integer :: i
      if (size(x) /= size(prob) .or. any(x < 0.0_dp) .or. any(prob < 0.0_dp) .or. &
          abs(sum(prob) - 1.0_dp) > 100.0_dp * epsilon(1.0_dp)) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      logans = log_gamma(sum(x) + 1.0_dp)
      do i = 1, size(x)
         logans = logans - log_gamma(x(i) + 1.0_dp)
         if (x(i) > 0.0_dp) then
            if (prob(i) <= 0.0_dp) then
               if (give_log) then
                  ans = -ieee_value(ans, ieee_positive_inf)
               else
                  ans = 0.0_dp
               end if
               return
            end if
            logans = logans + x(i) * log(prob(i))
         end if
      end do
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dmultinom

   pure elemental real(dp) function dshasho(x, mu, sigma, nu, tau, give_log) result(ans)
      real(dp), intent(in) :: x !! Evaluation point for the sinh-asinh density.
      real(dp), intent(in) :: mu !! Location parameter.
      real(dp), intent(in) :: sigma !! Strictly positive scale parameter.
      real(dp), intent(in) :: nu !! Skewness parameter.
      real(dp), intent(in) :: tau !! Strictly positive tail-weight parameter.
      logical, intent(in) :: give_log !! Return log density when true.
      real(dp) :: z, c, r, logans
      if (sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      z = (x - mu) / sigma
      c = cosh(tau * asinh(z) - nu)
      r = sinh(tau * asinh(z) - nu)
      logans = -log(sigma) + log(tau) - 0.5_dp * log(2.0_dp * pi) - &
               0.5_dp * log(1.0_dp + z * z) + log(c) - 0.5_dp * r * r
      if (give_log) then
         ans = logans
      else
         ans = exp(logans)
      end if
   end function dshasho

   pure elemental real(dp) function pshasho(q, mu, sigma, nu, tau, give_log) result(ans)
      real(dp), intent(in) :: q !! Quantile at which to evaluate the sinh-asinh CDF.
      real(dp), intent(in) :: mu !! Location parameter.
      real(dp), intent(in) :: sigma !! Strictly positive scale parameter.
      real(dp), intent(in) :: nu !! Skewness parameter.
      real(dp), intent(in) :: tau !! Strictly positive tail-weight parameter.
      logical, intent(in) :: give_log !! Return log CDF when true.
      real(dp) :: z, r, p
      if (sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      z = (q - mu) / sigma
      r = sinh(tau * asinh(z) - nu)
      p = pnorm(r, 0.0_dp, 1.0_dp)
      if (give_log) then
         ans = log(p)
      else
         ans = p
      end if
   end function pshasho

   pure elemental real(dp) function qshasho(p, mu, sigma, nu, tau, log_p) result(ans)
      real(dp), intent(in) :: p !! Probability, or log probability when log_p is true.
      real(dp), intent(in) :: mu !! Location parameter.
      real(dp), intent(in) :: sigma !! Strictly positive scale parameter.
      real(dp), intent(in) :: nu !! Skewness parameter.
      real(dp), intent(in) :: tau !! Strictly positive tail-weight parameter.
      logical, intent(in) :: log_p !! Interpret p as a natural-log probability when true.
      real(dp) :: probability, z
      if (sigma <= 0.0_dp .or. tau <= 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      end if
      if (log_p) then
         probability = exp(p)
      else
         probability = p
      end if
      z = qnorm(probability, 0.0_dp, 1.0_dp)
      ans = mu + sigma * sinh((asinh(z) + nu) / tau)
   end function qshasho

   pure elemental real(dp) function norm2shasho(x, mu, sigma, nu, tau, log_p) result(ans)
      real(dp), intent(in) :: x !! Standard-normal input variable.
      real(dp), intent(in) :: mu !! Location parameter of the transformed variable.
      real(dp), intent(in) :: sigma !! Strictly positive scale parameter of the transformed variable.
      real(dp), intent(in) :: nu !! Skewness parameter of the transformed variable.
      real(dp), intent(in) :: tau !! Strictly positive tail-weight parameter of the transformed variable.
      logical, intent(in) :: log_p !! Passed to qshasho for upstream-compatible semantics.
      ans = qshasho(pnorm(x, 0.0_dp, 1.0_dp), mu, sigma, nu, tau, log_p)
   end function norm2shasho

end module tmb_distributions
