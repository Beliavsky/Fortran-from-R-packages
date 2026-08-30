module learnbayes_models
   use learnbayes_distributions, only: dmnorm, gamma_logpdf, normal_logpdf, student_t_logpdf
   use learnbayes_kinds, only: dp
   use learnbayes_linalg, only: inverse_matrix
   use learnbayes_special, only: beta_log, logistic, normal_cdf
   implicit none
   private

   public :: betabinexch
   public :: betabinexch0
   public :: bfexch
   public :: bradley_terry_post
   public :: cauchyerrorpost
   public :: groupeddatapost
   public :: howardprior
   public :: lbinorm
   public :: logctablepost
   public :: logisticpost
   public :: logpoissgamma
   public :: logpoissnormal
   public :: mnormt_onesided
   public :: mnormt_twosided
   public :: normchi2post
   public :: normnormexch
   public :: poissgamexch
   public :: reg_gprior_post
   public :: transplantpost
   public :: weibullregpost

contains

   function betabinexch(theta, data) result(value)
      real(dp), intent(in) :: theta(2) !! Transformed parameters: logit population mean and log precision K.
      real(dp), intent(in) :: data(:, :) !! Binomial-group matrix with successes in column 1 and trials in column 2.
      real(dp) :: value
      real(dp) :: eta
      real(dp) :: kappa
      real(dp) :: y
      real(dp) :: n
      integer :: i

      eta = logistic(theta(1))
      kappa = exp(theta(2))
      value = 0.0_dp
      do i = 1, size(data, 1)
         y = data(i, 1)
         n = data(i, 2)
         value = value + beta_log(kappa*eta + y, kappa*(1.0_dp - eta) + n - y) - &
            beta_log(kappa*eta, kappa*(1.0_dp - eta))
      end do
      value = value + theta(2) - 2.0_dp*log(1.0_dp + exp(theta(2)))
   end function betabinexch

   function betabinexch0(theta, data) result(value)
      real(dp), intent(in) :: theta(2) !! Natural parameters: population mean eta in (0,1) and positive precision K.
      real(dp), intent(in) :: data(:, :) !! Binomial-group matrix with successes in column 1 and trials in column 2.
      real(dp) :: value
      real(dp) :: eta
      real(dp) :: kappa
      real(dp) :: y
      real(dp) :: n
      integer :: i

      eta = theta(1)
      kappa = theta(2)
      if (eta <= 0.0_dp .or. eta >= 1.0_dp .or. kappa <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do i = 1, size(data, 1)
         y = data(i, 1)
         n = data(i, 2)
         value = value + beta_log(kappa*eta + y, kappa*(1.0_dp - eta) + n - y) - &
            beta_log(kappa*eta, kappa*(1.0_dp - eta))
      end do
      value = value - 2.0_dp*log(1.0_dp + kappa) - log(eta) - log(1.0_dp - eta)
   end function betabinexch0

   function bfexch(theta, data, kappa) result(value)
      real(dp), intent(in) :: theta !! Logit population mean used in the exchangeable beta-binomial model.
      real(dp), intent(in) :: data(:, :) !! Binomial-group matrix with successes in column 1 and trials in column 2.
      real(dp), intent(in) :: kappa !! Positive fixed beta-binomial precision parameter K.
      real(dp) :: value
      real(dp) :: eta
      real(dp) :: sy
      real(dp) :: sf
      integer :: i

      eta = logistic(theta)
      value = 0.0_dp
      sy = 0.0_dp
      sf = 0.0_dp
      do i = 1, size(data, 1)
         value = value + beta_log(kappa*eta + data(i, 1), &
            kappa*(1.0_dp - eta) + data(i, 2) - data(i, 1)) - beta_log(kappa*eta, kappa*(1.0_dp - eta))
         sy = sy + data(i, 1)
         sf = sf + data(i, 2) - data(i, 1)
      end do
      value = value + log(eta*(1.0_dp - eta)) - beta_log(sy + 1.0_dp, sf + 1.0_dp)
   end function bfexch

   function bradley_terry_post(theta, data) result(value)
      real(dp), intent(in) :: theta(:) !! Player log-abilities followed by log prior standard deviation as the last element.
      real(dp), intent(in) :: data(:, :) !! Match matrix: player i, player j, i wins, then j wins, using one-based IDs.
      real(dp) :: value
      real(dp) :: p
      real(dp) :: sigma
      integer :: i
      integer :: j
      integer :: k
      integer :: m

      m = size(theta)
      sigma = exp(theta(m))
      value = 0.0_dp
      do k = 1, size(data, 1)
         i = nint(data(k, 1))
         j = nint(data(k, 2))
         p = logistic(theta(i) - theta(j))
         value = value + data(k, 3)*log(max(p, tiny(1.0_dp))) + &
            data(k, 4)*log(max(1.0_dp - p, tiny(1.0_dp)))
      end do
      do i = 1, m - 1
         value = value + normal_logpdf(theta(i), 0.0_dp, sigma)
      end do
   end function bradley_terry_post

   function cauchyerrorpost(theta, data) result(value)
      real(dp), intent(in) :: theta(2) !! Cauchy location and log scale parameters.
      real(dp), intent(in) :: data(:) !! Independent observations modeled with Cauchy errors.
      real(dp) :: value
      real(dp) :: scale
      integer :: i

      scale = exp(theta(2))
      value = 0.0_dp
      do i = 1, size(data)
         value = value + student_t_logpdf((data(i) - theta(1))/scale, 1.0_dp) - log(scale)
      end do
   end function cauchyerrorpost

   function groupeddatapost(theta, freq, int_lo, int_hi) result(value)
      real(dp), intent(in) :: theta(2) !! Normal mean and log standard deviation.
      real(dp), intent(in) :: freq(:) !! Nonnegative frequencies for grouped-data intervals.
      real(dp), intent(in) :: int_lo(:) !! Lower endpoints corresponding one-for-one with freq.
      real(dp), intent(in) :: int_hi(:) !! Upper endpoints corresponding one-for-one with freq.
      real(dp) :: value
      real(dp) :: mu
      real(dp) :: sigma
      real(dp) :: p
      integer :: i

      mu = theta(1)
      sigma = exp(theta(2))
      value = 0.0_dp
      do i = 1, size(freq)
         p = normal_cdf(int_hi(i), mu, sigma) - normal_cdf(int_lo(i), mu, sigma)
         value = value + freq(i)*log(max(p, tiny(1.0_dp)))
      end do
   end function groupeddatapost

   pure function howardprior(xy, par) result(value)
      real(dp), intent(in) :: xy(2) !! Two probabilities p1 and p2 in the open unit interval.
      real(dp), intent(in) :: par(5) !! Howard-prior parameters alpha, beta, gamma, delta, and positive sigma.
      real(dp) :: value
      real(dp) :: p1
      real(dp) :: p2
      real(dp) :: u

      p1 = xy(1)
      p2 = xy(2)
      if (p1 <= 0.0_dp .or. p1 >= 1.0_dp .or. p2 <= 0.0_dp .or. p2 >= 1.0_dp .or. par(5) <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      u = log(p1/(1.0_dp - p1)*(1.0_dp - p2)/p2)/par(5)
      value = -0.5_dp*u*u + (par(1) - 1.0_dp)*log(p1) + (par(2) - 1.0_dp)*log(1.0_dp - p1) + &
         (par(3) - 1.0_dp)*log(p2) + (par(4) - 1.0_dp)*log(1.0_dp - p2)
   end function howardprior

   function lbinorm(xy, mean, varcov) result(value)
      real(dp), intent(in) :: xy(2) !! Bivariate point at which the unnormalized Gaussian log density is evaluated.
      real(dp), intent(in) :: mean(2) !! Bivariate Gaussian mean.
      real(dp), intent(in) :: varcov(2, 2) !! Positive-definite bivariate covariance matrix.
      real(dp) :: value
      real(dp) :: r
      real(dp) :: zx
      real(dp) :: zy

      zx = (xy(1) - mean(1))/sqrt(varcov(1, 1))
      zy = (xy(2) - mean(2))/sqrt(varcov(2, 2))
      r = varcov(1, 2)/sqrt(varcov(1, 1)*varcov(2, 2))
      value = -0.5_dp/(1.0_dp - r*r)*(zx*zx - 2.0_dp*r*zx*zy + zy*zy)
   end function lbinorm

   pure function logctablepost(theta, data) result(value)
      real(dp), intent(in) :: theta(2) !! Log-odds contrast and overall log-odds parameters.
      real(dp), intent(in) :: data(4) !! Success/failure counts s1, f1, s2, f2 for two binomial groups.
      real(dp) :: value
      real(dp) :: logitp1
      real(dp) :: logitp2

      logitp1 = 0.5_dp*(theta(1) + theta(2))
      logitp2 = 0.5_dp*(theta(2) - theta(1))
      value = data(1)*logitp1 - (data(1) + data(2))*log(1.0_dp + exp(logitp1)) + &
         data(3)*logitp2 - (data(3) + data(4))*log(1.0_dp + exp(logitp2))
   end function logctablepost

   function logisticpost(beta, data) result(value)
      real(dp), intent(in) :: beta(2) !! Logistic-regression intercept and slope.
      real(dp), intent(in) :: data(:, :) !! Grouped logistic data with x, trials n, and successes y in columns 1:3.
      real(dp) :: value
      real(dp) :: p
      integer :: i

      value = 0.0_dp
      do i = 1, size(data, 1)
         p = logistic(beta(1) + beta(2)*data(i, 1))
         value = value + data(i, 3)*log(max(p, tiny(1.0_dp))) + &
            (data(i, 2) - data(i, 3))*log(max(1.0_dp - p, tiny(1.0_dp)))
      end do
   end function logisticpost

   function logpoissgamma(theta, y, prior_par) result(value)
      real(dp), intent(in) :: theta !! Log Poisson rate lambda.
      real(dp), intent(in) :: y(:) !! Observed Poisson counts, represented numerically as in the R package.
      real(dp), intent(in) :: prior_par(2) !! Gamma prior shape and rate for lambda.
      real(dp) :: value
      real(dp) :: lambda

      lambda = exp(theta)
      value = gamma_logpdf(lambda, sum(y) + 1.0_dp, real(size(y), dp)) + &
         gamma_logpdf(lambda, prior_par(1), prior_par(2)) + theta
   end function logpoissgamma

   function logpoissnormal(theta, y, prior_par) result(value)
      real(dp), intent(in) :: theta !! Log Poisson rate lambda.
      real(dp), intent(in) :: y(:) !! Observed Poisson counts, represented numerically as in the R package.
      real(dp), intent(in) :: prior_par(2) !! Normal prior mean and standard deviation for log(lambda).
      real(dp) :: value
      real(dp) :: lambda

      lambda = exp(theta)
      value = gamma_logpdf(lambda, sum(y) + 1.0_dp, real(size(y), dp)) + &
         normal_logpdf(theta, prior_par(1), prior_par(2))
   end function logpoissnormal

   pure function normchi2post(theta, data) result(value)
      real(dp), intent(in) :: theta(2) !! Normal mean and positive variance sigma squared.
      real(dp), intent(in) :: data(:) !! Independent Gaussian observations.
      real(dp) :: value
      integer :: i

      if (theta(2) <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      value = -log(theta(2))
      do i = 1, size(data)
         value = value - 0.5_dp*(data(i) - theta(1))**2/theta(2) - 0.5_dp*log(theta(2))
      end do
   end function normchi2post

   function normnormexch(theta, data) result(value)
      real(dp), intent(in) :: theta(2) !! Population mean and log between-group standard deviation tau.
      real(dp), intent(in) :: data(:, :) !! Group means and known variances in columns 1 and 2.
      real(dp) :: value
      real(dp) :: tau
      integer :: i

      tau = exp(theta(2))
      value = theta(2)
      do i = 1, size(data, 1)
         value = value + normal_logpdf(data(i, 1), theta(1), sqrt(data(i, 2) + tau*tau))
      end do
   end function normnormexch

   function poissgamexch(theta, data, z0) result(value)
      real(dp), intent(in) :: theta(2) !! Log gamma shape alpha and log mean mu for exchangeable Poisson rates.
      real(dp), intent(in) :: data(:, :) !! Exposure in column 1 and Poisson count in column 2.
      real(dp), intent(in) :: z0 !! Positive scale hyperparameter in the prior on alpha.
      real(dp) :: value
      real(dp) :: alpha
      real(dp) :: beta
      real(dp) :: mu
      real(dp) :: e
      real(dp) :: y
      integer :: i

      alpha = exp(theta(1))
      mu = exp(theta(2))
      beta = alpha/mu
      value = 0.0_dp
      do i = 1, size(data, 1)
         e = data(i, 1)
         y = data(i, 2)
         value = value + log_gamma(alpha + y) - (y + alpha)*log(e + beta) + &
            alpha*log(beta) - log_gamma(alpha)
      end do
      value = value + log(alpha) - 2.0_dp*log(alpha + z0)
   end function poissgamexch

   function reg_gprior_post(theta, y, x, c0, beta0) result(value)
      real(dp), intent(in) :: theta(:) !! Regression coefficients followed by log residual standard deviation.
      real(dp), intent(in) :: y(:) !! Gaussian regression response vector.
      real(dp), intent(in) :: x(:, :) !! Regression design matrix with one row per observation.
      real(dp), intent(in) :: c0 !! Positive g-prior multiplier applied to sigma squared times inverse(X'X).
      real(dp), intent(in) :: beta0(:) !! Prior mean vector for regression coefficients.
      real(dp) :: value
      real(dp), allocatable :: xtx(:, :)
      real(dp), allocatable :: xtx_inv(:, :)
      real(dp), allocatable :: cov(:, :)
      real(dp), allocatable :: beta(:)
      real(dp) :: sigma
      integer :: i
      integer :: info
      integer :: p

      p = size(x, 2)
      allocate(beta(p), xtx(p, p), xtx_inv(p, p), cov(p, p))
      beta = theta(1:p)
      sigma = exp(theta(p + 1))
      value = 0.0_dp
      do i = 1, size(y)
         value = value + normal_logpdf(y(i), dot_product(x(i, :), beta), sigma)
      end do
      xtx = matmul(transpose(x), x)
      call inverse_matrix(xtx, xtx_inv, info)
      if (info /= 0) then
         value = -huge(1.0_dp)
         return
      end if
      cov = c0*sigma*sigma*xtx_inv
      value = value + dmnorm(beta, beta0, cov, .true.)
   end function reg_gprior_post

   function transplantpost(theta, data) result(value)
      real(dp), intent(in) :: theta(3) !! Log treatment acceleration tau, log lambda, and log shape p.
      real(dp), intent(in) :: data(:, :) !! Matrix with x, y, treatment indicator t, and censoring indicator d.
      real(dp) :: value
      real(dp) :: tau
      real(dp) :: lambda
      real(dp) :: p
      real(dp) :: x
      real(dp) :: y
      real(dp) :: denom
      integer :: i
      integer :: tr
      integer :: d

      tau = exp(theta(1))
      lambda = exp(theta(2))
      p = exp(theta(3))
      value = 0.0_dp
      do i = 1, size(data, 1)
         x = data(i, 1)
         y = data(i, 2)
         tr = nint(data(i, 3))
         d = nint(data(i, 4))
         if (tr == 0) then
            denom = lambda + x
            if (d == 0) then
               value = value + p*log(lambda) + log(p) - (p + 1.0_dp)*log(denom)
            else
               value = value + p*log(lambda/denom)
            end if
         else
            denom = lambda + y + tau*x
            if (d == 0) then
               value = value + p*log(lambda) + log(p*tau) - (p + 1.0_dp)*log(denom)
            else
               value = value + p*log(lambda/denom)
            end if
         end if
      end do
      value = value + sum(theta)
   end function transplantpost

   function weibullregpost(theta, data) result(value)
      real(dp), intent(in) :: theta(:) !! Log Weibull scale sigma, intercept mu, then regression coefficients.
      real(dp), intent(in) :: data(:, :) !! Survival matrix: time, event indicator, then covariates.
      real(dp) :: value
      real(dp) :: sigma
      real(dp) :: mu
      real(dp) :: z
      real(dp) :: logf
      real(dp) :: logs
      integer :: i
      integer :: p

      p = size(data, 2) - 2
      sigma = exp(theta(1))
      mu = theta(2)
      value = 0.0_dp
      do i = 1, size(data, 1)
         z = (log(data(i, 1)) - mu - dot_product(data(i, 3:2 + p), theta(3:2 + p)))/sigma
         logf = -log(sigma) + z - exp(z)
         logs = -exp(z)
         value = value + data(i, 2)*logf + (1.0_dp - data(i, 2))*logs
      end do
   end function weibullregpost

   subroutine mnormt_onesided(m0, prior_mean, prior_sd, xbar, n, sampling_sd, bf, prior_odds, post_odds, post_h)
      real(dp), intent(in) :: m0 !! Mean value defining the one-sided null hypothesis mu <= m0.
      real(dp), intent(in) :: prior_mean !! Mean of the normal prior on mu.
      real(dp), intent(in) :: prior_sd !! Positive standard deviation of the normal prior on mu.
      real(dp), intent(in) :: xbar !! Observed sample mean.
      integer, intent(in) :: n !! Positive sample size associated with xbar.
      real(dp), intent(in) :: sampling_sd !! Known positive sampling standard deviation of individual observations.
      real(dp), intent(out) :: bf !! Bayes factor equal to posterior odds divided by prior odds.
      real(dp), intent(out) :: prior_odds !! Prior odds for mu <= m0 versus the complement.
      real(dp), intent(out) :: post_odds !! Posterior odds for mu <= m0 versus the complement.
      real(dp), intent(out) :: post_h !! Posterior probability that mu <= m0.
      real(dp) :: prior_h
      real(dp) :: post_mean
      real(dp) :: post_sd
      real(dp) :: post_precision

      prior_h = normal_cdf(m0, prior_mean, prior_sd)
      prior_odds = prior_h/(1.0_dp - prior_h)
      post_precision = 1.0_dp/(prior_sd*prior_sd) + real(n, dp)/(sampling_sd*sampling_sd)
      post_sd = sqrt(1.0_dp/post_precision)
      post_mean = (xbar*real(n, dp)/(sampling_sd*sampling_sd) + prior_mean/(prior_sd*prior_sd))/post_precision
      post_h = normal_cdf(m0, post_mean, post_sd)
      post_odds = post_h/(1.0_dp - post_h)
      bf = post_odds/prior_odds
   end subroutine mnormt_onesided

   subroutine mnormt_twosided(m0, prior_prob, prior_scale, xbar, n, sampling_sd, bf, post)
      real(dp), intent(in) :: m0 !! Point-null mean value being tested.
      real(dp), intent(in) :: prior_prob !! Prior probability assigned to the point-null hypothesis.
      real(dp), intent(in) :: prior_scale !! Positive normal scale t under the alternative hypothesis.
      real(dp), intent(in) :: xbar !! Observed sample mean.
      integer, intent(in) :: n !! Positive sample size associated with xbar.
      real(dp), intent(in) :: sampling_sd !! Known positive sampling standard deviation.
      real(dp), intent(out) :: bf !! Bayes factor for the point null versus the normal alternative.
      real(dp), intent(out) :: post !! Posterior probability assigned to the point-null hypothesis.
      real(dp) :: den
      real(dp) :: num
      real(dp) :: variance

      num = 0.5_dp*log(real(n, dp)) - log(sampling_sd) - &
         0.5_dp*real(n, dp)/(sampling_sd*sampling_sd)*(xbar - m0)**2
      variance = sampling_sd*sampling_sd/real(n, dp) + prior_scale*prior_scale
      den = -0.5_dp*log(variance) - 0.5_dp/variance*(xbar - m0)**2
      bf = exp(num - den)
      post = prior_prob*bf/(prior_prob*bf + 1.0_dp - prior_prob)
   end subroutine mnormt_twosided

end module learnbayes_models
