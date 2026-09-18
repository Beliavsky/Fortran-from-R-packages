! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_components
   use flexmix_kinds, only : dp
   use flexmix_numeric, only : weighted_least_squares, weighted_covariance, clip_probability, solve_linear_system
   use flexmix_numeric, only : normal_logpdf, poisson_logpmf, binomial_logpmf, mvnormal_logpdf_rows
   use flexmix_numeric, only : inverse_logdet_spd, symmetric_eigen_jacobi
   implicit none
   private
   public :: fit_gaussian_regression_component
   public :: fit_poisson_regression_component
   public :: fit_binomial_regression_component
   public :: fit_gamma_regression_component
   public :: fit_multinomial_regression_component
   public :: gaussian_regression_log_density
   public :: poisson_regression_log_density
   public :: binomial_regression_log_density
   public :: gamma_regression_log_density
   public :: multinomial_regression_probabilities
   public :: multinomial_regression_log_density
   public :: fit_conditional_logit_component
   public :: conditional_logit_log_density
   public :: fit_mvnormal_component
   public :: fit_factor_analysis_component
   public :: mvnormal_log_density
   public :: fit_mvbinary_component
   public :: mvbinary_log_density
   public :: fit_mvpois_component
   public :: mvpois_log_density
   public :: fit_mvcombi_component
   public :: mvcombi_log_density
   public :: fit_lognormal_component
   public :: lognormal_log_density
   public :: fit_exponential_component
   public :: exponential_log_density
   public :: fit_inverse_gaussian_component
   public :: inverse_gaussian_log_density
   public :: fit_gamma_component
   public :: gamma_log_density
   public :: fit_weibull_component
   public :: weibull_log_density

contains

   pure subroutine fit_gaussian_regression_component(x, y, weights, beta, sigma, df, info, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response vector, size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: beta(:) !! Fitted component regression coefficients, size `p`.
      real(dp), intent(out) :: sigma !! Fitted positive residual standard deviation.
      integer, intent(out) :: df !! Number of component parameters, matching `FLXMRglm` Gaussian semantics.
      integer, intent(out) :: info !! Zero on success; nonzero if weighted least squares is rank deficient.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: residuals(:), adjusted_y(:)
      real(dp) :: sw, rss
      integer :: rank_value, n
      n = size(x, 1)
      df = size(x, 2) + 1
      sigma = 1.0_dp
      allocate(residuals(n), adjusted_y(n))
      adjusted_y = y
      if (present(offset)) adjusted_y = adjusted_y - offset
      call weighted_least_squares(x, adjusted_y, weights, beta, rank_value, residuals, info)
      if (info /= 0) return
      sw = sum(weights)
      if (sw <= tiny(1.0_dp) .or. n <= rank_value) then
         info = 2
         return
      end if
      rss = sum(weights * residuals * residuals)
      sigma = sqrt(max(tiny(1.0_dp), rss * real(n, dp) / sw / real(n - rank_value, dp)))
   end subroutine fit_gaussian_regression_component

   pure subroutine gaussian_regression_log_density(x, y, beta, sigma, log_density, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix for density evaluation, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Gaussian response values, size `n`.
      real(dp), intent(in) :: beta(:) !! Component regression coefficients, size `p`.
      real(dp), intent(in) :: sigma !! Positive component residual standard deviation.
      real(dp), intent(out) :: log_density(:) !! Per-observation Gaussian log densities, size `n`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: mu(:)
      allocate(mu(size(y)))
      mu = matmul(x, beta)
      if (present(offset)) mu = mu + offset
      log_density = normal_logpdf(y, mu, sigma)
   end subroutine gaussian_regression_log_density

   pure subroutine fit_poisson_regression_component(x, y, weights, beta, df, info, start_beta, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson counts, size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: beta(:) !! Fitted log-link regression coefficients, size `p`.
      integer, intent(out) :: df !! Number of component regression parameters.
      integer, intent(out) :: info !! Zero on convergence; nonzero for rank or numerical failure.
      real(dp), intent(in), optional :: start_beta(:) !! Optional previous coefficient vector for IRLS initialization.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: eta(:), mu(:), work_weights(:), z(:), beta_new(:), residuals(:), ystart(:), off(:)
      real(dp) :: change
      integer :: iter, rank_value
      df = size(x, 2)
      info = 0
      allocate(eta(size(y)), mu(size(y)), work_weights(size(y)), z(size(y)))
      allocate(beta_new(size(beta)), residuals(size(y)), ystart(size(y)), off(size(y)))
      off = 0.0_dp
      if (present(offset)) off = offset
      if (present(start_beta)) then
         beta = start_beta
      else
         ystart = log(max(y, 0.0_dp) + 0.1_dp) - off
         call weighted_least_squares(x, ystart, weights, beta, rank_value, residuals, info)
         if (info /= 0) then
            beta = 0.0_dp
            info = 0
         end if
      end if
      do iter = 1, 60
         eta = max(-30.0_dp, min(30.0_dp, matmul(x, beta) + off))
         mu = exp(eta)
         work_weights = weights * mu
         z = eta + (y - mu) / max(mu, 1.0e-12_dp) - off
         call weighted_least_squares(x, z, work_weights, beta_new, rank_value, residuals, info)
         if (info /= 0) return
         change = maxval(abs(beta_new - beta))
         beta = beta_new
         if (change <= 1.0e-9_dp * (1.0_dp + maxval(abs(beta)))) exit
      end do
   end subroutine fit_poisson_regression_component

   pure subroutine poisson_regression_log_density(x, y, beta, log_density, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix for density evaluation, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Nonnegative Poisson counts, size `n`.
      real(dp), intent(in) :: beta(:) !! Component log-link coefficients, size `p`.
      real(dp), intent(out) :: log_density(:) !! Per-observation Poisson log probabilities, size `n`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: eta(:), lambda_value(:)
      allocate(eta(size(y)), lambda_value(size(y)))
      eta = matmul(x, beta)
      if (present(offset)) eta = eta + offset
      lambda_value = exp(max(-30.0_dp, min(30.0_dp, eta)))
      log_density = poisson_logpmf(y, lambda_value)
   end subroutine poisson_regression_log_density

   pure subroutine fit_binomial_regression_component(x, success, trials, weights, beta, df, info, start_beta, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`.
      real(dp), intent(in) :: success(:) !! Binomial success counts, size `n`.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals, size `n` and not smaller than successes.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: beta(:) !! Fitted logit-link regression coefficients, size `p`.
      integer, intent(out) :: df !! Number of component regression parameters.
      integer, intent(out) :: info !! Zero on convergence; nonzero for rank or numerical failure.
      real(dp), intent(in), optional :: start_beta(:) !! Optional previous coefficient vector for IRLS initialization.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: eta(:), prob(:), work_weights(:), z(:), beta_new(:), residuals(:), target(:), off(:)
      real(dp) :: change, denom, p0
      integer :: i, iter, rank_value
      df = size(x, 2)
      info = 0
      allocate(eta(size(success)), prob(size(success)), work_weights(size(success)), z(size(success)))
      allocate(beta_new(size(beta)), residuals(size(success)), target(size(success)), off(size(success)))
      off = 0.0_dp
      if (present(offset)) off = offset
      if (present(start_beta)) then
         beta = start_beta
      else
         do i = 1, size(success)
            denom = max(trials(i), 1.0_dp)
            p0 = clip_probability((success(i) + 0.5_dp) / (denom + 1.0_dp))
            target(i) = log(p0 / (1.0_dp - p0)) - off(i)
         end do
         call weighted_least_squares(x, target, weights * max(trials, 1.0_dp), beta, rank_value, residuals, info)
         if (info /= 0) then
            beta = 0.0_dp
            info = 0
         end if
      end if
      do iter = 1, 60
         eta = max(-30.0_dp, min(30.0_dp, matmul(x, beta) + off))
         prob = 1.0_dp / (1.0_dp + exp(-eta))
         work_weights = weights * trials * prob * (1.0_dp - prob)
         do i = 1, size(success)
            denom = max(trials(i) * prob(i) * (1.0_dp - prob(i)), 1.0e-12_dp)
            z(i) = eta(i) + (success(i) - trials(i) * prob(i)) / denom - off(i)
         end do
         call weighted_least_squares(x, z, work_weights, beta_new, rank_value, residuals, info)
         if (info /= 0) return
         change = maxval(abs(beta_new - beta))
         beta = beta_new
         if (change <= 1.0e-9_dp * (1.0_dp + maxval(abs(beta)))) exit
      end do
   end subroutine fit_binomial_regression_component

   pure subroutine binomial_regression_log_density(x, success, trials, beta, log_density, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix for density evaluation, shape `(n, p)`.
      real(dp), intent(in) :: success(:) !! Binomial success counts, size `n`.
      real(dp), intent(in) :: trials(:) !! Binomial trial totals, size `n`.
      real(dp), intent(in) :: beta(:) !! Component logit-link coefficients, size `p`.
      real(dp), intent(out) :: log_density(:) !! Per-observation binomial log probabilities, size `n`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive linear-predictor offset, size `n`.
      real(dp), allocatable :: eta(:), prob(:)
      allocate(eta(size(success)), prob(size(success)))
      eta = matmul(x, beta)
      if (present(offset)) eta = eta + offset
      eta = max(-30.0_dp, min(30.0_dp, eta))
      prob = 1.0_dp / (1.0_dp + exp(-eta))
      log_density = binomial_logpmf(success, trials, prob)
   end subroutine binomial_regression_log_density

   pure subroutine fit_gamma_regression_component(x, y, weights, beta, shape, df, info, start_beta, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Strictly positive Gamma responses, size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: beta(:) !! Fitted coefficients for R Gamma's default inverse link, size `p`.
      real(dp), intent(out) :: shape !! Gamma shape parameter estimated from FlexMix's weighted deviance rule.
      integer, intent(out) :: df !! Number of component parameters, including the Gamma shape parameter.
      integer, intent(out) :: info !! Zero on convergence; nonzero for invalid data or numerical failure.
      real(dp), intent(in), optional :: start_beta(:) !! Optional previous coefficient vector for IRLS initialization.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive inverse-link predictor offset, size `n`.
      real(dp), allocatable :: eta(:), mu(:), work_weights(:), z(:), beta_new(:), residuals(:), target(:), off(:)
      real(dp) :: change, deviance, sw, term
      integer :: i, iter, rank_value
      df = size(x, 2) + 1
      info = 0
      shape = 1.0_dp
      if (any(y <= 0.0_dp)) then
         info = 1
         beta = 0.0_dp
         return
      end if
      allocate(eta(size(y)), mu(size(y)), work_weights(size(y)), z(size(y)))
      allocate(beta_new(size(beta)), residuals(size(y)), target(size(y)), off(size(y)))
      off = 0.0_dp
      if (present(offset)) off = offset
      if (present(start_beta)) then
         beta = start_beta
      else
         target = 1.0_dp / y - off
         call weighted_least_squares(x, target, weights, beta, rank_value, residuals, info)
         if (info /= 0) then
            beta = 0.0_dp
            if (size(beta) > 0) beta(1) = 1.0_dp / max(sum(weights * y) / max(sum(weights), tiny(1.0_dp)), 1.0e-6_dp)
            info = 0
         end if
      end if
      do iter = 1, 80
         eta = matmul(x, beta) + off
         eta = max(1.0e-8_dp, min(1.0e8_dp, eta))
         mu = 1.0_dp / eta
         work_weights = weights * mu * mu
         z = eta - (y - mu) / max(mu * mu, 1.0e-20_dp) - off
         call weighted_least_squares(x, z, work_weights, beta_new, rank_value, residuals, info)
         if (info /= 0) return
         change = maxval(abs(beta_new - beta))
         beta = beta_new
         if (change <= 1.0e-9_dp * (1.0_dp + maxval(abs(beta)))) exit
      end do
      eta = max(1.0e-8_dp, min(1.0e8_dp, matmul(x, beta) + off))
      mu = 1.0_dp / eta
      sw = sum(weights)
      if (sw <= tiny(1.0_dp)) then
         info = 2
         return
      end if
      deviance = 0.0_dp
      do i = 1, size(y)
         if (weights(i) <= 0.0_dp) cycle
         term = (y(i) - mu(i)) / mu(i) - log(y(i) / mu(i))
         deviance = deviance + 2.0_dp * weights(i) * term
      end do
      shape = sw / max(deviance, 1.0e-12_dp * sw)
   end subroutine fit_gamma_regression_component

   pure subroutine gamma_regression_log_density(x, y, beta, shape, log_density, offset)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix for density evaluation, shape `(n, p)`.
      real(dp), intent(in) :: y(:) !! Strictly positive Gamma response values, size `n`.
      real(dp), intent(in) :: beta(:) !! Component inverse-link coefficients, size `p`.
      real(dp), intent(in) :: shape !! Positive Gamma shape parameter.
      real(dp), intent(out) :: log_density(:) !! Per-observation Gamma log densities, size `n`.
      real(dp), intent(in), optional :: offset(:) !! Optional fixed additive inverse-link predictor offset, size `n`.
      real(dp), allocatable :: eta(:), mu(:), rate_value(:)
      integer :: i
      allocate(eta(size(y)), mu(size(y)), rate_value(size(y)))
      eta = matmul(x, beta)
      if (present(offset)) eta = eta + offset
      eta = max(1.0e-8_dp, min(1.0e8_dp, eta))
      mu = 1.0_dp / eta
      rate_value = shape / mu
      do i = 1, size(y)
         if (y(i) <= 0.0_dp .or. shape <= 0.0_dp) then
            log_density(i) = -huge(1.0_dp)
         else
            log_density(i) = shape * log(rate_value(i)) - log_gamma(shape) + &
                             (shape - 1.0_dp) * log(y(i)) - rate_value(i) * y(i)
         end if
      end do
   end subroutine gamma_regression_log_density

   pure subroutine fit_multinomial_regression_component(x, y_class, nclass, weights, coef, df, info, start_coef)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`.
      integer, intent(in) :: y_class(:) !! One-based response class labels, size `n`, covering `1:nclass` in valid use.
      integer, intent(in) :: nclass !! Number of response categories, at least two.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: coef(:,:) !! Fitted nonreference coefficients, shape `(p, nclass-1)`; class one is the reference.
      integer, intent(out) :: df !! Number of free multinomial regression coefficients.
      integer, intent(out) :: info !! Zero on convergence; nonzero for invalid labels or a singular Newton system.
      real(dp), intent(in), optional :: start_coef(:,:) !! Optional previous coefficient matrix, shape `(p, nclass-1)`.
      real(dp), allocatable :: prob(:,:), grad(:), hess(:,:), delta(:), old_coef(:,:)
      real(dp) :: block_weight, change, ridge, wi
      integer :: a, b, c, d, i, iter, p, q, s, t, solve_info
      p = size(x, 2)
      q = p * max(0, nclass - 1)
      df = q
      info = 0
      coef = 0.0_dp
      if (nclass < 2 .or. size(coef,1) /= p .or. size(coef,2) /= nclass - 1 .or. &
          any(y_class < 1) .or. any(y_class > nclass)) then
         info = 1
         return
      end if
      if (present(start_coef)) coef = start_coef
      allocate(prob(size(x,1),nclass), grad(q), hess(q,q), delta(q), old_coef(p,nclass-1))
      do iter = 1, 50
         old_coef = coef
         call multinomial_regression_probabilities(x, coef, prob)
         grad = 0.0_dp
         hess = 0.0_dp
         do i = 1, size(x,1)
            wi = weights(i)
            if (wi <= 0.0_dp) cycle
            do a = 2, nclass
               c = a - 1
               do s = 1, p
                  if (y_class(i) == a) then
                     grad((c - 1) * p + s) = grad((c - 1) * p + s) + wi * x(i,s) * (1.0_dp - prob(i,a))
                  else
                     grad((c - 1) * p + s) = grad((c - 1) * p + s) - wi * x(i,s) * prob(i,a)
                  end if
               end do
               do b = 2, nclass
                  d = b - 1
                  if (a == b) then
                     block_weight = prob(i,a) * (1.0_dp - prob(i,a))
                  else
                     block_weight = -prob(i,a) * prob(i,b)
                  end if
                  do s = 1, p
                     do t = 1, p
                        hess((c - 1) * p + s, (d - 1) * p + t) = hess((c - 1) * p + s, (d - 1) * p + t) + &
                           wi * block_weight * x(i,s) * x(i,t)
                     end do
                  end do
               end do
            end do
         end do
         ridge = 1.0e-8_dp * max(1.0_dp, maxval(abs(hess)))
         do s = 1, q
            hess(s,s) = hess(s,s) + ridge
         end do
         call solve_linear_system(hess, grad, delta, solve_info)
         if (solve_info /= 0) then
            info = solve_info
            return
         end if
         do a = 2, nclass
            c = a - 1
            coef(:,c) = coef(:,c) + delta((c - 1) * p + 1:c * p)
         end do
         change = maxval(abs(coef - old_coef))
         if (change <= 1.0e-8_dp * (1.0_dp + maxval(abs(coef)))) exit
      end do
   end subroutine fit_multinomial_regression_component

   pure subroutine multinomial_regression_probabilities(x, coef, probability)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix, shape `(n, p)`.
      real(dp), intent(in) :: coef(:,:) !! Nonreference coefficient matrix, shape `(p, nclass-1)`.
      real(dp), intent(out) :: probability(:,:) !! Row-wise class probabilities `(n,nclass)` with class one as reference.
      real(dp), allocatable :: eta(:)
      real(dp) :: m, total
      integer :: a, i, nclass
      nclass = size(coef,2) + 1
      allocate(eta(nclass))
      do i = 1, size(x,1)
         eta(1) = 0.0_dp
         do a = 2, nclass
            eta(a) = dot_product(x(i,:), coef(:,a-1))
         end do
         m = maxval(eta)
         total = 0.0_dp
         do a = 1, nclass
            probability(i,a) = exp(max(-700.0_dp, min(700.0_dp, eta(a) - m)))
            total = total + probability(i,a)
         end do
         probability(i,:) = probability(i,:) / total
      end do
   end subroutine multinomial_regression_probabilities

   pure subroutine multinomial_regression_log_density(x, y_class, coef, log_density)
      real(dp), intent(in) :: x(:,:) !! Regression design matrix for density evaluation, shape `(n, p)`.
      integer, intent(in) :: y_class(:) !! One-based response class labels, size `n`.
      real(dp), intent(in) :: coef(:,:) !! Nonreference coefficient matrix, shape `(p, nclass-1)`.
      real(dp), intent(out) :: log_density(:) !! Per-observation log probability of the observed response class, size `n`.
      real(dp), allocatable :: probability(:,:)
      integer :: i, nclass
      nclass = size(coef,2) + 1
      allocate(probability(size(x,1),nclass))
      call multinomial_regression_probabilities(x, coef, probability)
      do i = 1, size(y_class)
         if (y_class(i) < 1 .or. y_class(i) > nclass) then
            log_density(i) = -huge(1.0_dp)
         else
            log_density(i) = log(max(probability(i,y_class(i)), tiny(1.0_dp)))
         end if
      end do
   end subroutine multinomial_regression_log_density

   pure subroutine fit_conditional_logit_component(x, y, strata, weights, beta, df, info, initial_beta)
      real(dp), intent(in) :: x(:,:) !! Conditional-logit covariate matrix without an intercept, shape `(n,p)`.
      real(dp), intent(in) :: y(:) !! Binary choice indicator, size `n`; exactly one row must equal one in every stratum.
      integer, intent(in) :: strata(:) !! Integer choice-set labels, size `n`; equal labels define one conditional-logit.
      real(dp), intent(in) :: weights(:) !! Nonnegative M-step weights, size `n`; weights must be constant within each stratum.
      real(dp), intent(out) :: beta(:) !! Fitted conditional-logit coefficients, size `p`.
      integer, intent(out) :: df !! Number of fitted coefficients, equal to `p`.
      integer, intent(out) :: info !! Zero on convergence; nonzero for invalid strata, weights, or Newton-system failure.
      real(dp), intent(in), optional :: initial_beta(:) !! Optional starting coefficient vector, size `p`; defaults to zero.
      real(dp), allocatable :: score(:), information(:,:), system(:,:), delta(:), candidate(:), probability(:)
      real(dp), allocatable :: mean_x(:), chosen_x(:), eta(:)
      real(dp) :: llh, llh_new, max_eta, denom, alpha, ridge, group_weight, scale, tolerance
      integer :: i, j, a, b, iter, line, n, p, group_size, event_count, solve_info
      logical :: is_first

      n = size(y)
      p = size(x,2)
      df = p
      beta = 0.0_dp
      info = 0
      if (size(x,1) /= n .or. size(strata) /= n .or. size(weights) /= n .or. size(beta) /= p .or. &
          n < 1 .or. p < 1 .or. any(weights < 0.0_dp)) then
         info = 1
         return
      end if
      if (any(abs(y) > 100.0_dp * epsilon(1.0_dp) .and. &
              abs(y - 1.0_dp) > 100.0_dp * epsilon(1.0_dp))) then
         info = 2
         return
      end if
      if (present(initial_beta)) then
         if (size(initial_beta) /= p) then
            info = 3
            return
         end if
         beta = initial_beta
      end if
      do i = 1, n
         is_first = .true.
         do j = 1, i - 1
            if (strata(j) == strata(i)) then
               is_first = .false.
               exit
            end if
         end do
         if (.not. is_first) cycle
         group_size = 0
         event_count = 0
         group_weight = weights(i)
         do j = 1, n
            if (strata(j) /= strata(i)) cycle
            group_size = group_size + 1
            if (abs(y(j) - 1.0_dp) <= 100.0_dp * epsilon(1.0_dp)) event_count = event_count + 1
            tolerance = 1.0e-10_dp * (1.0_dp + max(abs(group_weight),abs(weights(j))))
            if (abs(weights(j) - group_weight) > tolerance) then
               info = 4
               return
            end if
         end do
         if (group_size < 2 .or. event_count /= 1) then
            info = 5
            return
         end if
      end do
      allocate(score(p), information(p,p), system(p,p), delta(p), candidate(p), probability(n))
      allocate(mean_x(p), chosen_x(p), eta(n))
      llh = conditional_logit_weighted_loglik(x, y, strata, weights, beta)
      do iter = 1, 100
         eta = matmul(x, beta)
         score = 0.0_dp
         information = 0.0_dp
         do i = 1, n
            is_first = .true.
            do j = 1, i - 1
               if (strata(j) == strata(i)) then
                  is_first = .false.
                  exit
               end if
            end do
            if (.not. is_first) cycle
            group_weight = weights(i)
            if (group_weight <= tiny(1.0_dp)) cycle
            max_eta = -huge(1.0_dp)
            do j = 1, n
               if (strata(j) == strata(i)) max_eta = max(max_eta, eta(j))
            end do
            denom = 0.0_dp
            do j = 1, n
               if (strata(j) == strata(i)) denom = denom + exp(eta(j) - max_eta)
            end do
            mean_x = 0.0_dp
            chosen_x = 0.0_dp
            probability = 0.0_dp
            do j = 1, n
               if (strata(j) /= strata(i)) cycle
               probability(j) = exp(eta(j) - max_eta) / denom
               mean_x = mean_x + probability(j) * x(j,:)
               if (y(j) > 0.5_dp) chosen_x = chosen_x + x(j,:)
            end do
            score = score + group_weight * (chosen_x - mean_x)
            do j = 1, n
               if (strata(j) /= strata(i)) cycle
               do a = 1, p
                  do b = 1, p
                     information(a,b) = information(a,b) + group_weight * probability(j) * &
                                        (x(j,a) - mean_x(a)) * (x(j,b) - mean_x(b))
                  end do
               end do
            end do
         end do
         if (maxval(abs(score)) <= 1.0e-9_dp) exit
         scale = max(1.0_dp, maxval(abs(information)))
         ridge = 1.0e-10_dp * scale
         solve_info = 1
         do a = 1, 8
            system = information
            do b = 1, p
               system(b,b) = system(b,b) + ridge
            end do
            call solve_linear_system(system, score, delta, solve_info)
            if (solve_info == 0) exit
            ridge = 10.0_dp * ridge
         end do
         if (solve_info /= 0) then
            info = 6
            return
         end if
         alpha = 1.0_dp
         do line = 1, 40
            candidate = beta + alpha * delta
            llh_new = conditional_logit_weighted_loglik(x, y, strata, weights, candidate)
            if (llh_new >= llh + 1.0e-4_dp * alpha * dot_product(score, delta)) exit
            alpha = 0.5_dp * alpha
         end do
         if (line > 40) then
            info = 7
            return
         end if
         beta = candidate
         llh = llh_new
         if (maxval(abs(alpha * delta)) <= 1.0e-9_dp * (1.0_dp + maxval(abs(beta)))) exit
      end do
      info = 0
   end subroutine fit_conditional_logit_component

   pure subroutine conditional_logit_log_density(x, y, strata, beta, log_density, info)
      real(dp), intent(in) :: x(:,:) !! Conditional-logit covariate matrix without an intercept, shape `(n,p)`.
      real(dp), intent(in) :: y(:) !! Binary choice indicator, size `n`; exactly one selected row is expected in each stratum.
      integer, intent(in) :: strata(:) !! Integer choice-set labels, size `n`.
      real(dp), intent(in) :: beta(:) !! Conditional-logit coefficient vector, size `p`.
      real(dp), intent(out) :: log_density(:) !! Per-row stratum log likelihood divided by stratum size, matching upstream.
      integer, intent(out) :: info !! Zero on success; nonzero for incompatible dimensions or malformed strata.
      real(dp), allocatable :: eta(:)
      real(dp) :: max_eta, denom, event_eta, value
      integer :: i, j, n, group_size, event_count
      logical :: is_first

      n = size(y)
      log_density = 0.0_dp
      info = 0
      if (size(x,1) /= n .or. size(strata) /= n .or. size(log_density) /= n .or. size(x,2) /= size(beta)) then
         info = 1
         return
      end if
      allocate(eta(n))
      eta = matmul(x, beta)
      do i = 1, n
         is_first = .true.
         do j = 1, i - 1
            if (strata(j) == strata(i)) then
               is_first = .false.
               exit
            end if
         end do
         if (.not. is_first) cycle
         max_eta = -huge(1.0_dp)
         event_eta = 0.0_dp
         group_size = 0
         event_count = 0
         do j = 1, n
            if (strata(j) /= strata(i)) cycle
            max_eta = max(max_eta, eta(j))
            event_eta = event_eta + y(j) * eta(j)
            group_size = group_size + 1
            if (abs(y(j) - 1.0_dp) <= 100.0_dp * epsilon(1.0_dp)) event_count = event_count + 1
         end do
         if (group_size < 2 .or. event_count /= 1) then
            info = 2
            return
         end if
         denom = 0.0_dp
         do j = 1, n
            if (strata(j) == strata(i)) denom = denom + exp(eta(j) - max_eta)
         end do
         value = (event_eta - max_eta - log(denom)) / real(group_size,dp)
         do j = 1, n
            if (strata(j) == strata(i)) log_density(j) = value
         end do
      end do
   end subroutine conditional_logit_log_density

   pure function conditional_logit_weighted_loglik(x, y, strata, weights, beta) result(value)
      real(dp), intent(in) :: x(:,:) !! Conditional-logit covariate matrix without an intercept, shape `(n,p)`.
      real(dp), intent(in) :: y(:) !! Binary choice indicator, size `n`, with one selected alternative per stratum.
      integer, intent(in) :: strata(:) !! Integer choice-set labels, size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative stratum-constant case weights, size `n`.
      real(dp), intent(in) :: beta(:) !! Conditional-logit coefficient vector, size `p`.
      real(dp) :: value
      real(dp), allocatable :: eta(:)
      real(dp) :: max_eta, denom, event_eta
      integer :: i, j, n
      logical :: is_first

      n = size(y)
      allocate(eta(n))
      eta = matmul(x, beta)
      value = 0.0_dp
      do i = 1, n
         is_first = .true.
         do j = 1, i - 1
            if (strata(j) == strata(i)) then
               is_first = .false.
               exit
            end if
         end do
         if (.not. is_first .or. weights(i) <= tiny(1.0_dp)) cycle
         max_eta = -huge(1.0_dp)
         event_eta = 0.0_dp
         do j = 1, n
            if (strata(j) /= strata(i)) cycle
            max_eta = max(max_eta, eta(j))
            event_eta = event_eta + y(j) * eta(j)
         end do
         denom = 0.0_dp
         do j = 1, n
            if (strata(j) == strata(i)) denom = denom + exp(eta(j) - max_eta)
         end do
         value = value + weights(i) * (event_eta - max_eta - log(denom))
      end do
   end function conditional_logit_weighted_loglik

   pure subroutine fit_mvnormal_component(y, weights, center, covariance, diagonal, df, info)
      real(dp), intent(in) :: y(:,:) !! Multivariate observations, shape `(n, d)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: center(:) !! Weighted component center, size `d`.
      real(dp), intent(out) :: covariance(:,:) !! Weighted covariance matrix, shape `(d, d)`.
      logical, intent(in) :: diagonal !! If true, zero off-diagonal covariance entries as in `FLXMCmvnorm`.
      integer, intent(out) :: df !! Number of component parameters used by upstream FlexMix.
      integer, intent(out) :: info !! Zero on success; nonzero if weighted covariance cannot be estimated.
      logical :: ok
      integer :: i, j, d
      d = size(y, 2)
      call weighted_covariance(y, weights, center, covariance, ok)
      if (.not. ok) then
         info = 1
         return
      end if
      info = 0
      if (diagonal) then
         do j = 1, d
            do i = 1, d
               if (i /= j) covariance(i,j) = 0.0_dp
            end do
            covariance(j,j) = max(covariance(j,j), 1.0e-12_dp)
         end do
         df = 2 * d
      else
         do j = 1, d
            covariance(j,j) = max(covariance(j,j), 1.0e-12_dp)
         end do
         df = (3 * d + d * d) / 2
      end if
   end subroutine fit_mvnormal_component

   pure subroutine fit_factor_analysis_component(y, weights, factors, center, covariance, marginal_variance, &
                                                    loadings, uniqueness, df, info)
      real(dp), intent(in) :: y(:,:) !! Multivariate observations, shape `(n,d)`, for one weighted factor-analyzer component.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`, interpreted with `cov.wt` semantics.
      integer, intent(in) :: factors !! Number of latent factors; must satisfy the upstream factor-analysis degrees-of-freedom.
      real(dp), intent(out) :: center(:) !! Weighted component mean vector, size `d`.
      real(dp), intent(out) :: covariance(:,:) !! Fitted factor-model covariance on the original variable scale, shape `(d,d)`.
      real(dp), intent(out) :: marginal_variance(:) !! Diagonal of the unconstrained weighted covariance, size `d`.
      real(dp), intent(out) :: loadings(:,:) !! Maximum-likelihood correlation-scale factor loadings `(d,factors)`.
      real(dp), intent(out) :: uniqueness(:) !! Correlation-scale uniqueness parameters constrained to `[0.005,1]`, size `d`.
      integer, intent(out) :: df !! Upstream FlexMix component parameter count `(factors+2)*d`.
      integer, intent(out) :: info !! Zero on success; nonzero for dimensions, degenerate covariance, or failed optimization.
      real(dp), allocatable :: raw_covariance(:,:), correlation(:,:), inverse(:,:), sds(:), gradient(:), model_correlation(:,:)
      real(dp) :: logdet, objective
      logical :: ok
      integer :: d, i, inverse_info, opt_info

      d = size(y,2)
      df = (factors + 2) * d
      center = 0.0_dp
      covariance = 0.0_dp
      marginal_variance = 0.0_dp
      loadings = 0.0_dp
      uniqueness = 1.0_dp
      if (size(y,1) /= size(weights) .or. size(center) /= d .or. any(shape(covariance) /= [d,d]) .or. &
          size(marginal_variance) /= d .or. size(uniqueness) /= d .or. size(loadings,1) /= d .or. &
          size(loadings,2) /= factors .or. d < 3 .or. factors < 1) then
         info = 1
         return
      end if
      if ((d - factors)**2 - d - factors < 0) then
         info = 2
         return
      end if
      allocate(raw_covariance(d,d), correlation(d,d), inverse(d,d), sds(d), gradient(d), model_correlation(d,d))
      call weighted_covariance(y, weights, center, raw_covariance, ok)
      if (.not. ok) then
         info = 3
         return
      end if
      do i = 1, d
         marginal_variance(i) = raw_covariance(i,i)
      end do
      if (any(marginal_variance <= 100.0_dp * tiny(1.0_dp))) then
         info = 4
         return
      end if
      sds = sqrt(marginal_variance)
      do i = 1, d
         correlation(i,:) = raw_covariance(i,:) / (sds(i) * sds)
      end do
      call inverse_logdet_spd(correlation, inverse, logdet, inverse_info)
      if (inverse_info == 0) then
         do i = 1, d
            uniqueness(i) = (1.0_dp - 0.5_dp * real(factors,dp) / real(d,dp)) / max(inverse(i,i), tiny(1.0_dp))
         end do
      else
         uniqueness = 0.5_dp
      end if
      uniqueness = max(0.005_dp, min(1.0_dp, uniqueness))
      call optimize_factor_uniqueness(correlation, factors, uniqueness, loadings, objective, opt_info)
      if (opt_info /= 0) then
         info = 5
         return
      end if
      call factor_objective_gradient(uniqueness, correlation, factors, objective, gradient, loadings, opt_info)
      if (opt_info /= 0) then
         info = 6
         return
      end if
      model_correlation = matmul(loadings, transpose(loadings))
      do i = 1, d
         model_correlation(i,i) = model_correlation(i,i) + uniqueness(i)
      end do
      do i = 1, d
         covariance(i,:) = model_correlation(i,:) * sds(i) * sds
      end do
      info = 0
   end subroutine fit_factor_analysis_component

   pure subroutine optimize_factor_uniqueness(correlation, factors, uniqueness, loadings, objective, info)
      real(dp), intent(in) :: correlation(:,:) !! Positive-definite correlation matrix, shape `(d,d)`.
      integer, intent(in) :: factors !! Number of latent factors retained in the maximum-likelihood objective.
      real(dp), intent(inout) :: uniqueness(:) !! Starting uniquenesses on entry and bounded optimized uniquenesses on return.
      real(dp), intent(out) :: loadings(:,:) !! Correlation-scale loadings at optimized uniquenesses `(d,factors)`.
      real(dp), intent(out) :: objective !! Final `factanal.fit.mle` objective value.
      integer, intent(out) :: info !! Zero on a usable optimum; nonzero if the objective cannot be evaluated or line search.
      real(dp), allocatable :: gradient(:), gradient_new(:), projected(:), direction(:), candidate(:), step_vector(:)
      real(dp), allocatable :: delta_gradient(:)
      real(dp), allocatable :: hessian_inverse(:,:), identity(:,:), left(:,:), right(:,:), trial_loadings(:,:)
      real(dp) :: objective_new, alpha, descent, ys, rho, old_objective
      integer :: d, i, iter, eval_info, line

      d = size(uniqueness)
      allocate(gradient(d), gradient_new(d), projected(d), direction(d), candidate(d), step_vector(d))
      allocate(delta_gradient(d))
      allocate(hessian_inverse(d,d), identity(d,d), left(d,d), right(d,d), trial_loadings(d,factors))
      identity = 0.0_dp
      do i = 1, d
         identity(i,i) = 1.0_dp
      end do
      hessian_inverse = identity
      call factor_objective_gradient(uniqueness, correlation, factors, objective, gradient, loadings, eval_info)
      if (eval_info /= 0) then
         info = 1
         return
      end if
      do iter = 1, 250
         projected = gradient
         do i = 1, d
            if (uniqueness(i) <= 0.005_dp + 1.0e-10_dp .and. projected(i) > 0.0_dp) projected(i) = 0.0_dp
            if (uniqueness(i) >= 1.0_dp - 1.0e-10_dp .and. projected(i) < 0.0_dp) projected(i) = 0.0_dp
         end do
         if (maxval(abs(projected)) <= 2.0e-7_dp) exit
         direction = -matmul(hessian_inverse, projected)
         descent = dot_product(gradient, direction)
         if (descent >= -1.0e-14_dp * max(1.0_dp, sqrt(dot_product(direction,direction)))) then
            direction = -projected
            hessian_inverse = identity
         end if
         alpha = 1.0_dp
         old_objective = objective
         do line = 1, 40
            candidate = max(0.005_dp, min(1.0_dp, uniqueness + alpha * direction))
            step_vector = candidate - uniqueness
            if (maxval(abs(step_vector)) <= 1.0e-13_dp) then
               alpha = 0.5_dp * alpha
               cycle
            end if
            call factor_objective_gradient(candidate, correlation, factors, objective_new, gradient_new, trial_loadings, eval_info)
            if (eval_info == 0) then
               if (objective_new <= objective + 1.0e-4_dp * dot_product(gradient, step_vector)) exit
            end if
            alpha = 0.5_dp * alpha
         end do
         if (line > 40) then
            direction = -projected
            hessian_inverse = identity
            alpha = 1.0_dp
            do line = 1, 60
               candidate = max(0.005_dp, min(1.0_dp, uniqueness + alpha * direction))
               step_vector = candidate - uniqueness
               if (maxval(abs(step_vector)) <= 1.0e-13_dp) then
                  alpha = 0.5_dp * alpha
                  cycle
               end if
               call factor_objective_gradient(candidate, correlation, factors, objective_new, gradient_new, &
                                              trial_loadings, eval_info)
               if (eval_info == 0) then
                  if (objective_new <= objective + 1.0e-4_dp * dot_product(gradient, step_vector)) exit
               end if
               alpha = 0.5_dp * alpha
            end do
            if (line > 60) then
               info = 2
               return
            end if
         end if
         delta_gradient = gradient_new - gradient
         ys = dot_product(delta_gradient, step_vector)
         if (ys > 1.0e-12_dp * max(1.0_dp, sqrt(dot_product(delta_gradient,delta_gradient)) * &
                                      sqrt(dot_product(step_vector,step_vector)))) then
            rho = 1.0_dp / ys
            left = identity - rho * spread(step_vector,2,d) * spread(delta_gradient,1,d)
            right = identity - rho * spread(delta_gradient,2,d) * spread(step_vector,1,d)
            hessian_inverse = matmul(left, matmul(hessian_inverse, right)) + &
                              rho * spread(step_vector,2,d) * spread(step_vector,1,d)
         else
            hessian_inverse = identity
         end if
         uniqueness = candidate
         objective = objective_new
         gradient = gradient_new
         loadings = trial_loadings
         if (abs(old_objective - objective) <= 1.0e-11_dp * (1.0_dp + abs(objective))) exit
      end do
      info = 0
   end subroutine optimize_factor_uniqueness

   pure subroutine factor_objective_gradient(uniqueness, correlation, factors, objective, gradient, loadings, info)
      real(dp), intent(in) :: uniqueness(:) !! Positive correlation-scale uniqueness parameters, size `d`.
      real(dp), intent(in) :: correlation(:,:) !! Correlation matrix being factorized, shape `(d,d)`.
      integer, intent(in) :: factors !! Number of leading factors retained.
      real(dp), intent(out) :: objective !! Negative concentrated maximum-likelihood criterion used by `factanal.fit.mle`.
      real(dp), intent(out) :: gradient(:) !! Analytic gradient of the criterion with respect to uniquenesses, size `d`.
      real(dp), intent(out) :: loadings(:,:) !! Correlation-scale loadings implied by uniquenesses `(d,factors)`.
      integer, intent(out) :: info !! Zero when the eigensystem and objective are finite; nonzero otherwise.
      real(dp), allocatable :: scaled(:,:), eigenvalues(:), eigenvectors(:,:), model(:,:), scale(:)
      real(dp) :: value
      integer :: d, eigen_info, i, j

      d = size(uniqueness)
      objective = huge(1.0_dp)
      gradient = 0.0_dp
      loadings = 0.0_dp
      if (any(shape(correlation) /= [d,d]) .or. size(gradient) /= d .or. size(loadings,1) /= d .or. &
          size(loadings,2) /= factors .or. factors < 1 .or. factors >= d .or. any(uniqueness <= 0.0_dp)) then
         info = 1
         return
      end if
      allocate(scaled(d,d), eigenvalues(d), eigenvectors(d,d), model(d,d), scale(d))
      scale = 1.0_dp / sqrt(uniqueness)
      do i = 1, d
         scaled(i,:) = correlation(i,:) * scale(i) * scale
      end do
      call symmetric_eigen_jacobi(scaled, eigenvalues, eigenvectors, eigen_info)
      if (eigen_info /= 0 .or. any(eigenvalues(factors+1:d) <= 0.0_dp)) then
         info = 2
         return
      end if
      value = 0.0_dp
      do i = factors + 1, d
         value = value + log(eigenvalues(i)) - eigenvalues(i)
      end do
      objective = -value + real(factors-d,dp)
      do j = 1, factors
         loadings(:,j) = eigenvectors(:,j) * sqrt(max(eigenvalues(j) - 1.0_dp, 0.0_dp)) * sqrt(uniqueness)
      end do
      model = matmul(loadings, transpose(loadings))
      do i = 1, d
         model(i,i) = model(i,i) + uniqueness(i)
         gradient(i) = (model(i,i) - correlation(i,i)) / (uniqueness(i) * uniqueness(i))
      end do
      info = 0
   end subroutine factor_objective_gradient

   pure subroutine mvnormal_log_density(y, center, covariance, log_density, info)
      real(dp), intent(in) :: y(:,:) !! Multivariate observations, shape `(n, d)`.
      real(dp), intent(in) :: center(:) !! Component center, size `d`.
      real(dp), intent(in) :: covariance(:,:) !! Component covariance matrix, shape `(d, d)`.
      real(dp), intent(out) :: log_density(:) !! Per-observation multivariate-normal log densities.
      integer, intent(out) :: info !! Zero when covariance factorization succeeds.
      call mvnormal_logpdf_rows(y, center, covariance, log_density, info)
   end subroutine mvnormal_log_density

   pure subroutine fit_mvbinary_component(y, weights, probability, truncated, df, info)
      real(dp), intent(in) :: y(:,:) !! Binary observation matrix, shape `(n, d)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: probability(:) !! Fitted Bernoulli probabilities, size `d`.
      logical, intent(in) :: truncated !! If true, condition on at least one success in each row.
      integer, intent(out) :: df !! Number of Bernoulli component parameters, equal to `d`.
      integer, intent(out) :: info !! Zero on success; nonzero when total weight is zero.
      real(dp), allocatable :: rk(:), pnew(:)
      real(dp) :: r0, prodp, llh, llh_old, sw
      integer :: iter
      df = size(y, 2)
      sw = sum(weights)
      if (sw <= tiny(1.0_dp)) then
         probability = 0.5_dp
         info = 1
         return
      end if
      allocate(rk(size(probability)), pnew(size(probability)))
      rk = matmul(transpose(y), weights) / sw
      if (.not. truncated) then
         probability = clip_probability(rk)
         info = 0
         return
      end if
      r0 = 0.0_dp
      llh_old = -huge(1.0_dp)
      probability = clip_probability(rk)
      do iter = 1, 200
         pnew = clip_probability(rk / (1.0_dp + r0))
         llh = sum(merge(rk * log(pnew), 0.0_dp, rk > 0.0_dp)) + &
               sum(merge((1.0_dp - rk + r0) * log(1.0_dp - pnew), 0.0_dp, 1.0_dp - rk + r0 > 0.0_dp))
         probability = pnew
         if (abs(llh - llh_old) / (abs(llh) + 0.1_dp) < epsilon(1.0_dp)) exit
         llh_old = llh
         prodp = product(1.0_dp - probability)
         if (prodp >= 1.0_dp - 1.0e-15_dp) exit
         r0 = prodp / max(1.0e-15_dp, 1.0_dp - prodp)
      end do
      info = 0
   end subroutine fit_mvbinary_component

   pure subroutine mvbinary_log_density(y, probability, truncated, log_density)
      real(dp), intent(in) :: y(:,:) !! Binary observation matrix, shape `(n, d)`.
      real(dp), intent(in) :: probability(:) !! Component Bernoulli probabilities, size `d`.
      logical, intent(in) :: truncated !! If true, condition on rows with at least one success.
      real(dp), intent(out) :: log_density(:) !! Per-row Bernoulli product log probabilities.
      real(dp) :: normalizer, p
      integer :: i, j
      normalizer = 0.0_dp
      if (truncated) normalizer = log(max(1.0e-300_dp, 1.0_dp - product(1.0_dp - clip_probability(probability))))
      do i = 1, size(y, 1)
         log_density(i) = -normalizer
         do j = 1, size(y, 2)
            p = clip_probability(probability(j))
            log_density(i) = log_density(i) + y(i,j) * log(p) + (1.0_dp - y(i,j)) * log(1.0_dp - p)
         end do
      end do
   end subroutine mvbinary_log_density

   pure subroutine fit_mvpois_component(y, weights, lambda_value, df, info)
      real(dp), intent(in) :: y(:,:) !! Multivariate nonnegative count matrix, shape `(n, d)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: lambda_value(:) !! Fitted independent Poisson means, size `d`.
      integer, intent(out) :: df !! Number of Poisson mean parameters, equal to `d`.
      integer, intent(out) :: info !! Zero on success; nonzero when total weight is zero.
      real(dp) :: sw
      sw = sum(weights)
      df = size(y, 2)
      if (sw <= tiny(1.0_dp)) then
         lambda_value = 1.0_dp
         info = 1
         return
      end if
      lambda_value = max(1.0e-12_dp, matmul(transpose(y), weights) / sw)
      info = 0
   end subroutine fit_mvpois_component

   pure subroutine mvpois_log_density(y, lambda_value, log_density)
      real(dp), intent(in) :: y(:,:) !! Multivariate nonnegative counts, shape `(n, d)`.
      real(dp), intent(in) :: lambda_value(:) !! Independent Poisson means, size `d`.
      real(dp), intent(out) :: log_density(:) !! Per-row sum of Poisson log probabilities.
      integer :: i, j
      do i = 1, size(y, 1)
         log_density(i) = 0.0_dp
         do j = 1, size(y, 2)
            log_density(i) = log_density(i) + poisson_logpmf(y(i,j), lambda_value(j))
         end do
      end do
   end subroutine mvpois_log_density

   pure subroutine fit_mvcombi_component(y, weights, binary, center, variance, df, info)
      real(dp), intent(in) :: y(:,:) !! Mixed binary/Gaussian observations, shape `(n, d)`.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      logical, intent(in) :: binary(:) !! True for Bernoulli columns and false for Gaussian columns.
      real(dp), intent(out) :: center(:) !! Weighted centers/probabilities for every column, size `d`.
      real(dp), intent(out) :: variance(:) !! Gaussian variances for non-binary columns, packed in column order.
      integer, intent(out) :: df !! Upstream component parameter count `d + number_of_continuous_columns`.
      integer, intent(out) :: info !! Zero on success; nonzero if weighted covariance cannot be estimated.
      real(dp), allocatable :: covariance(:,:)
      logical :: ok
      integer :: j, m
      allocate(covariance(size(y,2), size(y,2)))
      call weighted_covariance(y, weights, center, covariance, ok)
      if (.not. ok) then
         info = 1
         return
      end if
      m = 0
      do j = 1, size(binary)
         if (binary(j)) then
            center(j) = clip_probability(center(j))
         else
            m = m + 1
            variance(m) = max(1.0e-12_dp, covariance(j,j))
         end if
      end do
      df = size(y, 2) + m
      info = 0
   end subroutine fit_mvcombi_component

   pure subroutine mvcombi_log_density(y, binary, center, variance, log_density)
      real(dp), intent(in) :: y(:,:) !! Mixed binary/Gaussian observations, shape `(n, d)`.
      logical, intent(in) :: binary(:) !! True for Bernoulli columns and false for Gaussian columns.
      real(dp), intent(in) :: center(:) !! Component centers/probabilities for all columns.
      real(dp), intent(in) :: variance(:) !! Packed variances for continuous columns.
      real(dp), intent(out) :: log_density(:) !! Per-row mixed independent log densities.
      real(dp) :: p
      integer :: i, j, m
      do i = 1, size(y, 1)
         log_density(i) = 0.0_dp
         m = 0
         do j = 1, size(y, 2)
            if (binary(j)) then
               p = clip_probability(center(j))
               log_density(i) = log_density(i) + y(i,j) * log(p) + (1.0_dp - y(i,j)) * log(1.0_dp - p)
            else
               m = m + 1
               log_density(i) = log_density(i) + normal_logpdf(y(i,j), center(j), sqrt(max(variance(m), 1.0e-12_dp)))
            end if
         end do
      end do
   end subroutine mvcombi_log_density

   pure subroutine fit_lognormal_component(y, weights, meanlog, sdlog, df, info)
      real(dp), intent(in) :: y(:) !! Strictly positive observations for log-normal fitting.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: meanlog !! Weighted mean of logarithms.
      real(dp), intent(out) :: sdlog !! Weighted maximum-likelihood standard deviation of logarithms.
      integer, intent(out) :: df !! Parameter count, equal to two.
      integer, intent(out) :: info !! Zero on success; nonzero for nonpositive data or zero total weight.
      real(dp), allocatable :: logy(:)
      real(dp) :: sw
      df = 2
      if (any(y <= 0.0_dp)) then
         info = 1
         meanlog = 0.0_dp
         sdlog = 1.0_dp
         return
      end if
      sw = sum(weights)
      if (sw <= tiny(1.0_dp)) then
         info = 2
         meanlog = 0.0_dp
         sdlog = 1.0_dp
         return
      end if
      allocate(logy(size(y)))
      logy = log(y)
      meanlog = dot_product(weights, logy) / sw
      sdlog = sqrt(max(1.0e-12_dp, dot_product(weights, (logy - meanlog)**2) / sw))
      info = 0
   end subroutine fit_lognormal_component

   pure subroutine lognormal_log_density(y, meanlog, sdlog, log_density)
      real(dp), intent(in) :: y(:) !! Positive observations for log-normal density evaluation.
      real(dp), intent(in) :: meanlog !! Mean on the logarithmic scale.
      real(dp), intent(in) :: sdlog !! Positive standard deviation on the logarithmic scale.
      real(dp), intent(out) :: log_density(:) !! Per-observation log-normal log densities.
      integer :: i
      do i = 1, size(y)
         if (y(i) <= 0.0_dp) then
            log_density(i) = -huge(1.0_dp)
         else
            log_density(i) = normal_logpdf(log(y(i)), meanlog, sdlog) - log(y(i))
         end if
      end do
   end subroutine lognormal_log_density

   pure subroutine fit_exponential_component(y, weights, rate, df, info)
      real(dp), intent(in) :: y(:) !! Nonnegative observations for exponential fitting.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: rate !! Fitted positive exponential rate.
      integer, intent(out) :: df !! Parameter count, equal to one.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid data or zero weighted mean.
      real(dp) :: sw, sy
      df = 1
      sw = sum(weights)
      sy = dot_product(weights, y)
      if (sw <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp) .or. any(y < 0.0_dp)) then
         rate = 1.0_dp
         info = 1
         return
      end if
      rate = sw / sy
      info = 0
   end subroutine fit_exponential_component

   pure subroutine exponential_log_density(y, rate, log_density)
      real(dp), intent(in) :: y(:) !! Nonnegative observations for exponential density evaluation.
      real(dp), intent(in) :: rate !! Positive exponential rate.
      real(dp), intent(out) :: log_density(:) !! Per-observation exponential log densities.
      integer :: i
      do i = 1, size(y)
         if (y(i) < 0.0_dp .or. rate <= 0.0_dp) then
            log_density(i) = -huge(1.0_dp)
         else
            log_density(i) = log(rate) - rate * y(i)
         end if
      end do
   end subroutine exponential_log_density

   pure subroutine fit_inverse_gaussian_component(y, weights, nu, lambda_value, df, info)
      real(dp), intent(in) :: y(:) !! Strictly positive observations for inverse-Gaussian fitting.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: nu !! Fitted inverse-Gaussian mean parameter.
      real(dp), intent(out) :: lambda_value !! Fitted positive inverse-Gaussian shape parameter.
      integer, intent(out) :: df !! Parameter count, equal to two.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid data or degenerate weighted moments.
      real(dp) :: sw, denom
      df = 2
      sw = sum(weights)
      if (sw <= tiny(1.0_dp) .or. any(y <= 0.0_dp)) then
         nu = 1.0_dp
         lambda_value = 1.0_dp
         info = 1
         return
      end if
      nu = dot_product(weights, y) / sw
      denom = dot_product(weights, 1.0_dp / y - 1.0_dp / nu)
      if (denom <= tiny(1.0_dp)) then
         lambda_value = huge(1.0_dp)
      else
         lambda_value = sw / denom
      end if
      info = 0
   end subroutine fit_inverse_gaussian_component

   pure subroutine inverse_gaussian_log_density(y, nu, lambda_value, log_density)
      real(dp), intent(in) :: y(:) !! Strictly positive observations for inverse-Gaussian density evaluation.
      real(dp), intent(in) :: nu !! Positive inverse-Gaussian mean parameter.
      real(dp), intent(in) :: lambda_value !! Positive inverse-Gaussian shape parameter.
      real(dp), intent(out) :: log_density(:) !! Per-observation inverse-Gaussian log densities.
      real(dp), parameter :: pi_local = acos(-1.0_dp)
      integer :: i
      do i = 1, size(y)
         if (y(i) <= 0.0_dp .or. nu <= 0.0_dp .or. lambda_value <= 0.0_dp) then
            log_density(i) = -huge(1.0_dp)
         else
            log_density(i) = 0.5_dp * (log(lambda_value) - log(2.0_dp * pi_local) - 3.0_dp * log(y(i))) - &
                             lambda_value * (y(i) - nu)**2 / (2.0_dp * nu * nu * y(i))
         end if
      end do
   end subroutine inverse_gaussian_log_density

   pure subroutine fit_gamma_component(y, weights, shape, rate, df, info)
      use flexmix_numeric, only : digamma_approx, trigamma_approx
      real(dp), intent(in) :: y(:) !! Strictly positive observations for gamma fitting.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: shape !! Fitted positive gamma shape.
      real(dp), intent(out) :: rate !! Fitted positive gamma rate.
      integer, intent(out) :: df !! Parameter count, equal to two.
      integer, intent(out) :: info !! Zero on convergence; nonzero for invalid data or a failed shape update.
      real(dp) :: sw, mean_y, mean_log_y, s, f, fp, shape_new, variance
      integer :: iter
      df = 2
      sw = sum(weights)
      if (sw <= tiny(1.0_dp) .or. any(y <= 0.0_dp)) then
         shape = 1.0_dp
         rate = 1.0_dp
         info = 1
         return
      end if
      mean_y = dot_product(weights, y) / sw
      mean_log_y = dot_product(weights, log(y)) / sw
      variance = dot_product(weights, (y - mean_y)**2) / sw
      shape = max(0.1_dp, mean_y * mean_y / max(variance, 1.0e-12_dp))
      s = log(mean_y) - mean_log_y
      do iter = 1, 100
         f = log(shape) - digamma_approx(shape) - s
         fp = 1.0_dp / shape - trigamma_approx(shape)
         if (abs(fp) <= 1.0e-14_dp) exit
         shape_new = shape - f / fp
         if (shape_new <= 0.0_dp) shape_new = 0.5_dp * shape
         if (abs(shape_new - shape) <= 1.0e-10_dp * (1.0_dp + shape)) then
            shape = shape_new
            exit
         end if
         shape = shape_new
      end do
      rate = shape / mean_y
      info = 0
   end subroutine fit_gamma_component

   pure subroutine gamma_log_density(y, shape, rate, log_density)
      real(dp), intent(in) :: y(:) !! Positive observations for gamma density evaluation.
      real(dp), intent(in) :: shape !! Positive gamma shape.
      real(dp), intent(in) :: rate !! Positive gamma rate.
      real(dp), intent(out) :: log_density(:) !! Per-observation gamma log densities.
      integer :: i
      do i = 1, size(y)
         if (y(i) <= 0.0_dp .or. shape <= 0.0_dp .or. rate <= 0.0_dp) then
            log_density(i) = -huge(1.0_dp)
         else
            log_density(i) = shape * log(rate) - log_gamma(shape) + (shape - 1.0_dp) * log(y(i)) - rate * y(i)
         end if
      end do
   end subroutine gamma_log_density

   pure subroutine fit_weibull_component(y, weights, shape, scale, df, info)
      real(dp), intent(in) :: y(:) !! Strictly positive observations for Weibull fitting.
      real(dp), intent(in) :: weights(:) !! Nonnegative posterior/case weights, size `n`.
      real(dp), intent(out) :: shape !! Fitted positive Weibull shape.
      real(dp), intent(out) :: scale !! Fitted positive Weibull scale.
      integer, intent(out) :: df !! Parameter count, equal to two.
      integer, intent(out) :: info !! Zero on convergence; nonzero for invalid data.
      real(dp), allocatable :: logy(:), yk(:)
      real(dp) :: sw, mean_log, variance_log, a, b, f, fp, shape_new
      integer :: iter
      df = 2
      sw = sum(weights)
      if (sw <= tiny(1.0_dp) .or. any(y <= 0.0_dp)) then
         shape = 1.0_dp
         scale = 1.0_dp
         info = 1
         return
      end if
      allocate(logy(size(y)), yk(size(y)))
      logy = log(y)
      mean_log = dot_product(weights, logy) / sw
      variance_log = dot_product(weights, (logy - mean_log)**2) / sw
      shape = max(0.1_dp, 1.2_dp / sqrt(max(variance_log, 1.0e-12_dp)))
      do iter = 1, 100
         yk = exp(min(700.0_dp, shape * logy))
         a = dot_product(weights, yk)
         b = dot_product(weights, yk * logy)
         f = 1.0_dp / shape + mean_log - b / a
         fp = -1.0_dp / (shape * shape) - &
              (dot_product(weights, yk * logy * logy) * a - b * b) / (a * a)
         if (abs(fp) <= 1.0e-14_dp) exit
         shape_new = shape - f / fp
         if (shape_new <= 0.0_dp) shape_new = 0.5_dp * shape
         if (abs(shape_new - shape) <= 1.0e-10_dp * (1.0_dp + shape)) then
            shape = shape_new
            exit
         end if
         shape = shape_new
      end do
      yk = exp(min(700.0_dp, shape * logy))
      scale = (dot_product(weights, yk) / sw)**(1.0_dp / shape)
      info = 0
   end subroutine fit_weibull_component

   pure subroutine weibull_log_density(y, shape, scale, log_density)
      real(dp), intent(in) :: y(:) !! Positive observations for Weibull density evaluation.
      real(dp), intent(in) :: shape !! Positive Weibull shape.
      real(dp), intent(in) :: scale !! Positive Weibull scale.
      real(dp), intent(out) :: log_density(:) !! Per-observation Weibull log densities.
      integer :: i
      do i = 1, size(y)
         if (y(i) <= 0.0_dp .or. shape <= 0.0_dp .or. scale <= 0.0_dp) then
            log_density(i) = -huge(1.0_dp)
         else
            log_density(i) = log(shape) - log(scale) + (shape - 1.0_dp) * (log(y(i)) - log(scale)) - &
                             (y(i) / scale)**shape
         end if
      end do
   end subroutine weibull_log_density

end module flexmix_components
