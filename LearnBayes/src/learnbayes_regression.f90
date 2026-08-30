module learnbayes_regression
   use learnbayes_distributions, only: dmnorm, rigamma, rmnorm, rtruncated_normal
   use learnbayes_kinds, only: dp
   use learnbayes_linalg, only: cholesky_lower, inverse_matrix, least_squares, solve_linear
   use learnbayes_models, only: reg_gprior_post
   use learnbayes_rng, only: rng_normal, rng_state
   use learnbayes_sampling, only: laplace
   use learnbayes_special, only: normal_cdf
   use learnbayes_types, only: blinreg_result, laplace_result, log_density_callback, model_selection_result, probit_result
   implicit none
   private

   public :: bayes_model_selection
   public :: bayes_probit
   public :: bayesresiduals
   public :: blinreg
   public :: blinregexpected
   public :: blinregpred
   public :: bprobit_probs

contains

   subroutine blinreg(rng, y, x, m, result, info, c0, beta0)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for inverse-gamma and regression-coefficient draws.
      real(dp), intent(in) :: y(:) !! Gaussian regression response vector.
      real(dp), intent(in) :: x(:, :) !! Design matrix with observations in rows and regression terms in columns.
      integer, intent(in) :: m !! Number of posterior draws to generate.
      type(blinreg_result), intent(out) :: result !! Posterior draws of regression coefficients and residual standard deviation.
      integer, intent(out) :: info !! Zero on success; nonzero for singular designs or inconsistent prior arguments.
      real(dp), intent(in), optional :: c0 !! Optional positive Zellner g-prior multiplier.
      real(dp), intent(in), optional :: beta0(:) !! Optional prior mean vector paired with c0; length must equal ncol(x).
      real(dp), allocatable :: bhat(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: xtx_inv(:, :)
      real(dp), allocatable :: xtx(:, :)
      real(dp), allocatable :: vbeta(:, :)
      real(dp), allocatable :: zero(:)
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: beta_mean(:)
      real(dp), allocatable :: delta(:)
      real(dp) :: s2
      real(dp) :: shape
      real(dp) :: rate
      real(dp) :: sigma
      integer :: i
      integer :: n
      integer :: p
      integer :: rank

      n = size(y)
      p = size(x, 2)
      info = 0
      allocate(bhat(p), residuals(n), xtx_inv(p, p), xtx(p, p), vbeta(p, p), zero(p), z(m, p), beta_mean(p))
      call least_squares(x, y, bhat, residuals, s2, xtx_inv, rank, info)
      if (info /= 0) return
      xtx = matmul(transpose(x), x)
      if (present(c0) .neqv. present(beta0)) then
         info = -3
         return
      end if
      if (.not. present(c0)) then
         shape = 0.5_dp*real(n - p, dp)
         rate = 0.5_dp*dot_product(residuals, residuals)
         beta_mean = bhat
         vbeta = xtx_inv
      else
         if (c0 <= 0.0_dp .or. size(beta0) /= p) then
            info = -4
            return
         end if
         allocate(delta(p))
         delta = beta0 - bhat
         shape = 0.5_dp*real(n, dp)
         rate = 0.5_dp*dot_product(residuals, residuals) + &
            0.5_dp*dot_product(delta, matmul(xtx, delta))/(c0 + 1.0_dp)
         beta_mean = (beta0 + c0*bhat)/(c0 + 1.0_dp)
         vbeta = xtx_inv*c0/(c0 + 1.0_dp)
      end if
      allocate(result%beta(m, p), result%sigma(m))
      zero = 0.0_dp
      call rmnorm(rng, m, zero, vbeta, z, info)
      if (info /= 0) return
      do i = 1, m
         sigma = sqrt(rigamma(rng, shape, rate))
         result%sigma(i) = sigma
         result%beta(i, :) = beta_mean + sigma*z(i, :)
      end do
   end subroutine blinreg

   subroutine blinregexpected(x1, theta_sample, expected)
      real(dp), intent(in) :: x1(:, :) !! New design matrix with the same number of columns as posterior beta draws.
      type(blinreg_result), intent(in) :: theta_sample !! Posterior regression draws returned by blinreg.
      real(dp), intent(out) :: expected(:, :) !! Expected-response draws shaped (number of posterior draws, nrow(x1)).
      integer :: j

      do j = 1, size(x1, 1)
         expected(:, j) = matmul(theta_sample%beta, x1(j, :))
      end do
   end subroutine blinregexpected

   subroutine blinregpred(rng, x1, theta_sample, prediction)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for posterior predictive Gaussian errors.
      real(dp), intent(in) :: x1(:, :) !! New design matrix with the same number of columns as posterior beta draws.
      type(blinreg_result), intent(in) :: theta_sample !! Posterior regression draws returned by blinreg.
      real(dp), intent(out) :: prediction(:, :) !! Posterior predictive draws shaped (number of posterior draws, nrow(x1)).
      integer :: i
      integer :: j

      do j = 1, size(x1, 1)
         prediction(:, j) = matmul(theta_sample%beta, x1(j, :))
         do i = 1, size(theta_sample%sigma)
            prediction(i, j) = prediction(i, j) + rng_normal(rng)*theta_sample%sigma(i)
         end do
      end do
   end subroutine blinregpred

   subroutine bayesresiduals(y, x, post, cutoff, probability, info)
      real(dp), intent(in) :: y(:) !! Regression response vector used to compute ordinary least-squares residuals.
      real(dp), intent(in) :: x(:, :) !! Design matrix used both for OLS residuals and leverage values.
      type(blinreg_result), intent(in) :: post !! Posterior regression simulation containing residual-standard-deviation draws.
      real(dp), intent(in) :: cutoff !! Positive standardized-residual cutoff k defining an outlier.
      real(dp), intent(out) :: probability(:) !! Posterior probability that each Bayesian residual exceeds the cutoff in magnitude.
      integer, intent(out) :: info !! Zero on success; nonzero if the design is singular or dimensions are inconsistent.
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: xtx_inv(:, :)
      real(dp) :: sigma2
      real(dp) :: leverage
      real(dp) :: z1
      real(dp) :: z2
      integer :: i
      integer :: j
      integer :: rank

      allocate(beta(size(x, 2)), residuals(size(y)), xtx_inv(size(x, 2), size(x, 2)))
      call least_squares(x, y, beta, residuals, sigma2, xtx_inv, rank, info)
      if (info /= 0) return
      probability = 0.0_dp
      do i = 1, size(y)
         leverage = dot_product(x(i, :), matmul(xtx_inv, x(i, :)))
         leverage = max(leverage, epsilon(1.0_dp))
         do j = 1, size(post%sigma)
            z1 = (cutoff - residuals(i)/post%sigma(j))/sqrt(leverage)
            z2 = (-cutoff - residuals(i)/post%sigma(j))/sqrt(leverage)
            probability(i) = probability(i) + 1.0_dp - normal_cdf(z1) + normal_cdf(z2)
         end do
         probability(i) = probability(i)/real(size(post%sigma), dp)
      end do
   end subroutine bayesresiduals

   subroutine probit_mle(y, x, beta, fitted, info)
      integer, intent(in) :: y(:) !! Binary response vector containing only zeros and ones.
      real(dp), intent(in) :: x(:, :) !! Probit-regression design matrix.
      real(dp), intent(out) :: beta(:) !! Maximum-likelihood probit coefficient estimate.
      real(dp), intent(out) :: fitted(:) !! Fitted Bernoulli probabilities at beta.
      integer, intent(out) :: info !! Zero on convergence; nonzero for dimensions, singular systems, or iteration failure.
      real(dp), allocatable :: score(:)
      real(dp), allocatable :: fisher(:, :)
      real(dp), allocatable :: step(:)
      real(dp) :: eta
      real(dp) :: p
      real(dp) :: phi
      real(dp) :: w
      integer :: i
      integer :: iter
      integer :: pcols

      pcols = size(x, 2)
      info = 0
      beta = 0.0_dp
      allocate(score(pcols), fisher(pcols, pcols), step(pcols))
      do iter = 1, 100
         score = 0.0_dp
         fisher = 0.0_dp
         do i = 1, size(y)
            eta = dot_product(x(i, :), beta)
            p = min(1.0_dp - 1.0e-12_dp, max(1.0e-12_dp, normal_cdf(eta)))
            phi = exp(-0.5_dp*eta*eta)/sqrt(2.0_dp*acos(-1.0_dp))
            score = score + x(i, :)*real(y(i) - p, dp)*phi/(p*(1.0_dp - p))
            w = phi*phi/(p*(1.0_dp - p))
            fisher = fisher + w*outer_product(x(i, :), x(i, :))
         end do
         call solve_linear(fisher, score, step, info)
         if (info /= 0) return
         beta = beta + step
         if (maxval(abs(step)) < 1.0e-9_dp) exit
      end do
      if (iter > 100) info = 2
      do i = 1, size(y)
         fitted(i) = normal_cdf(dot_product(x(i, :), beta))
      end do
   end subroutine probit_mle

   pure function outer_product(a, b) result(c)
      real(dp), intent(in) :: a(:) !! Left vector of the outer product.
      real(dp), intent(in) :: b(:) !! Right vector of the outer product.
      real(dp) :: c(size(a), size(b))
      integer :: i
      integer :: j

      do j = 1, size(b)
         do i = 1, size(a)
            c(i, j) = a(i)*b(j)
         end do
      end do
   end function outer_product

   subroutine bayes_probit(rng, y, x, m, result, info, prior_beta, prior_precision)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for latent-normal and coefficient draws.
      integer, intent(in) :: y(:) !! Binary response vector containing zeros and ones.
      real(dp), intent(in) :: x(:, :) !! Probit-regression design matrix with observations in rows.
      integer, intent(in) :: m !! Number of Gibbs iterations to simulate.
      type(probit_result), intent(out) :: result !! Posterior beta draws and optional Chib-style log marginal likelihood estimate.
      integer, intent(out) :: info !! Zero on success; nonzero for dimensions or singular matrix calculations.
      real(dp), intent(in), optional :: prior_beta(:) !! Optional normal prior mean vector; defaults to zero.
      real(dp), intent(in), optional :: prior_precision(:, :) !! Optional prior precision matrix; defaults to zero.
      real(dp), allocatable :: beta_s(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: beta0(:)
      real(dp), allocatable :: bi(:, :)
      real(dp), allocatable :: postvar(:, :)
      real(dp), allocatable :: postchol(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: mn(:)
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: noise(:)
      real(dp), allocatable :: priorcov(:, :)
      real(dp), allocatable :: system(:, :)
      real(dp) :: post_ord
      real(dp) :: logf
      real(dp) :: logg
      real(dp) :: lo
      real(dp) :: hi
      integer :: i
      integer :: j
      integer :: n
      integer :: p
      logical :: informative

      n = size(y)
      p = size(x, 2)
      info = 0
      allocate(beta_s(p), fitted(n), beta(p), beta0(p), bi(p, p), postvar(p, p), postchol(p, p), &
         rhs(p), mn(p), z(n), noise(p), system(p, p), result%beta(m, p))
      call probit_mle(y, x, beta_s, fitted, info)
      if (info /= 0) return
      beta0 = 0.0_dp
      if (present(prior_beta)) beta0 = prior_beta
      bi = 0.0_dp
      if (present(prior_precision)) bi = prior_precision
      system = bi + matmul(transpose(x), x)
      call inverse_matrix(system, postvar, info)
      if (info /= 0) return
      call cholesky_lower(postvar, postchol, info)
      if (info /= 0) return
      beta = beta_s
      post_ord = 0.0_dp
      do i = 1, m
         do j = 1, n
            if (y(j) == 0) then
               lo = -huge(1.0_dp)
               hi = 0.0_dp
            else
               lo = 0.0_dp
               hi = huge(1.0_dp)
            end if
            z(j) = rtruncated_normal(rng, lo, hi, dot_product(x(j, :), beta), 1.0_dp)
         end do
         rhs = matmul(bi, beta0) + matmul(transpose(x), z)
         call solve_linear(system, rhs, mn, info)
         if (info /= 0) return
         do j = 1, p
            noise(j) = rng_normal(rng)
         end do
         beta = mn + matmul(postchol, noise)
         post_ord = post_ord + dmnorm(beta_s, mn, postvar)
         result%beta(i, :) = beta
      end do
      informative = sum(abs(bi)) > 0.0_dp
      result%has_log_marginal = informative
      if (informative) then
         logf = 0.0_dp
         do i = 1, n
            logf = logf + real(y(i), dp)*log(max(fitted(i), tiny(1.0_dp))) + &
               real(1 - y(i), dp)*log(max(1.0_dp - fitted(i), tiny(1.0_dp)))
         end do
         allocate(priorcov(p, p))
         call inverse_matrix(bi, priorcov, info)
         if (info /= 0) return
         logg = dmnorm(beta_s, beta0, priorcov, .true.)
         result%log_marginal = logf + logg - log(post_ord/real(m, dp))
      end if
   end subroutine bayes_probit

   subroutine bprobit_probs(x1, beta_draws, probability)
      real(dp), intent(in) :: x1(:, :) !! New probit-regression design matrix.
      real(dp), intent(in) :: beta_draws(:, :) !! Posterior beta draws with one draw per row.
      real(dp), intent(out) :: probability(:, :) !! Probit probabilities shaped (number of draws, nrow(x1)).
      integer :: i
      integer :: j

      do j = 1, size(x1, 1)
         do i = 1, size(beta_draws, 1)
            probability(i, j) = normal_cdf(dot_product(x1(j, :), beta_draws(i, :)))
         end do
      end do
   end subroutine bprobit_probs

   subroutine bayes_model_selection(y, x_input, c, result, constant_in_x, info)
      real(dp), intent(in) :: y(:) !! Gaussian response vector.
      real(dp), intent(in) :: x_input(:, :) !! Covariate matrix; optionally already includes a leading intercept column.
      real(dp), intent(in) :: c !! Positive Zellner g-prior multiplier used for every candidate model.
      type(model_selection_result), intent(out) :: result !! Candidate models, log marginals, probabilities, and flags.
      logical, intent(in), optional :: constant_in_x !! True if x_input already contains a leading constant; defaults to true.
      integer, intent(out) :: info !! Zero on success; nonzero if a candidate model is singular or dimensions are invalid.
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: x0(:, :)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: xtx_inv(:, :)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: packed(:, :)
      real(dp) :: sigma2
      real(dp) :: mx
      integer :: rank
      integer :: nmodels
      integer :: p
      integer :: j
      integer :: k
      integer :: nsel
      logical :: has_constant
      type(log_density_callback) :: cb
      type(laplace_result) :: fit

      has_constant = .true.
      if (present(constant_in_x)) has_constant = constant_in_x
      if (has_constant) then
         allocate(x(size(x_input, 1), size(x_input, 2)))
         x = x_input
      else
         allocate(x(size(x_input, 1), size(x_input, 2) + 1))
         x(:, 1) = 1.0_dp
         x(:, 2:) = x_input
      end if
      p = size(x, 2) - 1
      nmodels = 2**p
      allocate(result%included(nmodels, p + 1), result%log_marginal(nmodels), result%probability(nmodels), &
         result%converged(nmodels))
      result%included = .false.
      result%included(:, 1) = .true.
      info = 0
      do j = 0, nmodels - 1
         do k = 1, p
            result%included(j + 1, k + 1) = btest(j, k - 1)
         end do
         nsel = count(result%included(j + 1, :))
         allocate(x0(size(x, 1), nsel), beta(nsel), residuals(size(y)), xtx_inv(nsel, nsel), theta(nsel + 1))
         x0 = pack_columns(x, result%included(j + 1, :))
         call least_squares(x0, y, beta, residuals, sigma2, xtx_inv, rank, info)
         if (info /= 0) return
         theta(1:nsel) = beta
         theta(nsel + 1) = log(sqrt(sigma2))
         allocate(packed(size(y), nsel + 1))
         packed(:, 1) = y
         packed(:, 2:) = x0
         cb%eval => regression_model_callback
         cb%data = packed
         cb%params = [c]
         call laplace(cb, theta, fit)
         result%log_marginal(j + 1) = fit%log_integral
         result%converged(j + 1) = fit%converged
         deallocate(x0, beta, residuals, xtx_inv, theta, packed)
      end do
      mx = maxval(result%log_marginal)
      result%probability = exp(result%log_marginal - mx)
      result%probability = result%probability/sum(result%probability)
   end subroutine bayes_model_selection

   function regression_model_callback(theta, data, params) result(value)
      real(dp), intent(in) :: theta(:) !! Candidate regression coefficients followed by log residual standard deviation.
      real(dp), intent(in) :: data(:, :) !! Packed matrix with response in column 1 and candidate design matrix thereafter.
      real(dp), intent(in) :: params(:) !! Numeric constants with Zellner g-prior multiplier c in the first element.
      real(dp) :: value
      real(dp), allocatable :: beta0(:)

      allocate(beta0(size(theta) - 1))
      beta0 = 0.0_dp
      value = reg_gprior_post(theta, data(:, 1), data(:, 2:), params(1), beta0)
   end function regression_model_callback

   function pack_columns(x, keep) result(out)
      real(dp), intent(in) :: x(:, :) !! Source matrix whose selected columns are copied in order.
      logical, intent(in) :: keep(:) !! Logical selector with one element per source column.
      real(dp) :: out(size(x, 1), count(keep))
      integer :: j
      integer :: k

      k = 0
      do j = 1, size(x, 2)
         if (keep(j)) then
            k = k + 1
            out(:, k) = x(:, j)
         end if
      end do
   end function pack_columns

end module learnbayes_regression
