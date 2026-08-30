module learnbayes_hierarchical
   use learnbayes_distributions, only: rdirichlet, rigamma, rmnorm, rtruncated_normal
   use learnbayes_kinds, only: dp
   use learnbayes_linalg, only: cholesky_lower, inverse_matrix, solve_linear
   use learnbayes_rng, only: rng_gamma, rng_normal, rng_state
   use learnbayes_special, only: beta_log, quantile_type7
   implicit none
   private

   public :: bayes_influence
   public :: bfindep
   public :: hiergibbs
   public :: normpostpred
   public :: normpostsim
   public :: ordergibbs
   public :: robustt

contains

   subroutine normpostsim(rng, data, m, mu_draw, sigma2_draw, prior_sigma2, prior_mu)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for inverse-gamma and Gaussian posterior draws.
      real(dp), intent(in) :: data(:) !! Independent normal observations.
      integer, intent(in) :: m !! Number of posterior draws requested.
      real(dp), intent(out) :: mu_draw(:) !! Posterior draws of the normal mean, with length m.
      real(dp), intent(out) :: sigma2_draw(:) !! Posterior draws of the normal variance, with length m.
      real(dp), intent(in), optional :: prior_sigma2(2) !! Optional inverse-gamma shape and rate for sigma squared.
      real(dp), intent(in), optional :: prior_mu(2) !! Optional normal prior mean and variance for mu.
      real(dp) :: s
      real(dp) :: xbar
      real(dp) :: sigma2
      real(dp) :: precision
      real(dp) :: post_mean
      real(dp) :: post_var
      real(dp) :: a1
      real(dp) :: b1
      real(dp) :: mu
      integer :: j
      integer :: n

      n = size(data)
      xbar = sum(data)/real(n, dp)
      s = sum((data - xbar)**2)
      if (.not. present(prior_sigma2) .and. .not. present(prior_mu)) then
         do j = 1, m
            sigma2_draw(j) = s/(2.0_dp*rng_gamma(rng, 0.5_dp*real(n - 1, dp), 1.0_dp))
            mu_draw(j) = rng_normal(rng, xbar, sqrt(sigma2_draw(j)/real(n, dp)))
         end do
         return
      end if
      if (.not. present(prior_sigma2) .or. .not. present(prior_mu)) then
         mu_draw = 0.0_dp
         sigma2_draw = 0.0_dp
         return
      end if
      sigma2 = s/real(n, dp)
      do j = 1, m
         precision = real(n, dp)/sigma2 + 1.0_dp/prior_mu(2)
         post_mean = (xbar*real(n, dp)/sigma2 + prior_mu(1)/prior_mu(2))/precision
         post_var = 1.0_dp/precision
         mu = rng_normal(rng, post_mean, sqrt(post_var))
         a1 = prior_sigma2(1) + 0.5_dp*real(n, dp)
         b1 = prior_sigma2(2) + 0.5_dp*sum((data - mu)**2)
         sigma2 = rigamma(rng, a1, b1)
         sigma2_draw(j) = sigma2
         mu_draw(j) = mu
      end do
   end subroutine normpostsim

   subroutine normpostpred(rng, mu, sigma2, sample_size, statistic, value)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for posterior predictive normal samples.
      real(dp), intent(in) :: mu(:) !! Posterior draws of the normal mean.
      real(dp), intent(in) :: sigma2(:) !! Posterior draws of the positive normal variance, paired with mu.
      integer, intent(in) :: sample_size !! Number of observations in each replicated posterior predictive dataset.
      character(len=*), intent(in), optional :: statistic !! Statistic name: min, max, mean, or sd; defaults to min.
      real(dp), intent(out) :: value(:) !! One posterior predictive statistic for each posterior parameter draw.
      real(dp), allocatable :: sample(:)
      character(len=16) :: stat
      real(dp) :: mean_value
      integer :: i
      integer :: j

      stat = 'min'
      if (present(statistic)) stat = adjustl(statistic)
      allocate(sample(sample_size))
      do i = 1, size(mu)
         do j = 1, sample_size
            sample(j) = rng_normal(rng, mu(i), sqrt(sigma2(i)))
         end do
         select case (trim(stat))
         case ('max')
            value(i) = maxval(sample)
         case ('mean')
            value(i) = sum(sample)/real(sample_size, dp)
         case ('sd')
            mean_value = sum(sample)/real(sample_size, dp)
            if (sample_size > 1) then
               value(i) = sqrt(sum((sample - mean_value)**2)/real(sample_size - 1, dp))
            else
               value(i) = 0.0_dp
            end if
         case default
            value(i) = minval(sample)
         end select
      end do
   end subroutine normpostpred

   subroutine robustt(rng, y, df, m, mu_draw, sigma2_draw, lambda_draw)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for all Student-t mixture Gibbs updates.
      real(dp), intent(in) :: y(:) !! Observations modeled with a Student-t location/scale distribution.
      real(dp), intent(in) :: df !! Positive Student-t degrees of freedom v.
      integer, intent(in) :: m !! Number of Gibbs iterations to simulate.
      real(dp), intent(out) :: mu_draw(:) !! Simulated location parameters, with length m.
      real(dp), intent(out) :: sigma2_draw(:) !! Simulated scale variances, with length m.
      real(dp), intent(out) :: lambda_draw(:, :) !! Latent gamma precisions shaped (m,size(y)).
      real(dp), allocatable :: lambda(:)
      real(dp) :: mu
      real(dp) :: sigma2
      real(dp) :: mean_y
      real(dp) :: sd_y
      integer :: i
      integer :: j
      integer :: n

      n = size(y)
      mean_y = sum(y)/real(n, dp)
      if (n > 1) then
         sd_y = sqrt(sum((y - mean_y)**2)/real(n - 1, dp))
      else
         sd_y = 1.0_dp
      end if
      mu = mean_y
      sigma2 = sd_y*sd_y
      allocate(lambda(n))
      lambda = 1.0_dp
      do i = 1, m
         do j = 1, n
            lambda(j) = rng_gamma(rng, 0.5_dp*(df + 1.0_dp), &
               0.5_dp*df + 0.5_dp*(y(j) - mu)**2/sigma2)
         end do
         mu = rng_normal(rng, sum(y*lambda)/sum(lambda), sqrt(sigma2/sum(lambda)))
         sigma2 = rigamma(rng, 0.5_dp*real(n, dp), 0.5_dp*sum(lambda*(y - mu)**2))
         mu_draw(i) = mu
         sigma2_draw(i) = sigma2
         lambda_draw(i, :) = lambda
      end do
   end subroutine robustt

   subroutine hiergibbs(rng, data, m, beta_draw, mu_draw, sigma2_pi_draw, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for all hierarchical Gibbs updates.
      real(dp), intent(in) :: data(:, :) !! Matrix with observed means, sample sizes, and two covariates in columns 1:4.
      integer, intent(in) :: m !! Number of Gibbs iterations to simulate.
      real(dp), intent(out) :: beta_draw(:, :) !! Regression-coefficient draws shaped (m,3).
      real(dp), intent(out) :: mu_draw(:, :) !! Cell-mean draws shaped (m,nrow(data)).
      real(dp), intent(out) :: sigma2_pi_draw(:) !! Second-stage variance draws, with length m.
      integer, intent(out) :: info !! Zero on success; nonzero if required covariance inversions fail.
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: nobs(:)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: s2(:)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: pvar(:, :)
      real(dp), allocatable :: precision(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: pmean(:)
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: postvar(:)
      real(dp), allocatable :: postmean(:)
      real(dp) :: b1(3)
      real(dp) :: bvar(3, 3)
      real(dp) :: ibvar(3, 3)
      real(dp) :: sigma2_pi
      real(dp), parameter :: prior_s = 0.02_dp
      real(dp), parameter :: prior_v = 16.0_dp
      integer :: i
      integer :: j
      integer :: n

      n = size(data, 1)
      allocate(y(n), nobs(n), x(n, 3), s2(n), mu(n), pvar(3, 3), precision(3, 3), rhs(3), pmean(3), &
         l(3, 3), z(3), postvar(n), postmean(n))
      y = data(:, 1)
      nobs = data(:, 2)
      x(:, 1) = 1.0_dp
      x(:, 2) = data(:, 3)
      x(:, 3) = data(:, 4)
      s2 = 0.65_dp**2/nobs
      b1 = [0.55_dp, 0.018_dp, 0.033_dp]
      bvar = reshape([ &
         8.49e-03_dp, -1.94e-05_dp, -2.88e-04_dp, &
         -1.94e-05_dp, 7.34e-07_dp, -1.52e-06_dp, &
         -2.88e-04_dp, -1.52e-06_dp, 1.71e-05_dp], [3, 3])
      call inverse_matrix(bvar, ibvar, info)
      if (info /= 0) return
      mu = y
      sigma2_pi = 0.006_dp
      do j = 1, m
         precision = ibvar + matmul(transpose(x), x)/sigma2_pi
         call inverse_matrix(precision, pvar, info)
         if (info /= 0) return
         rhs = matmul(ibvar, b1) + matmul(transpose(x), mu)/sigma2_pi
         pmean = matmul(pvar, rhs)
         call cholesky_lower(pvar, l, info)
         if (info /= 0) return
         do i = 1, 3
            z(i) = rng_normal(rng)
         end do
         beta_draw(j, :) = pmean + matmul(l, z)
         sigma2_pi = (0.5_dp*sum((mu - matmul(x, beta_draw(j, :)))**2) + 0.5_dp*prior_s)/ &
            rng_gamma(rng, 0.5_dp*(real(n, dp) + prior_v), 1.0_dp)
         postvar = 1.0_dp/(1.0_dp/s2 + 1.0_dp/sigma2_pi)
         postmean = (y/s2 + matmul(x, beta_draw(j, :))/sigma2_pi)*postvar
         do i = 1, n
            mu(i) = rng_normal(rng, postmean(i), sqrt(postvar(i)))
         end do
         mu_draw(j, :) = mu
         sigma2_pi_draw(j) = sigma2_pi
      end do
   end subroutine hiergibbs

   subroutine ordergibbs(rng, data, m, mu_draw, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for truncated-normal order-restricted Gibbs updates.
      real(dp), intent(in) :: data(:, :) !! Forty-row matrix of sample means and sample sizes from the order example.
      integer, intent(in) :: m !! Number of Gibbs sweeps to simulate.
      real(dp), intent(out) :: mu_draw(:, :) !! Simulated ordered means shaped (m,40), flattened in Fortran array order.
      integer, intent(out) :: info !! Zero on success; nonzero when the specialized 40-by-2 input shape is not supplied.
      real(dp) :: y(8, 5)
      real(dp) :: nobs(8, 5)
      real(dp) :: mu(10, 7)
      real(dp) :: lo
      real(dp) :: hi
      real(dp) :: sd
      integer :: i
      integer :: j
      integer :: k
      integer :: ii
      integer :: jj
      real(dp), parameter :: s = 0.65_dp
      real(dp), parameter :: init(8, 5) = reshape([ &
         1.59_dp, 1.85_dp, 1.85_dp, 2.04_dp, 2.31_dp, 2.37_dp, 2.37_dp, 2.64_dp, &
         1.59_dp, 1.85_dp, 1.85_dp, 2.11_dp, 2.33_dp, 2.47_dp, 2.63_dp, 3.02_dp, &
         1.59_dp, 1.85_dp, 1.85_dp, 2.11_dp, 2.33_dp, 2.64_dp, 2.74_dp, 3.02_dp, &
         1.67_dp, 1.88_dp, 2.10_dp, 2.33_dp, 2.33_dp, 2.66_dp, 2.76_dp, 3.07_dp, &
         1.88_dp, 1.88_dp, 2.10_dp, 2.33_dp, 2.33_dp, 2.66_dp, 2.91_dp, 3.34_dp], [8, 5])

      info = 0
      if (size(data, 1) /= 40 .or. size(data, 2) < 2) then
         info = 1
         return
      end if
      do j = 1, 5
         do i = 1, 8
            ii = 9 - i
            jj = (j - 1)*8 + i
            y(ii, j) = data(jj, 1)
            nobs(ii, j) = data(jj, 2)
         end do
      end do
      mu = huge(1.0_dp)
      mu(1, :) = -huge(1.0_dp)
      mu(:, 1) = -huge(1.0_dp)
      mu(1, 1) = huge(1.0_dp)
      mu(2:9, 2:6) = init
      do k = 1, m
         do i = 2, 9
            do j = 2, 6
               lo = max(mu(i - 1, j), mu(i, j - 1))
               hi = min(mu(i + 1, j), mu(i, j + 1))
               sd = s/sqrt(nobs(i - 1, j - 1))
               mu(i, j) = rtruncated_normal(rng, lo, hi, y(i - 1, j - 1), sd)
            end do
         end do
         mu_draw(k, :) = reshape(mu(2:9, 2:6), [40])
      end do
   end subroutine ordergibbs

   subroutine bfindep(rng, y, kappa, m, bf, nse)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for Dirichlet Monte Carlo draws.
      real(dp), intent(in) :: y(:, :) !! Nonnegative two-way contingency-table counts.
      real(dp), intent(in) :: kappa !! Positive Dirichlet precision parameter K.
      integer, intent(in) :: m !! Number of Monte Carlo iterations.
      real(dp), intent(out) :: bf !! Monte Carlo estimate of the Bayes factor against independence.
      real(dp), intent(out) :: nse !! Estimated Monte Carlo standard error of the Bayes-factor mean.
      real(dp), allocatable :: row_total(:)
      real(dp), allocatable :: col_total(:)
      real(dp), allocatable :: eta_a(:, :)
      real(dp), allocatable :: eta_b(:, :)
      real(dp), allocatable :: logint(:)
      real(dp), allocatable :: value(:)
      real(dp) :: mean_value
      real(dp) :: alpha
      integer :: i
      integer :: j
      integer :: r
      integer :: ncol
      integer :: col

      r = size(y, 1)
      ncol = size(y, 2)
      allocate(row_total(r), col_total(ncol), eta_a(m, r), eta_b(m, ncol), logint(m), value(m))
      row_total = sum(y, dim=2)
      col_total = sum(y, dim=1)
      call rdirichlet(rng, m, row_total + 1.0_dp, eta_a)
      call rdirichlet(rng, m, col_total + 1.0_dp, eta_b)
      logint = 0.0_dp
      do i = 1, m
         do j = 1, r
            logint(i) = logint(i) - row_total(j)*log(eta_a(i, j))
         end do
         do j = 1, ncol
            logint(i) = logint(i) - col_total(j)*log(eta_b(i, j))
         end do
         do j = 1, r
            do col = 1, ncol
               alpha = kappa*eta_a(i, j)*eta_b(i, col)
               logint(i) = logint(i) + log_gamma(alpha + y(j, col)) - log_gamma(alpha)
            end do
         end do
         logint(i) = logint(i) - (log_gamma(kappa + sum(y)) - log_gamma(kappa))
      end do
      value = exp(logint)
      bf = sum(value)/real(m, dp)
      mean_value = bf
      if (m > 1) then
         nse = sqrt(sum((value - mean_value)**2)/real(m - 1, dp))/sqrt(real(m, dp))
      else
         nse = 0.0_dp
      end if
   end subroutine bfindep

   subroutine bayes_influence(rng, theta, data, summary, summary_obs)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for weighted bootstrap resampling.
      real(dp), intent(in) :: theta(:, :) !! Posterior draws with logit eta in column 1 and log K in column 2.
      real(dp), intent(in) :: data(:, :) !! Binomial groups with successes in column 1 and trials in column 2.
      real(dp), intent(out) :: summary(3) !! Overall 5%, 50%, and 95% quantiles of log K.
      real(dp), intent(out) :: summary_obs(:, :) !! Leave-one-group influence quantiles shaped (nrow(data),3).
      real(dp), allocatable :: weight(:)
      real(dp), allocatable :: prob(:)
      real(dp), allocatable :: sample(:)
      real(dp), allocatable :: kappa(:)
      real(dp), allocatable :: eta(:)
      real(dp) :: y
      real(dp) :: n
      integer :: i
      integer :: j
      integer :: index
      integer :: m

      m = size(theta, 1)
      allocate(weight(m), prob(m), sample(m), kappa(m), eta(m))
      kappa = exp(theta(:, 2))
      eta = exp(theta(:, 1))/(1.0_dp + exp(theta(:, 1)))
      summary = [quantile_type7(theta(:, 2), 0.05_dp), quantile_type7(theta(:, 2), 0.5_dp), &
         quantile_type7(theta(:, 2), 0.95_dp)]
      do i = 1, size(data, 1)
         y = data(i, 1)
         n = data(i, 2)
         do j = 1, m
            weight(j) = exp(beta_log(kappa(j)*eta(j), kappa(j)*(1.0_dp - eta(j))) - &
               beta_log(kappa(j)*eta(j) + y, kappa(j)*(1.0_dp - eta(j)) + n - y))
         end do
         prob = weight/sum(weight)
         do j = 1, m
            index = discrete_draw_local(rng, prob)
            sample(j) = theta(index, 2)
         end do
         summary_obs(i, :) = [quantile_type7(sample, 0.05_dp), quantile_type7(sample, 0.5_dp), &
            quantile_type7(sample, 0.95_dp)]
      end do
   end subroutine bayes_influence

   function discrete_draw_local(rng, prob) result(index)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for one categorical draw.
      real(dp), intent(in) :: prob(:) !! Category probabilities or nonnegative weights.
      integer :: index
      real(dp) :: u
      real(dp) :: acc
      integer :: i

      u = rng_uniform_local(rng)*sum(prob)
      acc = 0.0_dp
      do i = 1, size(prob)
         acc = acc + prob(i)
         if (u <= acc) then
            index = i
            return
         end if
      end do
      index = size(prob)
   end function discrete_draw_local

   function rng_uniform_local(rng) result(value)
      use learnbayes_rng, only: rng_uniform
      type(rng_state), intent(inout) :: rng !! Mutable generator state advanced by one uniform draw.
      real(dp) :: value

      value = rng_uniform(rng)
   end function rng_uniform_local

end module learnbayes_hierarchical
