module learnbayes_discrete
   use learnbayes_distributions, only: beta_logpdf, binomial_logpmf, gamma_logpdf, normal_logpdf, poisson_logpmf
   use learnbayes_kinds, only: dp
   use learnbayes_special, only: beta_log, normal_quantile, regularized_beta
   use learnbayes_types, only: discrete_summary, interval_result, mixture_beta_result, mixture_normal_result
   implicit none
   private

   public :: beta_select
   public :: binomial_beta_mix
   public :: ctable
   public :: discint
   public :: histprior
   public :: normal_normal_mix
   public :: normal_select
   public :: pbetap
   public :: pbetat
   public :: pdisc
   public :: pdiscp
   public :: poisson_gamma_mix
   public :: prior_two_parameters
   public :: summarize_discrete

contains

   subroutine beta_select(p1, x1, p2, x2, ab, info)
      real(dp), intent(in) :: p1 !! First elicited beta-prior cumulative probability in (0,1).
      real(dp), intent(in) :: x1 !! First elicited beta-prior quantile location in (0,1).
      real(dp), intent(in) :: p2 !! Second elicited beta-prior cumulative probability in (0,1).
      real(dp), intent(in) :: x2 !! Second elicited beta-prior quantile location in (0,1).
      real(dp), intent(out) :: ab(2) !! Elicited beta shape parameters, rounded to two decimals as in LearnBayes.
      integer, intent(out) :: info !! Zero on success; nonzero when elicitation inputs are invalid.
      real(dp) :: logk(100)
      real(dp) :: prob2(100)
      real(dp) :: mean0(100)
      real(dp) :: k0
      real(dp) :: m0
      real(dp) :: target
      integer :: i
      integer :: j

      info = 0
      if (p1 <= 0.0_dp .or. p1 >= 1.0_dp .or. p2 <= 0.0_dp .or. p2 >= 1.0_dp .or. &
         x1 <= 0.0_dp .or. x1 >= 1.0_dp .or. x2 <= 0.0_dp .or. x2 >= 1.0_dp) then
         info = 1
         ab = 0.0_dp
         return
      end if
      do i = 1, 100
         logk(i) = -3.0_dp + 11.0_dp*real(i - 1, dp)/99.0_dp
         k0 = exp(logk(i))
         mean0(i) = beta_mean_for_quantile(k0, x1, p1)
         prob2(i) = regularized_beta(x2, k0*mean0(i), k0*(1.0_dp - mean0(i)))
      end do
      target = p2
      j = 1
      do i = 1, 99
         if ((prob2(i) - target)*(prob2(i + 1) - target) <= 0.0_dp) then
            j = i
            exit
         end if
      end do
      if (abs(prob2(j + 1) - prob2(j)) <= epsilon(1.0_dp)) then
         k0 = exp(logk(j))
      else
         k0 = exp(logk(j) + (target - prob2(j))*(logk(j + 1) - logk(j))/ &
            (prob2(j + 1) - prob2(j)))
      end if
      m0 = beta_mean_for_quantile(k0, x1, p1)
      ab(1) = real(nint(100.0_dp*k0*m0), dp)/100.0_dp
      ab(2) = real(nint(100.0_dp*k0*(1.0_dp - m0)), dp)/100.0_dp
   end subroutine beta_select

   function beta_mean_for_quantile(k, x, p) result(mean_value)
      real(dp), intent(in) :: k !! Positive beta precision K = alpha + beta.
      real(dp), intent(in) :: x !! Quantile location in (0,1).
      real(dp), intent(in) :: p !! Target cumulative probability at x.
      real(dp) :: mean_value
      real(dp) :: lo
      real(dp) :: hi
      real(dp) :: mid
      real(dp) :: value
      integer :: iter

      lo = 1.0e-10_dp
      hi = 1.0_dp - 1.0e-10_dp
      mid = 0.5_dp
      do iter = 1, 200
         mid = 0.5_dp*(lo + hi)
         value = regularized_beta(x, k*mid, k*(1.0_dp - mid))
         if (value < p) then
            hi = mid
         else
            lo = mid
         end if
         if (abs(value - p) < 1.0e-4_dp) exit
      end do
      mean_value = mid
   end function beta_mean_for_quantile

   subroutine normal_select(p1, x1, p2, x2, mu, sigma)
      real(dp), intent(in) :: p1 !! First elicited Gaussian cumulative probability in (0,1).
      real(dp), intent(in) :: x1 !! First elicited Gaussian quantile location.
      real(dp), intent(in) :: p2 !! Second elicited Gaussian cumulative probability in (0,1).
      real(dp), intent(in) :: x2 !! Second elicited Gaussian quantile location.
      real(dp), intent(out) :: mu !! Gaussian mean satisfying the two elicited quantiles.
      real(dp), intent(out) :: sigma !! Positive Gaussian standard deviation satisfying the quantiles.
      real(dp) :: z1
      real(dp) :: z2

      z1 = normal_quantile(p1)
      z2 = normal_quantile(p2)
      sigma = (x1 - x2)/(z1 - z2)
      mu = x1 - sigma*z1
   end subroutine normal_select

   pure function pbetap(a, b, n, s) result(prob)
      real(dp), intent(in) :: a !! Positive first beta prior shape parameter.
      real(dp), intent(in) :: b !! Positive second beta prior shape parameter.
      integer, intent(in) :: n !! Future binomial sample size.
      integer, intent(in) :: s !! Future success count in [0,n].
      real(dp) :: prob

      if (s < 0 .or. s > n) then
         prob = 0.0_dp
      else
         prob = exp(log_gamma(real(n + 1, dp)) - log_gamma(real(s + 1, dp)) - &
            log_gamma(real(n - s + 1, dp)) + beta_log(real(s, dp) + a, real(n - s, dp) + b) - beta_log(a, b))
      end if
   end function pbetap

   subroutine pbetat(p0, prior_prob, a, b, successes, failures, bf, post)
      real(dp), intent(in) :: p0 !! Point-null binomial probability being tested.
      real(dp), intent(in) :: prior_prob !! Prior probability assigned to the point-null hypothesis.
      real(dp), intent(in) :: a !! First beta shape under the alternative hypothesis.
      real(dp), intent(in) :: b !! Second beta shape under the alternative hypothesis.
      integer, intent(in) :: successes !! Observed number of binomial successes.
      integer, intent(in) :: failures !! Observed number of binomial failures.
      real(dp), intent(out) :: bf !! Bayes factor in favor of the point null versus the beta alternative.
      real(dp), intent(out) :: post !! Posterior probability of the point-null hypothesis.
      real(dp) :: lbf

      lbf = real(successes, dp)*log(p0) + real(failures, dp)*log(1.0_dp - p0) + beta_log(a, b) - &
         beta_log(a + real(successes, dp), b + real(failures, dp))
      bf = exp(lbf)
      post = prior_prob*bf/(prior_prob*bf + 1.0_dp - prior_prob)
   end subroutine pbetat

   subroutine pdisc(p, prior, successes, failures, post)
      real(dp), intent(in) :: p(:) !! Discrete candidate probabilities in [0,1].
      real(dp), intent(in) :: prior(:) !! Prior masses corresponding one-for-one with p.
      integer, intent(in) :: successes !! Observed Bernoulli successes.
      integer, intent(in) :: failures !! Observed Bernoulli failures.
      real(dp), intent(out) :: post(:) !! Normalized posterior masses corresponding one-for-one with p.
      real(dp), allocatable :: loglike(:)
      real(dp) :: p1
      real(dp) :: mx
      integer :: i

      allocate(loglike(size(p)))
      do i = 1, size(p)
         p1 = p(i)
         if (p1 <= 0.0_dp) p1 = 0.5_dp
         if (p1 >= 1.0_dp) p1 = 0.5_dp
         loglike(i) = real(successes, dp)*log(p1) + real(failures, dp)*log(1.0_dp - p1)
         if (p(i) <= 0.0_dp .and. successes > 0) loglike(i) = -999.0_dp
         if (p(i) >= 1.0_dp .and. failures > 0) loglike(i) = -999.0_dp
      end do
      mx = maxval(loglike)
      post = prior*exp(loglike - mx)
      post = post/sum(post)
   end subroutine pdisc

   subroutine pdiscp(p, probs, n, successes, pred)
      real(dp), intent(in) :: p(:) !! Discrete binomial probabilities in [0,1].
      real(dp), intent(in) :: probs(:) !! Normalized masses on p.
      integer, intent(in) :: n !! Future binomial sample size.
      integer, intent(in) :: successes(:) !! Success counts for which predictive masses are requested.
      real(dp), intent(out) :: pred(:) !! Predictive masses corresponding to each requested success count.
      integer :: i
      integer :: j

      pred = 0.0_dp
      do i = 1, size(p)
         do j = 1, size(successes)
            pred(j) = pred(j) + probs(i)*exp(binomial_logpmf(successes(j), n, p(i)))
         end do
      end do
   end subroutine pdiscp

   function ctable(y, a) result(bf)
      real(dp), intent(in) :: y(:, :) !! Observed nonnegative two-way contingency-table counts.
      real(dp), intent(in) :: a(:, :) !! Positive Dirichlet prior parameters with the same shape as y.
      real(dp) :: bf
      real(dp), allocatable :: ac(:)
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: yc(:)
      real(dp), allocatable :: yr(:)
      real(dp) :: lbf
      integer :: i
      integer :: j
      integer :: nr
      integer :: nc

      nr = size(y, 1)
      nc = size(y, 2)
      allocate(ac(nc), ar(nr), yc(nc), yr(nr))
      ac = sum(a, dim=1)
      ar = sum(a, dim=2)
      yc = sum(y, dim=1)
      yr = sum(y, dim=2)
      lbf = 0.0_dp
      do j = 1, nc
         do i = 1, nr
            lbf = lbf + log_gamma(y(i, j) + a(i, j)) - log_gamma(a(i, j))
         end do
      end do
      lbf = lbf - log_gamma(sum(y + a)) + log_gamma(sum(a))
      lbf = lbf + dirichlet_log_norm(ar - real(nc - 1, dp)) + dirichlet_log_norm(ac - real(nr - 1, dp))
      lbf = lbf - dirichlet_log_norm(yr + ar - real(nc - 1, dp))
      lbf = lbf - dirichlet_log_norm(yc + ac - real(nr - 1, dp))
      bf = exp(lbf)
   end function ctable

   pure function dirichlet_log_norm(alpha) result(value)
      real(dp), intent(in) :: alpha(:) !! Positive Dirichlet parameter vector.
      real(dp) :: value
      integer :: i

      value = -log_gamma(sum(alpha))
      do i = 1, size(alpha)
         value = value + log_gamma(alpha(i))
      end do
   end function dirichlet_log_norm

   function discint(values, probs, coverage) result(ans)
      real(dp), intent(in) :: values(:) !! Discrete support values.
      real(dp), intent(in) :: probs(:) !! Probability masses associated with values.
      real(dp), intent(in) :: coverage !! Requested minimum probability content in (0,1].
      type(interval_result) :: ans
      integer, allocatable :: idx(:)
      real(dp) :: cumulative
      integer :: i
      integer :: j
      integer :: n

      n = size(values)
      allocate(idx(n))
      do i = 1, n
         idx(i) = i
      end do
      call sort_indices_desc(probs, idx)
      cumulative = 0.0_dp
      j = 0
      do i = 1, n
         cumulative = cumulative + probs(idx(i))
         j = i
         if (cumulative >= coverage) exit
      end do
      allocate(ans%set(j))
      do i = 1, j
         ans%set(i) = values(idx(i))
      end do
      call sort_values(ans%set)
      ans%probability = cumulative
   end function discint

   function summarize_discrete(values, probs, coverage) result(ans)
      real(dp), intent(in) :: values(:) !! Discrete support values.
      real(dp), intent(in) :: probs(:) !! Probability masses associated with values.
      real(dp), intent(in), optional :: coverage !! Requested highest-mass set coverage; defaults to 0.9.
      type(discrete_summary) :: ans
      type(interval_result) :: hpd
      real(dp) :: cov

      cov = 0.9_dp
      if (present(coverage)) cov = coverage
      ans%mean = sum(values*probs)
      ans%sd = sqrt(sum((values - ans%mean)**2*probs))
      hpd = discint(values, probs, cov)
      ans%coverage = hpd%probability
      allocate(ans%set(size(hpd%set)))
      ans%set = hpd%set
   end function summarize_discrete

   subroutine normal_normal_mix(probs, normalpar, y, sigma2, result)
      real(dp), intent(in) :: probs(:) !! Prior mixture weights.
      real(dp), intent(in) :: normalpar(:, :) !! Prior normal parameters by row: mean then variance.
      real(dp), intent(in) :: y !! Observed normal response.
      real(dp), intent(in) :: sigma2 !! Known positive observation variance.
      type(mixture_normal_result), intent(out) :: result !! Updated mixture weights and posterior normal parameters.
      real(dp), allocatable :: marginal(:)
      real(dp), allocatable :: post_precision(:)
      integer :: i
      integer :: n

      n = size(probs)
      allocate(result%probs(n), result%par(n, 2), marginal(n), post_precision(n))
      post_precision = 1.0_dp/normalpar(:, 2) + 1.0_dp/sigma2
      result%par(:, 2) = 1.0_dp/post_precision
      result%par(:, 1) = (y/sigma2 + normalpar(:, 1)/normalpar(:, 2))/post_precision
      marginal = exp([(normal_logpdf(y, normalpar(i, 1), sqrt(sigma2 + normalpar(i, 2))), i=1,n)])
      result%probs = probs*marginal/sum(probs*marginal)
   end subroutine normal_normal_mix

   subroutine binomial_beta_mix(probs, betapar, successes, failures, result)
      real(dp), intent(in) :: probs(:) !! Prior mixture weights.
      real(dp), intent(in) :: betapar(:, :) !! Prior beta parameters by row: alpha then beta.
      integer, intent(in) :: successes !! Observed binomial successes.
      integer, intent(in) :: failures !! Observed binomial failures.
      type(mixture_beta_result), intent(out) :: result !! Updated mixture weights and posterior beta parameters.
      real(dp), allocatable :: marginal(:)
      real(dp) :: p
      integer :: i
      integer :: n

      n = size(probs)
      allocate(result%probs(n), result%par(n, 2), marginal(n))
      result%par(:, 1) = betapar(:, 1) + real(successes, dp)
      result%par(:, 2) = betapar(:, 2) + real(failures, dp)
      do i = 1, n
         p = result%par(i, 1)/sum(result%par(i, :))
         marginal(i) = exp(binomial_logpmf(successes, successes + failures, p) + &
            beta_logpdf(p, betapar(i, 1), betapar(i, 2)) - beta_logpdf(p, result%par(i, 1), result%par(i, 2)))
      end do
      result%probs = probs*marginal/sum(probs*marginal)
   end subroutine binomial_beta_mix

   subroutine poisson_gamma_mix(probs, gammapar, y, exposure, result)
      real(dp), intent(in) :: probs(:) !! Prior mixture weights.
      real(dp), intent(in) :: gammapar(:, :) !! Prior gamma parameters by row: shape then rate.
      integer, intent(in) :: y(:) !! Observed Poisson counts.
      real(dp), intent(in) :: exposure(:) !! Positive exposure values paired with y.
      type(mixture_beta_result), intent(out) :: result !! Updated weights and posterior gamma shape/rate parameters.
      real(dp), allocatable :: marginal(:)
      real(dp) :: lambda
      real(dp) :: loglike
      integer :: i
      integer :: j
      integer :: n

      n = size(probs)
      allocate(result%probs(n), result%par(n, 2), marginal(n))
      result%par(:, 1) = gammapar(:, 1) + real(sum(y), dp)
      result%par(:, 2) = gammapar(:, 2) + sum(exposure)
      do i = 1, n
         lambda = result%par(i, 1)/result%par(i, 2)
         loglike = 0.0_dp
         do j = 1, size(y)
            loglike = loglike + poisson_logpmf(y(j), lambda*exposure(j))
         end do
         marginal(i) = exp(loglike + gamma_logpdf(lambda, gammapar(i, 1), gammapar(i, 2)) - &
            gamma_logpdf(lambda, result%par(i, 1), result%par(i, 2)))
      end do
      result%probs = probs*marginal/sum(probs*marginal)
   end subroutine poisson_gamma_mix

   subroutine prior_two_parameters(parameter1, parameter2, prior)
      real(dp), intent(in) :: parameter1(:) !! First parameter grid; values define prior matrix rows.
      real(dp), intent(in) :: parameter2(:) !! Second parameter grid; values define prior matrix columns.
      real(dp), intent(out) :: prior(:, :) !! Uniform joint prior matrix over the Cartesian product of the grids.

      prior = 1.0_dp/real(size(parameter1)*size(parameter2), dp)
   end subroutine prior_two_parameters

   subroutine histprior(p, midpts, prob, value)
      real(dp), intent(in) :: p(:) !! Points at which the piecewise-constant histogram prior is evaluated.
      real(dp), intent(in) :: midpts(:) !! Equally spaced bin midpoints.
      real(dp), intent(in) :: prob(:) !! Bin heights or masses returned for points in each bin.
      real(dp), intent(out) :: value(:) !! Histogram-prior values corresponding one-for-one with p.
      real(dp), allocatable :: lo(:)
      real(dp) :: binwidth
      integer :: i
      integer :: j

      binwidth = midpts(2) - midpts(1)
      allocate(lo(size(midpts)))
      lo = real(nint(10000.0_dp*(midpts - 0.5_dp*binwidth)), dp)/10000.0_dp
      do i = 1, size(p)
         j = count(p(i) >= lo)
         j = max(1, min(size(prob), j))
         value(i) = prob(j)
      end do
   end subroutine histprior

   subroutine sort_indices_desc(x, idx)
      real(dp), intent(in) :: x(:) !! Values whose descending order determines the index permutation.
      integer, intent(inout) :: idx(:) !! Index permutation sorted so x(idx) is nonincreasing.
      integer :: i
      integer :: j
      integer :: key

      do i = 2, size(idx)
         key = idx(i)
         j = i - 1
         do while (j >= 1)
            if (x(idx(j)) >= x(key)) exit
            idx(j + 1) = idx(j)
            j = j - 1
         end do
         idx(j + 1) = key
      end do
   end subroutine sort_indices_desc

   subroutine sort_values(x)
      real(dp), intent(inout) :: x(:) !! Real vector sorted into nondecreasing order in place.
      real(dp) :: key
      integer :: i
      integer :: j

      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_values

end module learnbayes_discrete
