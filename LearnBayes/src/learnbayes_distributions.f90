module learnbayes_distributions
   use learnbayes_kinds, only: dp
   use learnbayes_linalg, only: cholesky_lower, inverse_matrix
   use learnbayes_rng, only: rng_chisq, rng_gamma, rng_normal, rng_state, rng_uniform
   use learnbayes_special, only: beta_log, normal_cdf, normal_quantile
   implicit none
   private

   real(dp), parameter :: pi_dp = acos(-1.0_dp)

   public :: beta_logpdf
   public :: binomial_logpmf
   public :: dmnorm
   public :: dmt
   public :: gamma_logpdf
   public :: normal_logpdf
   public :: poisson_logpmf
   public :: rdirichlet
   public :: rigamma
   public :: rmnorm
   public :: rmt
   public :: rtruncated_normal
   public :: student_t_logpdf

contains

   pure function normal_logpdf(x, mean, sd) result(value)
      real(dp), intent(in) :: x !! Observation whose Gaussian log density is requested.
      real(dp), intent(in) :: mean !! Gaussian location parameter.
      real(dp), intent(in) :: sd !! Positive Gaussian standard deviation.
      real(dp) :: value
      real(dp) :: z

      if (sd <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      z = (x - mean)/sd
      value = -0.5_dp*log(2.0_dp*pi_dp) - log(sd) - 0.5_dp*z*z
   end function normal_logpdf

   pure function gamma_logpdf(x, shape, rate) result(value)
      real(dp), intent(in) :: x !! Positive observation whose gamma log density is requested.
      real(dp), intent(in) :: shape !! Positive gamma shape parameter.
      real(dp), intent(in) :: rate !! Positive gamma rate parameter.
      real(dp) :: value

      if (x <= 0.0_dp .or. shape <= 0.0_dp .or. rate <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      value = shape*log(rate) - log_gamma(shape) + (shape - 1.0_dp)*log(x) - rate*x
   end function gamma_logpdf

   pure function beta_logpdf(x, a, b) result(value)
      real(dp), intent(in) :: x !! Probability-valued observation in the open interval (0, 1).
      real(dp), intent(in) :: a !! Positive first beta shape parameter.
      real(dp), intent(in) :: b !! Positive second beta shape parameter.
      real(dp) :: value

      if (x <= 0.0_dp .or. x >= 1.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      value = (a - 1.0_dp)*log(x) + (b - 1.0_dp)*log(1.0_dp - x) - beta_log(a, b)
   end function beta_logpdf

   pure function poisson_logpmf(y, lambda) result(value)
      integer, intent(in) :: y !! Nonnegative Poisson count.
      real(dp), intent(in) :: lambda !! Nonnegative Poisson mean.
      real(dp) :: value

      if (y < 0 .or. lambda < 0.0_dp) then
         value = -huge(1.0_dp)
      else if (lambda <= 0.0_dp) then
         value = merge(0.0_dp, -huge(1.0_dp), y == 0)
      else
         value = real(y, dp)*log(lambda) - lambda - log_gamma(real(y + 1, dp))
      end if
   end function poisson_logpmf

   pure function binomial_logpmf(s, n, p) result(value)
      integer, intent(in) :: s !! Number of successes in [0, n].
      integer, intent(in) :: n !! Nonnegative number of Bernoulli trials.
      real(dp), intent(in) :: p !! Success probability in [0, 1].
      real(dp) :: value

      if (n < 0 .or. s < 0 .or. s > n .or. p < 0.0_dp .or. p > 1.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      if (p <= 0.0_dp) then
         value = merge(0.0_dp, -huge(1.0_dp), s == 0)
         return
      else if (p >= 1.0_dp) then
         value = merge(0.0_dp, -huge(1.0_dp), s == n)
         return
      end if
      value = log_gamma(real(n + 1, dp)) - log_gamma(real(s + 1, dp)) - &
         log_gamma(real(n - s + 1, dp)) + real(s, dp)*log(p) + real(n - s, dp)*log(1.0_dp - p)
   end function binomial_logpmf

   pure function student_t_logpdf(x, df) result(value)
      real(dp), intent(in) :: x !! Standardized observation whose Student-t log density is requested.
      real(dp), intent(in) :: df !! Positive Student-t degrees of freedom.
      real(dp) :: value

      if (df <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      value = log_gamma(0.5_dp*(df + 1.0_dp)) - log_gamma(0.5_dp*df) - &
         0.5_dp*log(df*pi_dp) - 0.5_dp*(df + 1.0_dp)*log(1.0_dp + x*x/df)
   end function student_t_logpdf

   function dmnorm(x, mean, varcov, log_density) result(value)
      real(dp), intent(in) :: x(:) !! Point at which the multivariate-normal density is evaluated.
      real(dp), intent(in) :: mean(:) !! Mean vector with the same length as x.
      real(dp), intent(in) :: varcov(:, :) !! Symmetric positive-definite covariance matrix.
      logical, intent(in), optional :: log_density !! If true return log density; otherwise return density.
      real(dp) :: value
      real(dp), allocatable :: inv(:, :)
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: l(:, :)
      real(dp) :: logdet
      real(dp) :: logpdf
      integer :: i
      integer :: info
      integer :: d
      logical :: want_log

      d = size(x)
      if (size(mean) /= d .or. size(varcov, 1) /= d .or. size(varcov, 2) /= d) then
         value = 0.0_dp
         return
      end if
      allocate(inv(d, d), delta(d), l(d, d))
      call inverse_matrix(varcov, inv, info)
      if (info /= 0) then
         value = 0.0_dp
         return
      end if
      call cholesky_lower(varcov, l, info)
      if (info /= 0) then
         value = 0.0_dp
         return
      end if
      logdet = 0.0_dp
      do i = 1, d
         logdet = logdet + 2.0_dp*log(l(i, i))
      end do
      delta = x - mean
      logpdf = -0.5_dp*(real(d, dp)*log(2.0_dp*pi_dp) + logdet + dot_product(delta, matmul(inv, delta)))
      want_log = .false.
      if (present(log_density)) want_log = log_density
      if (want_log) then
         value = logpdf
      else
         value = exp(logpdf)
      end if
   end function dmnorm

   function dmt(x, mean, s, df, log_density) result(value)
      real(dp), intent(in) :: x(:) !! Point at which the multivariate-t density is evaluated.
      real(dp), intent(in) :: mean(:) !! Location vector with the same length as x.
      real(dp), intent(in) :: s(:, :) !! Symmetric positive-definite multivariate-t scale matrix.
      real(dp), intent(in) :: df !! Positive degrees of freedom; huge values approach the Gaussian density.
      logical, intent(in), optional :: log_density !! If true return log density; otherwise return density.
      real(dp) :: value
      real(dp), allocatable :: inv(:, :)
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: l(:, :)
      real(dp) :: logdet
      real(dp) :: logpdf
      real(dp) :: q
      integer :: d
      integer :: i
      integer :: info
      logical :: want_log

      if (df >= 0.5_dp*huge(1.0_dp)) then
         value = dmnorm(x, mean, s, log_density)
         return
      end if
      d = size(x)
      if (df <= 0.0_dp .or. size(mean) /= d .or. size(s, 1) /= d .or. size(s, 2) /= d) then
         value = 0.0_dp
         return
      end if
      allocate(inv(d, d), delta(d), l(d, d))
      call inverse_matrix(s, inv, info)
      if (info /= 0) then
         value = 0.0_dp
         return
      end if
      call cholesky_lower(s, l, info)
      if (info /= 0) then
         value = 0.0_dp
         return
      end if
      logdet = 0.0_dp
      do i = 1, d
         logdet = logdet + 2.0_dp*log(l(i, i))
      end do
      delta = x - mean
      q = dot_product(delta, matmul(inv, delta))
      logpdf = log_gamma(0.5_dp*(df + real(d, dp))) - log_gamma(0.5_dp*df) - &
         0.5_dp*(real(d, dp)*log(pi_dp*df) + logdet) - &
         0.5_dp*(df + real(d, dp))*log(1.0_dp + q/df)
      want_log = .false.
      if (present(log_density)) want_log = log_density
      if (want_log) then
         value = logpdf
      else
         value = exp(logpdf)
      end if
   end function dmt

   subroutine rmnorm(rng, n, mean, varcov, draws, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for all normal draws.
      integer, intent(in) :: n !! Number of independent multivariate-normal draws requested.
      real(dp), intent(in) :: mean(:) !! Mean vector of length d.
      real(dp), intent(in) :: varcov(:, :) !! Symmetric positive-definite d-by-d covariance matrix.
      real(dp), intent(out) :: draws(:, :) !! Output matrix shaped (n,d), one draw per row.
      integer, intent(out) :: info !! Zero on success; nonzero if dimensions or Cholesky factorization fail.
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: z(:)
      integer :: d
      integer :: i
      integer :: j

      d = size(mean)
      info = 0
      if (size(varcov, 1) /= d .or. size(varcov, 2) /= d) then
         info = -1
         return
      end if
      if (size(draws, 1) /= n .or. size(draws, 2) /= d) then
         info = -2
         return
      end if
      allocate(l(d, d), z(d))
      call cholesky_lower(varcov, l, info)
      if (info /= 0) return
      do i = 1, n
         do j = 1, d
            z(j) = rng_normal(rng)
         end do
         draws(i, :) = mean + matmul(l, z)
      end do
   end subroutine rmnorm

   subroutine rmt(rng, n, mean, s, df, draws, info)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for Gaussian and chi-square variates.
      integer, intent(in) :: n !! Number of independent multivariate-t draws requested.
      real(dp), intent(in) :: mean(:) !! Location vector of length d.
      real(dp), intent(in) :: s(:, :) !! Symmetric positive-definite d-by-d scale matrix.
      real(dp), intent(in) :: df !! Positive degrees of freedom; huge values select Gaussian draws.
      real(dp), intent(out) :: draws(:, :) !! Output matrix shaped (n,d), one draw per row.
      integer, intent(out) :: info !! Zero on success; nonzero if dimensions or covariance factorization fail.
      real(dp), allocatable :: zeros(:)
      real(dp), allocatable :: z(:, :)
      real(dp) :: scale
      integer :: d
      integer :: i

      d = size(mean)
      allocate(zeros(d), z(n, d))
      zeros = 0.0_dp
      call rmnorm(rng, n, zeros, s, z, info)
      if (info /= 0) return
      do i = 1, n
         if (df >= 0.5_dp*huge(1.0_dp)) then
            scale = 1.0_dp
         else
            scale = sqrt(rng_chisq(rng, df)/df)
         end if
         draws(i, :) = mean + z(i, :)/scale
      end do
   end subroutine rmt

   subroutine rdirichlet(rng, n, par, draws)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the component gamma draws.
      integer, intent(in) :: n !! Number of independent Dirichlet draws requested.
      real(dp), intent(in) :: par(:) !! Positive Dirichlet concentration parameters.
      real(dp), intent(out) :: draws(:, :) !! Output matrix shaped (n,size(par)), each row summing to one.
      integer :: i
      integer :: j
      real(dp) :: total

      do i = 1, n
         total = 0.0_dp
         do j = 1, size(par)
            draws(i, j) = rng_gamma(rng, par(j), 1.0_dp)
            total = total + draws(i, j)
         end do
         if (total > 0.0_dp) draws(i, :) = draws(i, :)/total
      end do
   end subroutine rdirichlet

   function rigamma(rng, a, b) result(value)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the inverse-gamma draw.
      real(dp), intent(in) :: a !! Positive inverse-gamma shape parameter.
      real(dp), intent(in) :: b !! Positive inverse-gamma rate parameter appearing as exp(-b/x).
      real(dp) :: value

      value = 1.0_dp/rng_gamma(rng, a, b)
   end function rigamma

   function rtruncated_normal(rng, lo, hi, mean, sd) result(value)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for inverse-CDF sampling.
      real(dp), intent(in) :: lo !! Lower truncation bound; may be a large negative sentinel for minus infinity.
      real(dp), intent(in) :: hi !! Upper truncation bound; may be a large positive sentinel for plus infinity.
      real(dp), intent(in) :: mean !! Mean of the untruncated Gaussian distribution.
      real(dp), intent(in) :: sd !! Positive standard deviation of the untruncated Gaussian distribution.
      real(dp) :: value
      real(dp) :: plo
      real(dp) :: phi
      real(dp) :: u

      plo = normal_cdf(lo, mean, sd)
      phi = normal_cdf(hi, mean, sd)
      u = plo + rng_uniform(rng)*(phi - plo)
      value = normal_quantile(u, mean, sd)
   end function rtruncated_normal

end module learnbayes_distributions
