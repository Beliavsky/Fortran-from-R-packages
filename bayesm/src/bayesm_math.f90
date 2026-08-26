module bayesm_math
  use bayesm_kinds, only: dp, pi, sqrt2
  use r_transforms, only: r_log1pexp, r_logistic
  use r_distributions, only: r_qnorm
  use r_stability, only: r_log_sum_exp
  use r_special, only: r_log_beta, r_log_choose, r_log_factorial
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile, logsumexp, logistic, log1pexp
  public :: log_beta, log_multivariate_gamma, log_factorial, log_choose, clamp_probability
contains
  pure elemental real(dp) function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  pure elemental real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp*erfc(-x/sqrt2)
  end function normal_cdf

  pure elemental real(dp) function clamp_probability(p) result(q)
    real(dp), intent(in) :: p
    q = min(1.0_dp, max(0.0_dp, p))
  end function clamp_probability

  pure elemental real(dp) function logistic(x) result(p)
    real(dp), intent(in) :: x
    p = r_logistic(x)
  end function logistic

  pure elemental real(dp) function log1pexp(x) result(y)
    real(dp), intent(in) :: x
    y = r_log1pexp(x)
  end function log1pexp

  pure real(dp) function logsumexp(x) result(v)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      v = -huge(1.0_dp)
    else
      v = r_log_sum_exp(x)
    end if
  end function logsumexp

  pure elemental real(dp) function log_beta(a,b) result(v)
    real(dp), intent(in) :: a,b
    v = r_log_beta(a,b)
  end function log_beta

  pure real(dp) function log_multivariate_gamma(a,p) result(v)
    real(dp), intent(in) :: a
    integer, intent(in) :: p
    integer :: j
    v = 0.25_dp*real(p*(p-1),dp)*log(pi)
    do j=1,p
      v = v + log_gamma(a + 0.5_dp*real(1-j,dp))
    end do
  end function log_multivariate_gamma

  pure elemental real(dp) function log_factorial(n) result(v)
    integer, intent(in) :: n
    if (n < 0) then
      v = huge(1.0_dp)
    else
      v = r_log_factorial(n)
    end if
  end function log_factorial

  pure elemental real(dp) function log_choose(n,k) result(v)
    integer, intent(in) :: n,k
    if (k < 0 .or. k > n .or. n < 0) then
      v = -huge(1.0_dp)
    else
      v = r_log_choose(n, k)
    end if
  end function log_choose

  pure elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: pp
    pp = min(1.0_dp-epsilon(1.0_dp), max(tiny(1.0_dp), p))
    x = r_qnorm(pp)
  end function normal_quantile
end module bayesm_math
