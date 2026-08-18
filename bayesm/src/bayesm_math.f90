module bayesm_math
  use bayesm_kinds, only: dp, pi, sqrt2
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
    if (x >= 0.0_dp) then
      p = 1.0_dp/(1.0_dp + exp(-x))
    else
      p = exp(x)/(1.0_dp + exp(x))
    end if
  end function logistic

  pure elemental real(dp) function log1pexp(x) result(y)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      y = x + log(1.0_dp + exp(-x))
    else
      y = log(1.0_dp + exp(x))
    end if
  end function log1pexp

  pure real(dp) function logsumexp(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) == 0) then
      v = -huge(1.0_dp)
      return
    end if
    m = maxval(x)
    if (m <= -0.5_dp*huge(1.0_dp)) then
      v = m
    else
      v = m + log(sum(exp(x-m)))
    end if
  end function logsumexp

  pure elemental real(dp) function log_beta(a,b) result(v)
    real(dp), intent(in) :: a,b
    v = log_gamma(a)+log_gamma(b)-log_gamma(a+b)
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
      v = log_gamma(real(n+1,dp))
    end if
  end function log_factorial

  pure elemental real(dp) function log_choose(n,k) result(v)
    integer, intent(in) :: n,k
    if (k < 0 .or. k > n .or. n < 0) then
      v = -huge(1.0_dp)
    else
      v = log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp))
    end if
  end function log_choose

  pure elemental real(dp) function normal_quantile(p) result(x)
    ! Acklam inverse-normal approximation with one Halley correction.
    real(dp), intent(in) :: p
    real(dp), parameter :: a1=-3.969683028665376e1_dp, a2=2.209460984245205e2_dp
    real(dp), parameter :: a3=-2.759285104469687e2_dp, a4=1.383577518672690e2_dp
    real(dp), parameter :: a5=-3.066479806614716e1_dp, a6=2.506628277459239_dp
    real(dp), parameter :: b1=-5.447609879822406e1_dp, b2=1.615858368580409e2_dp
    real(dp), parameter :: b3=-1.556989798598866e2_dp, b4=6.680131188771972e1_dp
    real(dp), parameter :: b5=-1.328068155288572e1_dp
    real(dp), parameter :: c1=-7.784894002430293e-3_dp, c2=-3.223964580411365e-1_dp
    real(dp), parameter :: c3=-2.400758277161838_dp, c4=-2.549732539343734_dp
    real(dp), parameter :: c5=4.374664141464968_dp, c6=2.938163982698783_dp
    real(dp), parameter :: d1=7.784695709041462e-3_dp, d2=3.224671290700398e-1_dp
    real(dp), parameter :: d3=2.445134137142996_dp, d4=3.754408661907416_dp
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q,r,e,u,pp
    pp = min(1.0_dp-epsilon(1.0_dp), max(tiny(1.0_dp), p))
    if (pp < plow) then
      q = sqrt(-2.0_dp*log(pp))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (pp > phigh) then
      q = sqrt(-2.0_dp*log(1.0_dp-pp))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else
      q = pp-0.5_dp
      r = q*q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    end if
    e = normal_cdf(x)-pp
    u = e*sqrt(2.0_dp*pi)*exp(0.5_dp*x*x)
    x = x-u/(1.0_dp+0.5_dp*x*u)
  end function normal_quantile
end module bayesm_math
