! SPDX-License-Identifier: GPL-3.0-only
module svdnf_stats
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use svdnf_kinds, only : dp, pi
  implicit none
  private
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: gamma_pdf, gamma_cdf, poisson_pmf, binomial_pmf
  public :: seed_random, random_normal, random_gamma, random_poisson
  public :: random_bernoulli, sample_discrete, mean_value, standard_deviation

contains

  pure real(dp) function normal_pdf(x, mean, sd) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mean, sd
    real(dp) :: mu, sigma, z
    mu = 0.0_dp
    sigma = 1.0_dp
    if (present(mean)) mu = mean
    if (present(sd)) sigma = sd
    if (sigma <= 0.0_dp) then
      value = 0.0_dp
      return
    end if
    z = (x - mu) / sigma
    value = exp(-0.5_dp * z * z) / (sqrt(2.0_dp * pi) * sigma)
  end function normal_pdf

  pure real(dp) function normal_cdf(x, mean, sd) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mean, sd
    real(dp) :: mu, sigma
    mu = 0.0_dp
    sigma = 1.0_dp
    if (present(mean)) mu = mean
    if (present(sd)) sigma = sd
    if (sigma <= 0.0_dp) then
      if (x < mu) then
        value = 0.0_dp
      else
        value = 1.0_dp
      end if
      return
    end if
    value = 0.5_dp * erfc(-(x - mu) / (sigma * sqrt(2.0_dp)))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, &
      2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, &
      2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ 7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
      2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q * q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
  end function normal_quantile

  pure real(dp) function gamma_pdf(x, shape, scale) result(value)
    real(dp), intent(in) :: x, shape, scale
    if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. x < 0.0_dp) then
      value = 0.0_dp
    else if (x <= tiny(1.0_dp) .and. shape < 1.0_dp) then
      value = huge(1.0_dp)
    else if (x <= tiny(1.0_dp) .and. shape > 1.0_dp) then
      value = 0.0_dp
    else
      value = exp((shape-1.0_dp)*log(max(x,tiny(1.0_dp))) - x/scale - &
        log_gamma(shape) - shape*log(scale))
    end if
  end function gamma_pdf

  pure real(dp) function gamma_cdf(x, shape, scale) result(value)
    real(dp), intent(in) :: x, shape, scale
    integer :: n
    real(dp) :: ap, del, sum_value, b, c, d, h, an, z, gln
    integer, parameter :: max_iter = 300
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
    if (x <= 0.0_dp) then
      value = 0.0_dp
      return
    end if
    if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
      value = 0.0_dp
      return
    end if
    z = x / scale
    gln = log_gamma(shape)
    if (z < shape + 1.0_dp) then
      ap = shape
      del = 1.0_dp / shape
      sum_value = del
      do n = 1, max_iter
        ap = ap + 1.0_dp
        del = del * z / ap
        sum_value = sum_value + del
        if (abs(del) <= abs(sum_value)*eps) exit
      end do
      value = sum_value * exp(-z + shape*log(z) - gln)
    else
      b = z + 1.0_dp - shape
      c = 1.0_dp / fpmin
      d = 1.0_dp / b
      h = d
      do n = 1, max_iter
        an = -real(n,dp) * (real(n,dp)-shape)
        b = b + 2.0_dp
        d = an*d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      value = 1.0_dp - exp(-z + shape*log(z) - gln)*h
    end if
    value = min(max(value,0.0_dp),1.0_dp)
  end function gamma_cdf

  pure real(dp) function poisson_pmf(n, lambda) result(value)
    integer, intent(in) :: n
    real(dp), intent(in) :: lambda
    if (n < 0 .or. lambda < 0.0_dp) then
      value = 0.0_dp
    else if (lambda <= tiny(1.0_dp)) then
      if (n == 0) then
        value = 1.0_dp
      else
        value = 0.0_dp
      end if
    else
      value = exp(-lambda + real(n,dp)*log(lambda) - log_gamma(real(n+1,dp)))
    end if
  end function poisson_pmf

  pure real(dp) function binomial_pmf(n, size, probability) result(value)
    integer, intent(in) :: n, size
    real(dp), intent(in) :: probability
    if (n < 0 .or. n > size .or. probability < 0.0_dp .or. probability > 1.0_dp) then
      value = 0.0_dp
    else if (probability <= tiny(1.0_dp)) then
      value = merge(1.0_dp, 0.0_dp, n == 0)
    else if (probability >= 1.0_dp-tiny(1.0_dp)) then
      value = merge(1.0_dp, 0.0_dp, n == size)
    else
      value = exp(log_gamma(real(size+1,dp)) - log_gamma(real(n+1,dp)) - &
        log_gamma(real(size-n+1,dp)) + real(n,dp)*log(probability) + &
        real(size-n,dp)*log(1.0_dp-probability))
    end if
  end function binomial_pmf

  subroutine seed_random(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: values(:)
    call random_seed(size=n)
    allocate(values(n))
    do i = 1, n
      values(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
      if (values(i) == 0) values(i) = i
    end do
    call random_seed(put=values)
  end subroutine seed_random

  real(dp) function random_normal() result(value)
    real(dp) :: u1, u2
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    value = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function random_normal

  recursive real(dp) function random_gamma(shape, scale) result(value)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d, c, x, v, u
    if (shape <= 0.0_dp .or. scale < 0.0_dp) then
      value = 0.0_dp
      return
    end if
    if (scale <= tiny(1.0_dp)) then
      value = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      call random_number(u)
      value = random_gamma(shape+1.0_dp,scale)*u**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      x = random_normal()
      v = 1.0_dp + c*x
      if (v <= 0.0_dp) cycle
      v = v*v*v
      call random_number(u)
      if (u < 1.0_dp - 0.0331_dp*x**4) exit
      if (log(u) < 0.5_dp*x*x + d*(1.0_dp-v+log(v))) exit
    end do
    value = scale*d*v
  end function random_gamma

  integer function random_poisson(lambda) result(value)
    real(dp), intent(in) :: lambda
    real(dp) :: threshold, product_value, u, z
    integer :: k
    if (lambda <= 0.0_dp) then
      value = 0
    else if (lambda < 30.0_dp) then
      threshold = exp(-lambda)
      product_value = 1.0_dp
      k = 0
      do
        k = k + 1
        call random_number(u)
        product_value = product_value*u
        if (product_value <= threshold) exit
      end do
      value = k - 1
    else
      do
        z = lambda + sqrt(lambda)*random_normal()
        if (z >= 0.0_dp) exit
      end do
      value = nint(z)
    end if
  end function random_poisson

  integer function random_bernoulli(probability) result(value)
    real(dp), intent(in) :: probability
    real(dp) :: u
    call random_number(u)
    if (u < probability) then
      value = 1
    else
      value = 0
    end if
  end function random_bernoulli

  integer function sample_discrete(probabilities) result(index_value)
    real(dp), intent(in) :: probabilities(:)
    real(dp) :: u, cumulative, total
    integer :: i
    total = sum(max(probabilities,0.0_dp))
    if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) then
      index_value = 1
      return
    end if
    call random_number(u)
    u = u*total
    cumulative = 0.0_dp
    do i = 1, size(probabilities)
      cumulative = cumulative + max(probabilities(i),0.0_dp)
      if (u <= cumulative) then
        index_value = i
        return
      end if
    end do
    index_value = size(probabilities)
  end function sample_discrete

  pure real(dp) function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x)/real(size(x),dp)
    end if
  end function mean_value

  pure real(dp) function standard_deviation(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: mean_x
    if (size(x) <= 1) then
      value = 0.0_dp
    else
      mean_x = sum(x)/real(size(x),dp)
      value = sqrt(sum((x-mean_x)**2)/real(size(x)-1,dp))
    end if
  end function standard_deviation

end module svdnf_stats
