! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module dowd_math
  use dowd_kinds, only: dp, pi, sqrt_two_pi
  implicit none
  private

  public :: mean_value, sample_sd, sample_skewness, sample_kurtosis
  public :: sort_in_place, sorted_copy, quantile_linear
  public :: normal_pdf, normal_cdf, normal_quantile
  public :: student_t_pdf, student_t_cdf, student_t_quantile
  public :: chi_square_cdf, regularized_beta, regularized_gamma_p
  public :: binomial_cdf, clamp_probability
  public :: random_normal, random_lognormal, random_bernoulli
  public :: set_random_seed_scalar

contains

  pure real(dp) function clamp_probability(p) result(q)
    real(dp), intent(in) :: p
    q = min(max(p, 8.0_dp*epsilon(1.0_dp)), 1.0_dp - 8.0_dp*epsilon(1.0_dp))
  end function clamp_probability

  pure real(dp) function mean_value(x) result(m)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      m = 0.0_dp
    else
      m = sum(x)/real(size(x), dp)
    end if
  end function mean_value

  pure real(dp) function sample_sd(x) result(s)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) < 2) then
      s = 0.0_dp
    else
      m = mean_value(x)
      s = sqrt(sum((x-m)**2)/real(size(x)-1, dp))
    end if
  end function sample_sd

  pure real(dp) function sample_skewness(x) result(skew)
    real(dp), intent(in) :: x(:)
    real(dp) :: m, s
    integer :: n
    n = size(x)
    if (n < 3) then
      skew = 0.0_dp
      return
    end if
    m = mean_value(x)
    s = sample_sd(x)
    if (s <= 0.0_dp) then
      skew = 0.0_dp
    else
      skew = real(n,dp)/real((n-1)*(n-2),dp) * sum(((x-m)/s)**3)
    end if
  end function sample_skewness

  pure real(dp) function sample_kurtosis(x) result(kurt)
    real(dp), intent(in) :: x(:)
    real(dp) :: m, s, a, b
    integer :: n
    n = size(x)
    if (n < 4) then
      kurt = 3.0_dp
      return
    end if
    m = mean_value(x)
    s = sample_sd(x)
    if (s <= 0.0_dp) then
      kurt = 3.0_dp
    else
      a = real(n*(n+1),dp)/real((n-1)*(n-2)*(n-3),dp)
      b = 3.0_dp*real((n-1)*(n-1),dp)/real((n-2)*(n-3),dp)
      kurt = a*sum(((x-m)/s)**4) - b + 3.0_dp
    end if
  end function sample_kurtosis

  subroutine sort_in_place(x)
    real(dp), intent(inout) :: x(:)
    if (size(x) > 1) call quicksort(x, 1, size(x))
  end subroutine sort_in_place

  recursive subroutine quicksort(x, lo, hi)
    real(dp), intent(inout) :: x(:)
    integer, intent(in) :: lo, hi
    integer :: i, j
    real(dp) :: pivot, temp
    if (lo >= hi) return
    pivot = x((lo+hi)/2)
    i = lo
    j = hi
    do
      do while (x(i) < pivot)
        i = i + 1
      end do
      do while (x(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        temp = x(i)
        x(i) = x(j)
        x(j) = temp
        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do
    if (lo < j) call quicksort(x, lo, j)
    if (i < hi) call quicksort(x, i, hi)
  end subroutine quicksort

  function sorted_copy(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: y(:)
    y = x
    call sort_in_place(y)
  end function sorted_copy

  real(dp) function quantile_linear(x, p) result(q)
    real(dp), intent(in) :: x(:), p
    real(dp), allocatable :: y(:)
    real(dp) :: h, w
    integer :: lo, hi, n
    n = size(x)
    if (n == 0) error stop "quantile_linear: empty input"
    if (p < 0.0_dp .or. p > 1.0_dp) error stop "quantile_linear: p outside [0,1]"
    y = sorted_copy(x)
    if (n == 1) then
      q = y(1)
      return
    end if
    h = 1.0_dp + real(n-1,dp)*p
    lo = max(1, min(n, int(floor(h))))
    hi = max(1, min(n, int(ceiling(h))))
    w = h-real(lo,dp)
    q = (1.0_dp-w)*y(lo) + w*y(hi)
  end function quantile_linear

  elemental real(dp) function normal_pdf(x) result(f)
    real(dp), intent(in) :: x
    f = exp(-0.5_dp*x*x)/sqrt_two_pi
  end function normal_pdf

  elemental real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: q, r
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
       4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
       2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
    if (p > 0.0_dp .and. p < 1.0_dp) then
      x = x - (normal_cdf(x)-p)/max(normal_pdf(x), tiny(1.0_dp))
    end if
  end function normal_quantile

  pure real(dp) function beta_continued_fraction(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: max_iter = 300
    real(dp), parameter :: eps = 4.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
    real(dp) :: aa, c, d, del, h, qab, qam, qap
    integer :: m, m2
    qab = a+b
    qap = a+1.0_dp
    qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp-qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, max_iter
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
           ((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del-1.0_dp) <= eps) exit
    end do
    cf = h
  end function beta_continued_fraction

  pure real(dp) function regularized_beta(x, a, b) result(value)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) error stop "regularized_beta: invalid shape"
    if (x <= 0.0_dp) then
      value = 0.0_dp
    else if (x >= 1.0_dp) then
      value = 1.0_dp
    else
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
        value = bt*beta_continued_fraction(a,b,x)/a
      else
        value = 1.0_dp-bt*beta_continued_fraction(b,a,1.0_dp-x)/b
      end if
    end if
  end function regularized_beta

  elemental real(dp) function student_t_pdf(x, nu) result(f)
    real(dp), intent(in) :: x, nu
    if (nu <= 0.0_dp) then
      f = 0.0_dp
    else
      f = exp(log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)) / &
          sqrt(nu*pi) * (1.0_dp+x*x/nu)**(-0.5_dp*(nu+1.0_dp))
    end if
  end function student_t_pdf

  real(dp) function student_t_cdf(x, nu) result(p)
    real(dp), intent(in) :: x, nu
    real(dp) :: ib
    if (nu <= 0.0_dp) error stop "student_t_cdf: nu must be positive"
    if (abs(x) <= tiny(1.0_dp)) then
      p = 0.5_dp
    else
      ib = regularized_beta(nu/(nu+x*x), 0.5_dp*nu, 0.5_dp)
      if (x > 0.0_dp) then
        p = 1.0_dp-0.5_dp*ib
      else
        p = 0.5_dp*ib
      end if
    end if
  end function student_t_cdf

  real(dp) function student_t_quantile(p, nu) result(x)
    real(dp), intent(in) :: p, nu
    real(dp) :: lo, hi, mid, c
    integer :: iter
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    lo = -1.0_dp
    hi = 1.0_dp
    do while (student_t_cdf(lo,nu) > p)
      lo = 2.0_dp*lo
    end do
    do while (student_t_cdf(hi,nu) < p)
      hi = 2.0_dp*hi
    end do
    do iter = 1, 120
      mid = 0.5_dp*(lo+hi)
      c = student_t_cdf(mid,nu)
      if (c < p) then
        lo = mid
      else
        hi = mid
      end if
      if (hi-lo <= 32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(mid))) exit
    end do
    x = 0.5_dp*(lo+hi)
  end function student_t_quantile

  pure real(dp) function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: max_iter = 1000
    real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
    real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
    real(dp) :: ap, del, sumv, b, c, d, h, an
    integer :: i
    if (a <= 0.0_dp) error stop "regularized_gamma_p: a must be positive"
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x < a+1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do i = 1, max_iter
        ap = ap+1.0_dp
        del = del*x/ap
        sumv = sumv+del
        if (abs(del) < abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x+1.0_dp-a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do i = 1, max_iter
        an = -real(i,dp)*(real(i,dp)-a)
        b = b+2.0_dp
        d = an*d+b
        if (abs(d) < fpmin) d = fpmin
        c = b+an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) < eps) exit
      end do
      p = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p = min(max(p,0.0_dp),1.0_dp)
  end function regularized_gamma_p

  pure real(dp) function chi_square_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    if (x <= 0.0_dp) then
      p = 0.0_dp
    else
      p = regularized_gamma_p(0.5_dp*df, 0.5_dp*x)
    end if
  end function chi_square_cdf

  real(dp) function binomial_cdf(k, n, p) result(cdf)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: p
    integer :: i
    real(dp) :: term, q
    if (k < 0) then
      cdf = 0.0_dp
      return
    else if (k >= n) then
      cdf = 1.0_dp
      return
    end if
    if (p <= 0.0_dp) then
      cdf = 1.0_dp
      return
    else if (p >= 1.0_dp) then
      cdf = 0.0_dp
      return
    end if
    q = 1.0_dp-p
    term = q**n
    cdf = term
    do i = 1, k
      term = term*real(n-i+1,dp)/real(i,dp)*p/q
      cdf = cdf+term
    end do
    cdf = min(max(cdf,0.0_dp),1.0_dp)
  end function binomial_cdf

  subroutine set_random_seed_scalar(seed_value)
    integer, intent(in) :: seed_value
    integer :: n, i
    integer, allocatable :: seed(:)
    call random_seed(size=n)
    allocate(seed(n))
    do i = 1, n
      seed(i) = modulo(seed_value + 104729*i, huge(1)-1)
      if (seed(i) == 0) seed(i) = i
    end do
    call random_seed(put=seed)
  end subroutine set_random_seed_scalar

  real(dp) function random_normal(mu, sigma) result(x)
    real(dp), intent(in), optional :: mu, sigma
    real(dp) :: u1, u2, m, s
    m = 0.0_dp
    s = 1.0_dp
    if (present(mu)) m = mu
    if (present(sigma)) s = sigma
    call random_number(u1)
    call random_number(u2)
    u1 = max(u1, tiny(1.0_dp))
    x = m+s*sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function random_normal

  real(dp) function random_lognormal(mu, sigma) result(x)
    real(dp), intent(in) :: mu, sigma
    x = exp(random_normal(mu,sigma))
  end function random_lognormal

  integer function random_bernoulli(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: u
    call random_number(u)
    x = merge(1,0,u < p)
  end function random_bernoulli

end module dowd_math
