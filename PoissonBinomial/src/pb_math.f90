! SPDX-License-Identifier: GPL-3.0-only
module pb_math
  use pb_kinds, only : dp, pi
  implicit none
  private
  public :: normal_cdf, normal_pdf, reg_beta, reg_gamma_p, reg_gamma_q
  public :: binom_pmf, binom_cdf, poisson_pmf, poisson_cdf

contains

  pure real(dp) function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    y = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function beta_cf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 300
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = tiny(1.0_dp)/eps
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap

    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, maxit
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x / &
           ((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
           ((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del-1.0_dp) <= eps) exit
    end do
    cf = h
  end function beta_cf

  pure real(dp) function reg_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      p = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) + &
             a*log(x) + b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      p = bt*beta_cf(a,b,x)/a
    else
      p = 1.0_dp - bt*beta_cf(b,a,1.0_dp-x)/b
    end if
    p = max(0.0_dp, min(1.0_dp, p))
  end function reg_beta

  pure real(dp) function reg_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 3.0e-14_dp
    integer :: n
    real(dp) :: ap, del, summ
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x >= a + 1.0_dp) then
      p = 1.0_dp - reg_gamma_q(a,x)
      return
    end if
    ap = a
    summ = 1.0_dp/a
    del = summ
    do n = 1, maxit
      ap = ap + 1.0_dp
      del = del*x/ap
      summ = summ + del
      if (abs(del) <= abs(summ)*eps) exit
    end do
    p = summ*exp(-x+a*log(x)-log_gamma(a))
    p = max(0.0_dp, min(1.0_dp, p))
  end function reg_gamma_p

  pure real(dp) function reg_gamma_q(a, x) result(q)
    real(dp), intent(in) :: a, x
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = tiny(1.0_dp)/eps
    integer :: i
    real(dp) :: an, b, c, d, del, h
    if (x <= 0.0_dp) then
      q = 1.0_dp
      return
    end if
    if (x < a + 1.0_dp) then
      q = 1.0_dp - reg_gamma_p(a,x)
      return
    end if
    b = x + 1.0_dp - a
    c = 1.0_dp/fpmin
    d = 1.0_dp/b
    h = d
    do i = 1, maxit
      an = -real(i,dp)*(real(i,dp)-a)
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
    q = exp(-x+a*log(x)-log_gamma(a))*h
    q = max(0.0_dp, min(1.0_dp, q))
  end function reg_gamma_q

  pure real(dp) function binom_pmf(k, n, p0) result(p)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: p0
    if (k < 0 .or. k > n) then
      p = 0.0_dp
    else if (p0 <= 0.0_dp) then
      p = merge(1.0_dp, 0.0_dp, k == 0)
    else if (p0 >= 1.0_dp) then
      p = merge(1.0_dp, 0.0_dp, k == n)
    else
      p = exp(log_gamma(real(n+1,dp)) - log_gamma(real(k+1,dp)) - &
              log_gamma(real(n-k+1,dp)) + real(k,dp)*log(p0) + &
              real(n-k,dp)*log(1.0_dp-p0))
    end if
  end function binom_pmf

  pure real(dp) function binom_cdf(k, n, p0) result(p)
    integer, intent(in) :: k, n
    real(dp), intent(in) :: p0
    if (k < 0) then
      p = 0.0_dp
    else if (k >= n) then
      p = 1.0_dp
    else if (p0 <= 0.0_dp) then
      p = 1.0_dp
    else if (p0 >= 1.0_dp) then
      p = 0.0_dp
    else
      p = reg_beta(1.0_dp-p0, real(n-k,dp), real(k+1,dp))
    end if
  end function binom_cdf

  pure real(dp) function poisson_pmf(k, lambda) result(p)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    if (k < 0 .or. lambda < 0.0_dp) then
      p = 0.0_dp
    else if (lambda <= tiny(1.0_dp)) then
      p = merge(1.0_dp, 0.0_dp, k == 0)
    else
      p = exp(-lambda + real(k,dp)*log(lambda) - log_gamma(real(k+1,dp)))
    end if
  end function poisson_pmf

  pure real(dp) function poisson_cdf(k, lambda) result(p)
    integer, intent(in) :: k
    real(dp), intent(in) :: lambda
    if (k < 0) then
      p = 0.0_dp
    else if (lambda <= 0.0_dp) then
      p = 1.0_dp
    else
      p = reg_gamma_q(real(k+1,dp), lambda)
    end if
  end function poisson_cdf

end module pb_math
