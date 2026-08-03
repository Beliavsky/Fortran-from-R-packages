! SPDX-License-Identifier: GPL-3.0-only
module imputefin_math
  use imputefin_kinds, only : dp
  implicit none
  private
  public :: normal_cdf, student_t_cdf, golden_minimize, digamma_dp, log_gamma_dp
contains
  elemental function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  function betacf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    real(dp) :: cf, qab, qap, qam, c, d, h, aa, del
    integer :: m, m2
    real(dp), parameter :: fpmin = 1.0e-300_dp, eps = 3.0e-14_dp
    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab * x / qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d
    do m = 1, 300
      m2 = 2 * m
      aa = real(m,dp) * (b - real(m,dp)) * x / ((qam + real(m2,dp)) * (a + real(m2,dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      h = h * d * c
      aa = -(a + real(m,dp)) * (qab + real(m,dp)) * x / &
           ((a + real(m2,dp)) * (qap + real(m2,dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) <= eps) exit
    end do
    cf = h
  end function betacf

  function regularized_beta(x, a, b) result(v)
    real(dp), intent(in) :: x, a, b
    real(dp) :: v, bt
    if (x <= 0.0_dp) then
      v = 0.0_dp
    else if (x >= 1.0_dp) then
      v = 1.0_dp
    else
      bt = exp(log_gamma(a+b) - log_gamma(a) - log_gamma(b) + a*log(x) + b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
        v = bt * betacf(a,b,x) / a
      else
        v = 1.0_dp - bt * betacf(b,a,1.0_dp-x) / b
      end if
    end if
  end function regularized_beta

  function student_t_cdf(x, nu) result(p)
    real(dp), intent(in) :: x, nu
    real(dp) :: p, z, ib
    if (nu <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (nu > 1.0e5_dp) then
      p = normal_cdf(x)
      return
    end if
    z = nu / (nu + x*x)
    ib = regularized_beta(z, 0.5_dp*nu, 0.5_dp)
    if (x >= 0.0_dp) then
      p = 1.0_dp - 0.5_dp*ib
    else
      p = 0.5_dp*ib
    end if
  end function student_t_cdf

  pure function log_gamma_dp(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v
    v = log_gamma(x)
  end function log_gamma_dp

  function digamma_dp(xin) result(v)
    real(dp), intent(in) :: xin
    real(dp) :: v, x, inv, inv2
    x = xin
    v = 0.0_dp
    if (x <= 0.0_dp) then
      v = huge(1.0_dp)
      return
    end if
    do while (x < 8.0_dp)
      v = v - 1.0_dp/x
      x = x + 1.0_dp
    end do
    inv = 1.0_dp/x
    inv2 = inv*inv
    v = v + log(x) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp - &
        inv2*(1.0_dp/252.0_dp - inv2/240.0_dp)))
  end function digamma_dp

  subroutine golden_minimize(f, a, b, xmin, fmin, tol, maxiter)
    interface
      function f(x) result(y)
        import dp
        real(dp), intent(in) :: x
        real(dp) :: y
      end function f
    end interface
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: xmin, fmin
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    real(dp) :: left, right, x1, x2, f1, f2, eps, gr
    integer :: i, nmax
    eps = 1.0e-8_dp
    if (present(tol)) eps = tol
    nmax = 200
    if (present(maxiter)) nmax = maxiter
    gr = (sqrt(5.0_dp)-1.0_dp)/2.0_dp
    left = a
    right = b
    x1 = right - gr*(right-left)
    x2 = left + gr*(right-left)
    f1 = f(x1)
    f2 = f(x2)
    do i = 1, nmax
      if (abs(right-left) <= eps*(1.0_dp+abs(x1)+abs(x2))) exit
      if (f1 > f2) then
        left = x1
        x1 = x2
        f1 = f2
        x2 = left + gr*(right-left)
        f2 = f(x2)
      else
        right = x2
        x2 = x1
        f2 = f1
        x1 = right - gr*(right-left)
        f1 = f(x1)
      end if
    end do
    if (f1 <= f2) then
      xmin = x1
      fmin = f1
    else
      xmin = x2
      fmin = f2
    end if
  end subroutine golden_minimize
end module imputefin_math
