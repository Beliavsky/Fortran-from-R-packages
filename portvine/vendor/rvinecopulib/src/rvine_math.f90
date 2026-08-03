! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
module rvine_math
  use rvine_kinds, only : dp, pi, sqrt_two, eps_prob, clamp_prob
  implicit none
  private
  public :: normal_cdf, normal_pdf, normal_quantile
  public :: student_cdf, student_pdf, student_quantile
  public :: regularized_beta, integrate_unit, integrate_interval, gauss_legendre_rule
  public :: seed_rng, random_normal, random_uniform
  public :: sort_real

  abstract interface
    function scalar_fun(x) result(y)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: y
    end function scalar_fun
  end interface

contains

  pure elemental real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt_two)
  end function normal_cdf

  pure elemental real(dp) function normal_pdf(x) result(p)
    real(dp), intent(in) :: x
    p = exp(-0.5_dp * x * x) / sqrt(2.0_dp * pi)
  end function normal_pdf

  pure elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: q, r, pp
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
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow

    pp = clamp_prob(p)
    if (pp < plow) then
      q = sqrt(-2.0_dp * log(pp))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (pp <= phigh) then
      q = pp - 0.5_dp
      r = q * q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp - pp))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
    x = x - (normal_cdf(x) - pp) / max(normal_pdf(x), tiny(1.0_dp))
  end function normal_quantile

  pure real(dp) function betacf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: maxit = 300
    real(dp), parameter :: fpmin = 1.0e-300_dp, tol = 3.0e-14_dp
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap
    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab * x / qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d
    do m = 1, maxit
      m2 = 2 * m
      aa = real(m, dp) * (b - real(m, dp)) * x / &
           ((qam + real(m2, dp)) * (a + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      h = h * d * c
      aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
           ((a + real(m2, dp)) * (qap + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) <= tol) exit
    end do
    cf = h
  end function betacf

  pure real(dp) function regularized_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt, xx
    xx = min(1.0_dp, max(0.0_dp, x))
    if (xx <= 0.0_dp) then
      p = 0.0_dp
      return
    else if (xx >= 1.0_dp) then
      p = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a+b) - log_gamma(a) - log_gamma(b) + &
             a*log(xx) + b*log(1.0_dp-xx))
    if (xx < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      p = bt * betacf(a, b, xx) / a
    else
      p = 1.0_dp - bt * betacf(b, a, 1.0_dp - xx) / b
    end if
    p = min(1.0_dp, max(0.0_dp, p))
  end function regularized_beta

  pure elemental real(dp) function student_pdf(x, nu) result(p)
    real(dp), intent(in) :: x, nu
    p = exp(log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu)) / &
        sqrt(nu*pi) * (1.0_dp + x*x/nu)**(-0.5_dp*(nu+1.0_dp))
  end function student_pdf

  pure elemental real(dp) function student_cdf(x, nu) result(p)
    real(dp), intent(in) :: x, nu
    real(dp) :: z, ib
    if (abs(x) <= tiny(1.0_dp)) then
      p = 0.5_dp
      return
    end if
    z = nu / (nu + x*x)
    ib = regularized_beta(z, 0.5_dp*nu, 0.5_dp)
    if (x > 0.0_dp) then
      p = 1.0_dp - 0.5_dp * ib
    else
      p = 0.5_dp * ib
    end if
  end function student_cdf

  real(dp) function student_quantile(p, nu) result(x)
    real(dp), intent(in) :: p, nu
    real(dp) :: lo, hi, mid, pp
    integer :: iter
    pp = clamp_prob(p)
    lo = -8.0_dp
    hi = 8.0_dp
    do while (student_cdf(lo, nu) > pp)
      lo = 2.0_dp * lo
    end do
    do while (student_cdf(hi, nu) < pp)
      hi = 2.0_dp * hi
    end do
    do iter = 1, 100
      mid = 0.5_dp * (lo + hi)
      if (student_cdf(mid, nu) < pp) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function student_quantile

  subroutine gauss_legendre_rule(nn_in, x, w)
    integer, intent(in) :: nn_in
    real(dp), allocatable, intent(out) :: x(:), w(:)
    integer :: nn, i, j, m
    real(dp) :: z, z1, p1, p2, p3, pp
    nn = max(8, nn_in)
    if (mod(nn, 2) /= 0) nn = nn + 1
    allocate(x(nn), w(nn))
    m = (nn + 1) / 2
    do i = 1, m
      z = cos(pi * (real(i,dp) - 0.25_dp) / (real(nn,dp) + 0.5_dp))
      do
        p1 = 1.0_dp
        p2 = 0.0_dp
        do j = 1, nn
          p3 = p2
          p2 = p1
          p1 = ((2.0_dp*real(j,dp)-1.0_dp)*z*p2 - &
                (real(j,dp)-1.0_dp)*p3) / real(j,dp)
        end do
        pp = real(nn,dp) * (z*p1 - p2) / (z*z - 1.0_dp)
        z1 = z
        z = z1 - p1 / pp
        if (abs(z-z1) <= 2.0e-15_dp) exit
      end do
      x(i) = -z
      x(nn+1-i) = z
      w(i) = 2.0_dp / ((1.0_dp-z*z)*pp*pp)
      w(nn+1-i) = w(i)
    end do
  end subroutine gauss_legendre_rule

  real(dp) function integrate_interval(f, a, b, n) result(ans)
    procedure(scalar_fun) :: f
    real(dp), intent(in) :: a, b
    integer, intent(in), optional :: n
    integer :: nn, i
    real(dp) :: xm, xl
    real(dp), allocatable :: x(:), w(:)
    nn = 64
    if (present(n)) nn = n
    call gauss_legendre_rule(nn, x, w)
    xm = 0.5_dp * (b + a)
    xl = 0.5_dp * (b - a)
    ans = 0.0_dp
    do i = 1, size(x)
      ans = ans + w(i) * f(xm + xl*x(i))
    end do
    ans = ans * xl
  end function integrate_interval

  real(dp) function integrate_unit(f, n) result(ans)
    procedure(scalar_fun) :: f
    integer, intent(in), optional :: n
    if (present(n)) then
      ans = integrate_interval(f, 0.0_dp, 1.0_dp, n)
    else
      ans = integrate_interval(f, 0.0_dp, 1.0_dp)
    end if
  end function integrate_unit

  subroutine seed_rng(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i = 1, n
      put(i) = modulo(seed + 104729*i + 7919*i*i, huge(1)-1)
      if (put(i) <= 0) put(i) = i
    end do
    call random_seed(put=put)
  end subroutine seed_rng

  real(dp) function random_uniform() result(u)
    call random_number(u)
    u = clamp_prob(u)
  end function random_uniform

  real(dp) function random_normal() result(z)
    real(dp) :: u1, u2
    u1 = random_uniform()
    u2 = random_uniform()
    z = sqrt(-2.0_dp*log(u1)) * cos(2.0_dp*pi*u2)
  end function random_normal

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key
    do i = 2, size(x)
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j - 1
      end do
      x(j+1) = key
    end do
  end subroutine sort_real

end module rvine_math
