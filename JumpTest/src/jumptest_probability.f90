! SPDX-License-Identifier: MIT
module jumptest_probability
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use jumptest_kinds, only : dp
  implicit none
  private

  public :: normal_cdf, normal_quantile, chi_square_cdf, beta_cdf
  public :: sample_quantile_type7

contains

  elemental function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p

    p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  elemental function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x
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
    real(dp), parameter :: plow = 0.02425_dp
    real(dp), parameter :: phigh = 1.0_dp - plow

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6))/ &
        ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6))*q/ &
        (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + 1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp - p))
      x = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6))/ &
        ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0_dp)
    end if

    if (ieee_is_finite(x)) then
      x = x - (normal_cdf(x) - p)/(exp(-0.5_dp*x*x)/sqrt(2.0_dp*acos(-1.0_dp)))
    end if
  end function normal_quantile

  pure function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: p
    real(dp) :: ap, del, sum, b, c, d, h, an
    real(dp) :: gln
    integer :: i
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 5.0e-15_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= tiny(1.0_dp)) then
      p = 0.0_dp
      return
    end if

    gln = log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sum = 1.0_dp/a
      del = sum
      do i = 1, maxit
        ap = ap + 1.0_dp
        del = del*x/ap
        sum = sum + del
        if (abs(del) <= abs(sum)*eps) exit
      end do
      p = sum*exp(-x + a*log(x) - gln)
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/max(abs(b), fpmin)
      if (b < 0.0_dp) d = -d
      h = d
      do i = 1, maxit
        an = -real(i, dp)*(real(i, dp) - a)
        b = b + 2.0_dp
        d = an*d + b
        if (abs(d) < fpmin) d = sign(fpmin, d)
        c = b + an/c
        if (abs(c) < fpmin) c = sign(fpmin, c)
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del - 1.0_dp) <= eps) exit
      end do
      p = 1.0_dp - exp(-x + a*log(x) - gln)*h
    end if
    p = min(1.0_dp, max(0.0_dp, p))
  end function regularized_gamma_p

  elemental function chi_square_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    real(dp) :: p

    if (df <= 0.0_dp .or. x < 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      p = regularized_gamma_p(0.5_dp*df, 0.5_dp*x)
    end if
  end function chi_square_cdf

  pure function beta_continued_fraction(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    real(dp) :: cf
    real(dp) :: qab, qap, qam, c, d, h, aa, del
    integer :: m, m2
    integer, parameter :: maxit = 10000
    real(dp), parameter :: eps = 5.0e-15_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp

    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab*x/qap
    if (abs(d) < fpmin) d = sign(fpmin, d)
    d = 1.0_dp/d
    h = d
    do m = 1, maxit
      m2 = 2*m
      aa = real(m, dp)*(b - real(m, dp))*x/ &
        ((qam + real(m2, dp))*(a + real(m2, dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = sign(fpmin, d)
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = sign(fpmin, c)
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a + real(m, dp))*(qab + real(m, dp))*x/ &
        ((a + real(m2, dp))*(qap + real(m2, dp)))
      d = 1.0_dp + aa*d
      if (abs(d) < fpmin) d = sign(fpmin, d)
      c = 1.0_dp + aa/c
      if (abs(c) < fpmin) c = sign(fpmin, c)
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del - 1.0_dp) <= eps) exit
    end do
    cf = h
  end function beta_continued_fraction

  elemental function beta_cdf(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: p
    real(dp) :: bt

    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x >= 1.0_dp) then
      p = 1.0_dp
      return
    end if

    bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
      a*log(x) + b*log(1.0_dp - x))
    if (x < (a + 1.0_dp)/(a + b + 2.0_dp)) then
      p = bt*beta_continued_fraction(a, b, x)/a
    else
      p = 1.0_dp - bt*beta_continued_fraction(b, a, 1.0_dp - x)/b
    end if
    p = min(1.0_dp, max(0.0_dp, p))
  end function beta_cdf

  function sample_quantile_type7(x, probability) result(q)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: probability
    real(dp) :: q
    real(dp), allocatable :: work(:)
    real(dp) :: h, frac
    integer :: n, lower

    n = size(x)
    if (n < 1 .or. probability < 0.0_dp .or. probability > 1.0_dp) then
      q = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    allocate(work(n))
    work = x
    call sort_real(work)
    if (n == 1) then
      q = work(1)
      return
    end if
    h = 1.0_dp + real(n - 1, dp)*probability
    lower = int(floor(h))
    frac = h - real(lower, dp)
    if (lower >= n) then
      q = work(n)
    else
      q = work(lower) + frac*(work(lower + 1) - work(lower))
    end if
  end function sample_quantile_type7

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    integer :: i, j
    real(dp) :: key

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
  end subroutine sort_real

end module jumptest_probability
