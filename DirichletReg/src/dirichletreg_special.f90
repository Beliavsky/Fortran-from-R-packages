! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_special
  use dirichletreg_kinds, only : dp
  implicit none
  private
  public :: digamma, trigamma, normal_cdf, normal_quantile, gamma_q, chi_square_sf

contains

  pure real(dp) function digamma(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: y, r, inv, inv2

    if (x <= 0.0_dp) then
      res = huge(1.0_dp)
      return
    end if

    y = x
    r = 0.0_dp
    do while (y < 8.0_dp)
      r = r - 1.0_dp/y
      y = y + 1.0_dp
    end do

    inv = 1.0_dp/y
    inv2 = inv*inv
    res = r + log(y) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp - &
          inv2*(1.0_dp/252.0_dp - inv2*(1.0_dp/240.0_dp - inv2*(5.0_dp/660.0_dp)))))
  end function digamma


  pure real(dp) function trigamma(x) result(res)
    real(dp), intent(in) :: x
    real(dp) :: y, r, inv, inv2

    if (x <= 0.0_dp) then
      res = huge(1.0_dp)
      return
    end if

    y = x
    r = 0.0_dp
    do while (y < 8.0_dp)
      r = r + 1.0_dp/(y*y)
      y = y + 1.0_dp
    end do

    inv = 1.0_dp/y
    inv2 = inv*inv
    res = r + inv + 0.5_dp*inv2 + inv*inv2*(1.0_dp/6.0_dp - inv2*(1.0_dp/30.0_dp - &
          inv2*(1.0_dp/42.0_dp - inv2*(1.0_dp/30.0_dp - inv2*(5.0_dp/66.0_dp)))))
  end function trigamma


  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf


  pure real(dp) function normal_quantile(p) result(x)
    ! Peter J. Acklam's rational approximation, followed by one Newton step.
    real(dp), intent(in) :: p
    real(dp), parameter :: a1=-3.969683028665376e+01_dp, a2= 2.209460984245205e+02_dp
    real(dp), parameter :: a3=-2.759285104469687e+02_dp, a4= 1.383577518672690e+02_dp
    real(dp), parameter :: a5=-3.066479806614716e+01_dp, a6= 2.506628277459239e+00_dp
    real(dp), parameter :: b1=-5.447609879822406e+01_dp, b2= 1.615858368580409e+02_dp
    real(dp), parameter :: b3=-1.556989798598866e+02_dp, b4= 6.680131188771972e+01_dp
    real(dp), parameter :: b5=-1.328068155288572e+01_dp
    real(dp), parameter :: c1=-7.784894002430293e-03_dp, c2=-3.223964580411365e-01_dp
    real(dp), parameter :: c3=-2.400758277161838e+00_dp, c4=-2.549732539343734e+00_dp
    real(dp), parameter :: c5= 4.374664141464968e+00_dp, c6= 2.938163982698783e+00_dp
    real(dp), parameter :: d1= 7.784695709041462e-03_dp, d2= 3.224671290700398e-01_dp
    real(dp), parameter :: d3= 2.445134137142996e+00_dp, d4= 3.754408661907416e+00_dp
    real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
    real(dp) :: q, r, e, pdf

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if

    if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if

    e = normal_cdf(x) - p
    pdf = exp(-0.5_dp*x*x)/sqrt(2.0_dp*acos(-1.0_dp))
    x = x - e/pdf
  end function normal_quantile


  pure real(dp) function gamma_q(a, x) result(q)
    ! Regularized upper incomplete gamma Q(a,x), Numerical Recipes style.
    real(dp), intent(in) :: a, x
    integer, parameter :: itmax = 10000
    real(dp), parameter :: eps = 5.0e-15_dp, fpmin = tiny(1.0_dp)/eps
    integer :: n
    real(dp) :: ap, del, sumv, b, c, d, h, an, gln, p

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      q = huge(1.0_dp)
      return
    end if
    if (x <= tiny(1.0_dp)) then
      q = 1.0_dp
      return
    end if

    gln = log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, itmax
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x + a*log(x) - gln)
      q = max(0.0_dp, min(1.0_dp, 1.0_dp-p))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do n = 1, itmax
        an = -real(n,dp)*(real(n,dp)-a)
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
      q = h*exp(-x + a*log(x) - gln)
      q = max(0.0_dp, min(1.0_dp, q))
    end if
  end function gamma_q


  pure real(dp) function chi_square_sf(x, df) result(p)
    real(dp), intent(in) :: x, df
    if (x < 0.0_dp .or. df <= 0.0_dp) then
      p = huge(1.0_dp)
    else
      p = gamma_q(0.5_dp*df, 0.5_dp*x)
    end if
  end function chi_square_sf

end module dirichletreg_special
