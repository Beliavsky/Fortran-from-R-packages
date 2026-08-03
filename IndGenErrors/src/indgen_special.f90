! SPDX-License-Identifier: GPL-3.0-only
module indgen_special
  use indgen_kinds, only : dp
  implicit none
  private

  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
  real(dp), parameter :: sqrt_two = 1.4142135623730950488016887242097_dp

  public :: normal_cdf, normal_quantile, chi_square_survival

contains

  elemental function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    real(dp) :: p
    p = 0.5_dp*erfc(-x/sqrt_two)
  end function normal_cdf

  elemental function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: x, q, r
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
    real(dp), parameter :: phigh = 1.0_dp-plow

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
        (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
        ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
  end function normal_quantile

  recursive function log_gamma_lanczos(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v, y, t, s
    real(dp), parameter :: cof(9) = [ &
      0.99999999999980993_dp, 676.5203681218851_dp, &
      -1259.1392167224028_dp, 771.32342877765313_dp, &
      -176.61502916214059_dp, 12.507343278686905_dp, &
      -0.13857109526572012_dp, 9.9843695780195716e-6_dp, &
      1.5056327351493116e-7_dp ]
    integer :: i

    if (x < 0.5_dp) then
      v = log(pi)-log(sin(pi*x))-log_gamma_lanczos(1.0_dp-x)
      return
    end if
    y = x-1.0_dp
    s = cof(1)
    do i = 2, size(cof)
      s = s+cof(i)/(y+real(i-1,dp))
    end do
    t = y+7.5_dp
    v = 0.5_dp*log(2.0_dp*pi)+(y+0.5_dp)*log(t)-t+log(s)
  end function log_gamma_lanczos

  function gamma_p(a,x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: p, ap, del, sumv, b, c, d, h, an
    integer :: n
    real(dp), parameter :: eps = 2.0e-15_dp
    real(dp), parameter :: tiny_value = 1.0e-300_dp

    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x < a+1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, 10000
        ap = ap+1.0_dp
        del = del*x/ap
        sumv = sumv+del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x+a*log(x)-log_gamma_lanczos(a))
    else
      b = x+1.0_dp-a
      c = 1.0_dp/tiny_value
      d = 1.0_dp/b
      h = d
      do n = 1, 10000
        an = -real(n,dp)*(real(n,dp)-a)
        b = b+2.0_dp
        d = an*d+b
        if (abs(d) < tiny_value) d = tiny_value
        c = b+an/c
        if (abs(c) < tiny_value) c = tiny_value
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      p = 1.0_dp-exp(-x+a*log(x)-log_gamma_lanczos(a))*h
    end if
    p = min(max(p,0.0_dp),1.0_dp)
  end function gamma_p

  function chi_square_survival(x,df) result(q)
    real(dp), intent(in) :: x, df
    real(dp) :: q
    if (x <= 0.0_dp) then
      q = 1.0_dp
    else if (df <= 0.0_dp) then
      q = 0.0_dp
    else
      q = 1.0_dp-gamma_p(0.5_dp*df,0.5_dp*x)
      q = min(max(q,0.0_dp),1.0_dp)
    end if
  end function chi_square_survival

end module indgen_special
