module pgnorm_special
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  public :: regularized_gamma_p
contains
  pure function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: p
    real(dp) :: ap, del, sumv, b, c, d, h, an, gln
    integer :: n
    real(dp), parameter :: eps = 1.0e-14_dp, fpmin = 1.0e-300_dp

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
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, 10000
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x + a*log(x) - gln)
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/b
      h = d
      do n = 1, 10000
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
      p = 1.0_dp - exp(-x + a*log(x) - gln)*h
    end if
    p = min(1.0_dp, max(0.0_dp, p))
  end function regularized_gamma_p
end module pgnorm_special
