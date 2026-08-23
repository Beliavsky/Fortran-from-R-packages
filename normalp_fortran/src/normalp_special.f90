module normalp_special
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  public :: reg_gamma_p, gamma_quantile
contains
  pure function reg_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    real(dp) :: p
    real(dp) :: ap, del, sumv, b, c, d, h, an, gln
    integer :: n
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = tiny(1.0_dp)/eps
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x == 0.0_dp) then
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
    p = max(0.0_dp, min(1.0_dp, p))
  end function reg_gamma_p

  function gamma_quantile(prob, shape, scale) result(x)
    real(dp), intent(in) :: prob, shape, scale
    real(dp) :: x, lo, hi, mid
    integer :: i
    if (prob <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (prob >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    lo = 0.0_dp
    hi = max(scale*shape, scale)
    do while (reg_gamma_p(shape, hi/scale) < prob)
      hi = 2.0_dp*hi
      if (hi > huge(1.0_dp)/4.0_dp) exit
    end do
    do i = 1, 120
      mid = 0.5_dp*(lo+hi)
      if (reg_gamma_p(shape, mid/scale) < prob) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo+hi)
  end function gamma_quantile
end module normalp_special
