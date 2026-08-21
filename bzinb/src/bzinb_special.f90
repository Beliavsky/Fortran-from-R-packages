module bzinb_special
  use bzinb_kinds, only : dp
  implicit none
  private
  public :: digamma_fn, trigamma_fn, inverse_digamma, logsumexp_pair
contains
  pure real(dp) function digamma_fn(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: y, inv, inv2
    v = 0.0_dp
    y = x
    if (y <= 0.0_dp) then
      v = huge(1.0_dp)
      return
    end if
    do while (y < 8.0_dp)
      v = v - 1.0_dp/y
      y = y + 1.0_dp
    end do
    inv = 1.0_dp/y
    inv2 = inv*inv
    v = v + log(y) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp - &
        inv2*(1.0_dp/252.0_dp - inv2*(1.0_dp/240.0_dp))))
  end function digamma_fn

  pure real(dp) function trigamma_fn(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: y, inv, inv2
    v = 0.0_dp
    y = x
    if (y <= 0.0_dp) then
      v = huge(1.0_dp)
      return
    end if
    do while (y < 8.0_dp)
      v = v + 1.0_dp/(y*y)
      y = y + 1.0_dp
    end do
    inv = 1.0_dp/y
    inv2 = inv*inv
    v = v + inv + 0.5_dp*inv2 + inv*inv2/6.0_dp - inv*inv2*inv2/30.0_dp + &
        inv*inv2*inv2*inv2/42.0_dp - inv*inv2*inv2*inv2*inv2/30.0_dp
  end function trigamma_fn

  pure real(dp) function inverse_digamma(y, x0) result(x)
    real(dp), intent(in) :: y
    real(dp), intent(in), optional :: x0
    real(dp) :: h, xx
    integer :: it
    if (present(x0)) then
      xx = max(x0, 1.0e-10_dp)
    else if (y >= -2.22_dp) then
      xx = exp(y) + 0.5_dp
    else
      xx = -1.0_dp/(y + 0.5772156649015328606_dp)
    end if
    do it = 1, 100
      h = (digamma_fn(xx) - y)/trigamma_fn(xx)
      do while (h >= xx)
        h = 0.5_dp*h
      end do
      xx = xx - h
      if (abs(h) < 1.0e-12_dp*(1.0_dp + xx)) exit
    end do
    x = xx
  end function inverse_digamma

  pure real(dp) function logsumexp_pair(a, b) result(v)
    real(dp), intent(in) :: a, b
    real(dp) :: m
    if (a <= -huge(1.0_dp)/4.0_dp) then
      v = b
    else if (b <= -huge(1.0_dp)/4.0_dp) then
      v = a
    else
      m = max(a,b)
      v = m + log(exp(a-m) + exp(b-m))
    end if
  end function logsumexp_pair
end module bzinb_special
