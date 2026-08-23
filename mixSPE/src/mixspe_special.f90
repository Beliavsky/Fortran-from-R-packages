module mixspe_special
  use mvtnorm_kinds, only : dp, pi
  implicit none
  private
  public :: normal_cdf, normal_logcdf, digamma, trigamma, logsumexp
contains
  pure real(dp) function normal_cdf(x) result(v)
    real(dp), intent(in) :: x
    v = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function

  pure real(dp) function normal_logcdf(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: c
    if (x > -8.0_dp) then
      c = normal_cdf(x)
      v = log(max(c, tiny(1.0_dp)))
    else
      v = -0.5_dp*x*x - log(-x) - 0.5_dp*log(2.0_dp*pi)
      v = v + log(max(1.0_dp - 1.0_dp/(x*x) + 3.0_dp/(x**4), 0.5_dp))
    end if
  end function

  pure real(dp) function digamma(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: z, r
    z=x; y=0.0_dp
    do while (z < 8.0_dp)
      y=y-1.0_dp/z; z=z+1.0_dp
    end do
    r=1.0_dp/z
    y=y+log(z)-0.5_dp*r-r*r*(1.0_dp/12.0_dp-r*r*(1.0_dp/120.0_dp-r*r/252.0_dp))
  end function

  pure real(dp) function trigamma(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: z, r, r2
    z=x; y=0.0_dp
    do while (z < 8.0_dp)
      y=y+1.0_dp/(z*z); z=z+1.0_dp
    end do
    r=1.0_dp/z; r2=r*r
    y=y+r+0.5_dp*r2+r*r2/6.0_dp-r*r2*r2/30.0_dp+r*r2*r2*r2/42.0_dp
  end function

  pure real(dp) function logsumexp(x) result(v)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    m=maxval(x)
    if (m <= -huge(1.0_dp)/2) then
      v=m
    else
      v=m+log(sum(exp(x-m)))
    end if
  end function
end module
