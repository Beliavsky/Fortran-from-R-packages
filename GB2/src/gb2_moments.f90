! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_moments
  use gb2_kinds, only : dp
  use gb2_special, only : digamma_fn, trigamma_fn, polygamma2_fn, polygamma3_fn, quiet_nan
  use gb2_distribution, only : pgb2
  implicit none
  private
  public :: moment_gb2, incomplete_moment_gb2, expected_log_gb2, variance_log_gb2
  public :: skewness_log_gb2, kurtosis_log_gb2
contains
  pure real(dp) function moment_gb2(k,shape1,scale,shape2,shape3) result(v)
    real(dp), intent(in) :: k,shape1,scale,shape2,shape3
    real(dp) :: pk,qk,lv
    if(shape1<=0.0_dp .or. scale<=0.0_dp .or. shape2<=0.0_dp .or. shape3<=0.0_dp) then
    v=quiet_nan()
    return
    end if
    pk=shape2+k/shape1
    qk=shape3-k/shape1
    if(pk<=0.0_dp .or. qk<=0.0_dp) then
    v=quiet_nan()
    return
    end if
    lv=k*log(scale)+log_gamma(pk)+log_gamma(qk)-log_gamma(shape2)-log_gamma(shape3)
    v=exp(lv)
  end function moment_gb2

  pure real(dp) function incomplete_moment_gb2(x,k,shape1,scale,shape2,shape3) result(v)
    real(dp), intent(in) :: x,k,shape1,scale,shape2,shape3
    real(dp) :: pk,qk
    pk=shape2+k/shape1
    qk=shape3-k/shape1
    if(pk<=0.0_dp .or. qk<=0.0_dp) then
    v=quiet_nan()
    return
    end if
    v=pgb2(x,shape1,scale,pk,qk)
  end function incomplete_moment_gb2

  pure real(dp) function expected_log_gb2(shape1,scale,shape2,shape3) result(v)
    real(dp), intent(in) :: shape1,scale,shape2,shape3
    v=log(scale)+(digamma_fn(shape2)-digamma_fn(shape3))/shape1
  end function expected_log_gb2

  pure real(dp) function variance_log_gb2(shape1,shape2,shape3) result(v)
    real(dp), intent(in) :: shape1,shape2,shape3
    v=(trigamma_fn(shape2)+trigamma_fn(shape3))/shape1**2
  end function variance_log_gb2

  pure real(dp) function skewness_log_gb2(shape2,shape3) result(v)
    real(dp), intent(in) :: shape2,shape3
    v=(polygamma2_fn(shape2)-polygamma2_fn(shape3))/(trigamma_fn(shape2)+trigamma_fn(shape3))**1.5_dp
  end function skewness_log_gb2

  pure real(dp) function kurtosis_log_gb2(shape2,shape3) result(v)
    real(dp), intent(in) :: shape2,shape3
    v=(polygamma3_fn(shape2)+polygamma3_fn(shape3))/(trigamma_fn(shape2)+trigamma_fn(shape3))**2
  end function kurtosis_log_gb2
end module gb2_moments
