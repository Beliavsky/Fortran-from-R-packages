! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_distribution
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf
  use gb2_kinds, only : dp
  use gb2_special, only : log_beta, reg_incomplete_beta, beta_quantile, random_beta, quiet_nan
  implicit none
  private
  public :: dgb2, pgb2, qgb2, rgb2
contains
  pure real(dp) function dgb2(x,shape1,scale,shape2,shape3) result(f)
    real(dp), intent(in) :: x,shape1,scale,shape2,shape3
    real(dp) :: lx,t,logf
    if(shape1<=0.0_dp .or. scale<=0.0_dp .or. shape2<=0.0_dp .or. shape3<=0.0_dp .or. x<0.0_dp) then
      f=quiet_nan()
      return
    end if
    if(.not.ieee_is_finite(x)) then
      f=0.0_dp
      return
    end if
    if(x<=0.0_dp) then
      if(shape1*shape2>1.0_dp+20.0_dp*epsilon(1.0_dp)) then
        f=0.0_dp
      else if(abs(shape1*shape2-1.0_dp)<=20.0_dp*epsilon(1.0_dp)) then
        f=shape1/(scale*exp(log_beta(shape2,shape3)))
      else
        f=ieee_value(0.0_dp,ieee_positive_inf)
      end if
      return
    end if
    lx=log(x/scale)
    t=shape1*lx
    if(t>40.0_dp) then
      logf=log(shape1/scale)-log_beta(shape2,shape3)+(shape1*shape2-1.0_dp)*lx-(shape2+shape3)*(t+log(1.0_dp+exp(-t)))
    else
      logf=log(shape1/scale)-log_beta(shape2,shape3)+(shape1*shape2-1.0_dp)*lx-(shape2+shape3)*log(1.0_dp+exp(t))
    end if
    if(logf<log(tiny(1.0_dp))) then
    f=0.0_dp
    else
    f=exp(logf)
    end if
  end function dgb2

  pure real(dp) function pgb2(x,shape1,scale,shape2,shape3) result(p)
    real(dp), intent(in) :: x,shape1,scale,shape2,shape3
    real(dp) :: t,z,e
    if(shape1<=0.0_dp .or. scale<=0.0_dp .or. shape2<=0.0_dp .or. shape3<=0.0_dp .or. x<0.0_dp) then
      p=quiet_nan()
      return
    end if
    if(x<=0.0_dp) then
    p=0.0_dp
    return
    end if
    if(.not.ieee_is_finite(x)) then
    p=1.0_dp
    return
    end if
    t=shape1*log(x/scale)
    if(t>=0.0_dp) then
      e=exp(-t)
      z=1.0_dp/(1.0_dp+e)
    else
      e=exp(t)
      z=e/(1.0_dp+e)
    end if
    p=reg_incomplete_beta(z,shape2,shape3)
  end function pgb2

  real(dp) function qgb2(prob,shape1,scale,shape2,shape3) result(q)
    real(dp), intent(in) :: prob,shape1,scale,shape2,shape3
    real(dp) :: z,ratio
    if(shape1<=0.0_dp .or. scale<=0.0_dp .or. shape2<=0.0_dp .or. shape3<=0.0_dp .or. prob<0.0_dp .or. prob>1.0_dp) then
      q=quiet_nan()
      return
    end if
    if(prob<=0.0_dp) then
    q=0.0_dp
    return
    end if
    if(prob>=1.0_dp) then
    q=huge(1.0_dp)
    return
    end if
    if(prob<=0.5_dp) then
      z=beta_quantile(prob,shape2,shape3)
      ratio=z/(1.0_dp-z)
    else
      z=beta_quantile(1.0_dp-prob,shape3,shape2)
      ratio=(1.0_dp-z)/z
    end if
    q=scale*ratio**(1.0_dp/shape1)
  end function qgb2

  subroutine rgb2(n,shape1,scale,shape2,shape3,x)
    integer, intent(in) :: n
    real(dp), intent(in) :: shape1,scale,shape2,shape3
    real(dp), intent(out) :: x(:)
    real(dp) :: z
    integer :: i
    if(size(x)/=n) error stop 'rgb2: size(x) must equal n'
    if(shape1<=0.0_dp .or. scale<=0.0_dp .or. shape2<=0.0_dp .or. shape3<=0.0_dp) error stop 'rgb2: parameters must be positive'
    do i=1,n
      call random_beta(shape2,shape3,z)
      x(i)=scale*(z/(1.0_dp-z))**(1.0_dp/shape1)
    end do
  end subroutine rgb2
end module gb2_distribution
