! SPDX-License-Identifier: GPL-2.0-or-later
module tmvtnorm_distributions
  use mvtnorm_kinds, only : dp
  use mvtnorm_types, only : probability_control, probability_result
  use mvtnorm_probabilities, only : mvn_prob => pmvnorm, mvt_prob => pmvt
  use mvtnorm_distributions, only : dmvnorm_one, dmvt_one
  use tmvtnorm_utils, only : in_box
  implicit none
  private
  public :: ptmvnorm, dtmvnorm_one, dtmvnorm
  public :: ptmvt, dtmvt_one, dtmvt

contains

  real(dp) function ptmvnorm(lowerx,upperx,mean,sigma,lower,upper,control) result(p)
    real(dp),intent(in)::lowerx(:),upperx(:),mean(:),sigma(:,:),lower(:),upper(:)
    type(probability_control),intent(in),optional::control
    real(dp),allocatable::lo(:),up(:)
    real(dp)::den
    type(probability_result) :: a,b
    allocate(lo(size(lower)),up(size(upper)))
    lo=max(lowerx,lower)
    up=min(upperx,upper)
    if(any(lo>=up)) then
    p=0.0_dp
    return
    end if
    if(present(control)) then
      a=mvn_prob(lo,up,mean,sigma,control)
      b=mvn_prob(lower,upper,mean,sigma,control)
    else
      a=mvn_prob(lo,up,mean,sigma)
      b=mvn_prob(lower,upper,mean,sigma)
    end if
    den=b%value
    if(den<=0.0_dp) then
    p=0.0_dp
    else
    p=max(0.0_dp,min(1.0_dp,a%value/den))
    end if
  end function ptmvnorm

  real(dp) function dtmvnorm_one(x,mean,sigma,lower,upper,log_density,control) result(v)
    real(dp),intent(in)::x(:),mean(:),sigma(:,:),lower(:),upper(:)
    logical,intent(in),optional::log_density
    type(probability_control),intent(in),optional::control
    logical::ll
    type(probability_result)::pr
    real(dp)::ld
    ll=.false.
    if(present(log_density)) ll=log_density
    if(.not.in_box(x,lower,upper)) then
      v=merge(-huge(1.0_dp),0.0_dp,ll)
      return
    end if
    if(present(control)) then
    pr=mvn_prob(lower,upper,mean,sigma,control)
    else
    pr=mvn_prob(lower,upper,mean,sigma)
    end if
    if(pr%value<=0.0_dp) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    ld=dmvnorm_one(x,mean,sigma,.true.)-log(pr%value)
    v=merge(ld,exp(ld),ll)
  end function dtmvnorm_one

  function dtmvnorm(x,mean,sigma,lower,upper,log_density,control) result(v)
    real(dp),intent(in)::x(:,:),mean(:),sigma(:,:),lower(:),upper(:)
    logical,intent(in),optional::log_density
    type(probability_control),intent(in),optional::control
    real(dp),allocatable::v(:)
    integer::i
    allocate(v(size(x,1)))
    do i=1,size(x,1)
      v(i)=dtmvnorm_one(x(i,:),mean,sigma,lower,upper,log_density,control)
    end do
  end function dtmvnorm

  real(dp) function ptmvt(lowerx,upperx,mean,sigma,df,lower,upper,control) result(p)
    real(dp),intent(in)::lowerx(:),upperx(:),mean(:),sigma(:,:),df,lower(:),upper(:)
    type(probability_control),intent(in),optional::control
    real(dp),allocatable::lo(:),up(:)
    type(probability_result)::a,b
    allocate(lo(size(lower)),up(size(upper)))
    lo=max(lowerx,lower)
    up=min(upperx,upper)
    if(any(lo>=up)) then
    p=0.0_dp
    return
    end if
    if(present(control)) then
      a=mvt_prob(lo,up,mean,sigma,df,control)
      b=mvt_prob(lower,upper,mean,sigma,df,control)
    else
      a=mvt_prob(lo,up,mean,sigma,df)
      b=mvt_prob(lower,upper,mean,sigma,df)
    end if
    if(b%value<=0.0_dp) then
    p=0.0_dp
    else
    p=max(0.0_dp,min(1.0_dp,a%value/b%value))
    end if
  end function ptmvt

  real(dp) function dtmvt_one(x,mean,sigma,df,lower,upper,log_density,control) result(v)
    real(dp),intent(in)::x(:),mean(:),sigma(:,:),df,lower(:),upper(:)
    logical,intent(in),optional::log_density
    type(probability_control),intent(in),optional::control
    logical::ll
    real(dp)::ld
    type(probability_result)::pr
    ll=.false.
    if(present(log_density)) ll=log_density
    if(.not.in_box(x,lower,upper)) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    if(present(control)) then
    pr=mvt_prob(lower,upper,mean,sigma,df,control)
    else
    pr=mvt_prob(lower,upper,mean,sigma,df)
    end if
    if(pr%value<=0.0_dp) then
    v=merge(-huge(1.0_dp),0.0_dp,ll)
    return
    end if
    ld=dmvt_one(x,mean,sigma,df,.true.)-log(pr%value)
    v=merge(ld,exp(ld),ll)
  end function dtmvt_one

  function dtmvt(x,mean,sigma,df,lower,upper,log_density,control) result(v)
    real(dp),intent(in)::x(:,:),mean(:),sigma(:,:),df,lower(:),upper(:)
    logical,intent(in),optional::log_density
    type(probability_control),intent(in),optional::control
    real(dp),allocatable::v(:)
    integer::i
    allocate(v(size(x,1)))
    do i=1,size(x,1)
      v(i)=dtmvt_one(x(i,:),mean,sigma,df,lower,upper,log_density,control)
    end do
  end function dtmvt

end module tmvtnorm_distributions
