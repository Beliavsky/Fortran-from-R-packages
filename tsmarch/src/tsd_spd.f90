! SPDX-License-Identifier: GPL-2.0-only
module tsd_spd
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use ghyp_kinds, only : dp
  use ghyp_rng, only : rng_state, uniform_rng
  use tsd_types, only : tsd_success, tsd_invalid_argument
  use tsd_math, only : sample_sd, quantile_type7, sort_real, normal_pdf, normal_cdf
  implicit none
  private

  type, public :: spd_specification
    real(dp), allocatable :: y(:)
    real(dp) :: lower_probability = 0.1_dp
    real(dp) :: upper_probability = 0.9_dp
    character(len=12) :: kernel_type = 'normal'
  end type spd_specification

  type, public :: spd_fit
    real(dp), allocatable :: y(:)
    real(dp) :: lower_probability = 0.1_dp
    real(dp) :: upper_probability = 0.9_dp
    real(dp) :: lower_threshold = 0.0_dp
    real(dp) :: upper_threshold = 0.0_dp
    real(dp) :: lower_scale = 1.0_dp
    real(dp) :: lower_shape = 0.0_dp
    real(dp) :: upper_scale = 1.0_dp
    real(dp) :: upper_shape = 0.0_dp
    real(dp) :: bandwidth = 1.0_dp
    character(len=12) :: kernel_type = 'normal'
    real(dp) :: log_likelihood = 0.0_dp
    real(dp) :: aic = 0.0_dp
    real(dp) :: bic = 0.0_dp
    real(dp) :: covariance(4,4) = 0.0_dp
    integer :: status = tsd_invalid_argument
    character(len=160) :: message = ''
  end type spd_fit

  public :: spd_modelspec, estimate_spd, dspd, pspd, qspd, rspd
  public :: gpd_density, gpd_cdf, gpd_quantile

contains

  function spd_modelspec(y,lower,upper,kernel_type) result(spec)
    real(dp),intent(in)::y(:)
    real(dp),intent(in),optional::lower,upper
    character(len=*),intent(in),optional::kernel_type
    type(spd_specification)::spec
    allocate(spec%y(size(y)))
    spec%y=y
    if(present(lower))spec%lower_probability=lower
    if(present(upper))spec%upper_probability=upper
    if(present(kernel_type))spec%kernel_type=adjustl(kernel_type)
  end function spd_modelspec

  function estimate_spd(spec) result(fit)
    type(spd_specification),intent(in)::spec
    type(spd_fit)::fit
    real(dp),allocatable::ys(:),lo_exc(:),up_exc(:)
    real(dp)::s,cov_lower(2,2),cov_upper(2,2)
    integer::n,i,nl,nu
    n=size(spec%y)
    if(n<10 .or. spec%lower_probability<=0.0_dp .or. spec%upper_probability>=1.0_dp .or. &
       spec%lower_probability>=spec%upper_probability .or. any(.not.ieee_is_finite(spec%y)))then
      fit%message='invalid data or tail probabilities'
      return
    end if
    ys=spec%y
    call sort_real(ys)
    allocate(fit%y(size(ys)))
    fit%y=ys
    fit%lower_probability=spec%lower_probability
    fit%upper_probability=spec%upper_probability
    fit%kernel_type=spec%kernel_type
    fit%lower_threshold=quantile_type7(ys,spec%lower_probability)
    fit%upper_threshold=quantile_type7(ys,spec%upper_probability)
    nl=count(ys<fit%lower_threshold)
    nu=count(ys>fit%upper_threshold)
    if(nl<3 .or. nu<3)then
    fit%message='too few tail exceedances'
    return
    end if
    allocate(lo_exc(nl),up_exc(nu))
    nl=0
    nu=0
    do i=1,n
      if(ys(i)<fit%lower_threshold)then
      nl=nl+1
      lo_exc(nl)=fit%lower_threshold-ys(i)
      end if
      if(ys(i)>fit%upper_threshold)then
      nu=nu+1
      up_exc(nu)=ys(i)-fit%upper_threshold
      end if
    end do
    call pwm_gpd(lo_exc,fit%lower_scale,fit%lower_shape,cov_lower)
    call pwm_gpd(up_exc,fit%upper_scale,fit%upper_shape,cov_upper)
    fit%covariance(1:2,1:2)=cov_lower
    fit%covariance(3:4,3:4)=cov_upper
    s=sample_sd(ys)
    fit%bandwidth=1.06_dp*s*real(n,dp)**(-0.2_dp)
    if(.not.ieee_is_finite(fit%bandwidth).or.fit%bandwidth<=0.0_dp)fit%bandwidth=max(s,1.0_dp)*0.1_dp
    fit%log_likelihood=0.0_dp
    do i=1,n
    fit%log_likelihood=fit%log_likelihood+log(max(dspd(ys(i),fit),tiny(1.0_dp)))
    end do
    fit%aic=-2.0_dp*fit%log_likelihood+8.0_dp
    fit%bic=-2.0_dp*fit%log_likelihood+4.0_dp*log(real(n,dp))
    fit%status=tsd_success
    fit%message='PWM tails and kernel interior fitted'
  end function estimate_spd

  real(dp) function dspd(x,fit,log_density) result(value)
    real(dp),intent(in)::x
    type(spd_fit),intent(in)::fit
    logical,intent(in),optional::log_density
    logical::lg
    real(dp)::mass,denom
    lg=.false.
    if(present(log_density))lg=log_density
    if(fit%status/=tsd_success)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    if(x<fit%lower_threshold)then
      value=fit%lower_probability*gpd_density(fit%lower_threshold-x,fit%lower_scale,fit%lower_shape)
    else if(x>fit%upper_threshold)then
      value=(1.0_dp-fit%upper_probability)*gpd_density(x-fit%upper_threshold,fit%upper_scale,fit%upper_shape)
    else
      mass=fit%upper_probability-fit%lower_probability
      denom=kde_cdf(fit%upper_threshold,fit)-kde_cdf(fit%lower_threshold,fit)
      value=mass*kde_density(x,fit)/max(denom,tiny(1.0_dp))
    end if
    if(lg)value=log(max(value,tiny(1.0_dp)))
  end function dspd

  real(dp) function pspd(q,fit,lower_tail) result(value)
    real(dp),intent(in)::q
    type(spd_fit),intent(in)::fit
    logical,intent(in),optional::lower_tail
    logical::lower
    real(dp)::denom
    lower=.true.
    if(present(lower_tail))lower=lower_tail
    if(fit%status/=tsd_success)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    if(q<fit%lower_threshold)then
      value=fit%lower_probability*(1.0_dp-gpd_cdf(fit%lower_threshold-q,fit%lower_scale,fit%lower_shape))
    else if(q>fit%upper_threshold)then
      value=fit%upper_probability+(1.0_dp-fit%upper_probability)*gpd_cdf(q-fit%upper_threshold,fit%upper_scale,fit%upper_shape)
    else
      denom=kde_cdf(fit%upper_threshold,fit)-kde_cdf(fit%lower_threshold,fit)
      value=fit%lower_probability+(fit%upper_probability-fit%lower_probability)* &
        (kde_cdf(q,fit)-kde_cdf(fit%lower_threshold,fit))/max(denom,tiny(1.0_dp))
    end if
    value=min(max(value,0.0_dp),1.0_dp)
    if(.not.lower)value=1.0_dp-value
  end function pspd

  real(dp) function qspd(p,fit,lower_tail) result(value)
    real(dp),intent(in)::p
    type(spd_fit),intent(in)::fit
    logical,intent(in),optional::lower_tail
    real(dp)::pr,lo,hi,mid
    logical::lower
    integer::iter
    lower=.true.
    if(present(lower_tail))lower=lower_tail
    pr=p
    if(.not.lower)pr=1.0_dp-pr
    if(pr<0.0_dp.or.pr>1.0_dp.or.fit%status/=tsd_success)then
    value=ieee_value(0.0_dp,ieee_quiet_nan)
    return
    end if
    if(pr<fit%lower_probability)then
      value=fit%lower_threshold-gpd_quantile(1.0_dp-pr/fit%lower_probability,fit%lower_scale,fit%lower_shape)
      return
    else if(pr>fit%upper_probability)then
      value=fit%upper_threshold+gpd_quantile( &
        (pr-fit%upper_probability)/(1.0_dp-fit%upper_probability), &
        fit%upper_scale,fit%upper_shape)
      return
    end if
    lo=fit%lower_threshold
    hi=fit%upper_threshold
    do iter=1,100
      mid=0.5_dp*(lo+hi)
      if(pspd(mid,fit)<pr)then
      lo=mid
      else
      hi=mid
      end if
      if(abs(hi-lo)<=1.0e-10_dp*(1.0_dp+abs(mid)))exit
    end do
    value=0.5_dp*(lo+hi)
  end function qspd

  function rspd(n,fit,rng) result(x)
    integer,intent(in)::n
    type(spd_fit),intent(in)::fit
    type(rng_state),intent(inout)::rng
    real(dp),allocatable::x(:)
    integer::i
    allocate(x(max(n,0)))
    do i=1,n
    x(i)=qspd(uniform_rng(rng),fit)
    end do
  end function rspd

  pure real(dp) function gpd_density(x,scale,shape) result(value)
    real(dp),intent(in)::x,scale,shape
    real(dp)::z
    if(x<0.0_dp.or.scale<=0.0_dp)then
    value=0.0_dp
    return
    end if
    if(abs(shape)<1.0e-10_dp)then
    value=exp(-x/scale)/scale
    return
    end if
    z=1.0_dp+shape*x/scale
    if(z<=0.0_dp)then
    value=0.0_dp
    else
    value=z**(-1.0_dp/shape-1.0_dp)/scale
    end if
  end function gpd_density

  pure real(dp) function gpd_cdf(x,scale,shape) result(value)
    real(dp),intent(in)::x,scale,shape
    real(dp)::z
    if(x<=0.0_dp)then
    value=0.0_dp
    return
    end if
    if(abs(shape)<1.0e-10_dp)then
    value=1.0_dp-exp(-x/scale)
    return
    end if
    z=1.0_dp+shape*x/scale
    if(z<=0.0_dp)then
    value=1.0_dp
    else
    value=1.0_dp-z**(-1.0_dp/shape)
    end if
  end function gpd_cdf

  pure real(dp) function gpd_quantile(p,scale,shape) result(value)
    real(dp),intent(in)::p,scale,shape
    real(dp)::pr
    pr=min(max(p,0.0_dp),1.0_dp-epsilon(1.0_dp))
    if(abs(shape)<1.0e-10_dp)then
    value=-scale*log(1.0_dp-pr)
    else
    value=scale*((1.0_dp-pr)**(-shape)-1.0_dp)/shape
    end if
  end function gpd_quantile

  subroutine pwm_gpd(excess,scale,shape,covariance)
    real(dp),intent(in)::excess(:)
    real(dp),intent(out)::scale,shape,covariance(2,2)
    real(dp),allocatable::x(:)
    real(dp)::mu,a1,pv,denom,oneone,twotwo,cdiag
    integer::i,n
    n=size(excess)
    x=excess
    call sort_real(x)
    mu=sum(x)/real(n,dp)
    a1=0.0_dp
    do i=1,n
    pv=(real(i,dp)-0.35_dp)/real(n,dp)
    a1=a1+x(i)*(1.0_dp-pv)
    end do
    a1=a1/real(n,dp)
    if(abs(mu-2.0_dp*a1)<=tiny(1.0_dp))then
    shape=0.0_dp
    scale=max(mu,tiny(1.0_dp))
    else
      shape=2.0_dp-mu/(mu-2.0_dp*a1)
      scale=(2.0_dp*mu*a1)/(mu-2.0_dp*a1)
    end if
    shape=min(max(shape,-0.49_dp),0.49_dp)
    scale=max(scale,1.0e-12_dp)
    covariance=0.0_dp
    denom=real(n,dp)*(1.0_dp-2.0_dp*shape)*(3.0_dp-2.0_dp*shape)
    if(denom>0.0_dp)then
      oneone=(7.0_dp-18.0_dp*shape+11.0_dp*shape**2-2.0_dp*shape**3)*scale**2/denom
      twotwo=(1.0_dp-shape)*(1.0_dp-shape+2.0_dp*shape**2)*(2.0_dp-shape)**2/denom
      cdiag=scale*(2.0_dp-shape)*(2.0_dp-6.0_dp*shape+7.0_dp*shape**2-2.0_dp*shape**3)/denom
      covariance=reshape([oneone,cdiag,cdiag,twotwo],[2,2])
    end if
  end subroutine pwm_gpd

  real(dp) function kde_density(x,fit) result(value)
    real(dp),intent(in)::x
    type(spd_fit),intent(in)::fit
    integer::i,n
    real(dp)::u
    n=size(fit%y)
    value=0.0_dp
    do i=1,n
    u=(x-fit%y(i))/fit%bandwidth
    value=value+kernel_pdf(u,fit%kernel_type)
    end do
    value=value/(real(n,dp)*fit%bandwidth)
  end function kde_density

  real(dp) function kde_cdf(x,fit) result(value)
    real(dp),intent(in)::x
    type(spd_fit),intent(in)::fit
    integer::i,n
    real(dp)::u
    n=size(fit%y)
    value=0.0_dp
    do i=1,n
    u=(x-fit%y(i))/fit%bandwidth
    value=value+kernel_cdf(u,fit%kernel_type)
    end do
    value=value/real(n,dp)
  end function kde_cdf

  pure real(dp) function kernel_pdf(u,name) result(value)
    real(dp),intent(in)::u
    character(len=*),intent(in)::name
    select case(trim(adjustl(name)))
    case('box')
    if(abs(u)<=1.0_dp)then
    value=0.5_dp
    else
    value=0.0_dp
    end if
    case('epanech','epanechnikov')
    if(abs(u)<=1.0_dp)then
    value=0.75_dp*(1.0_dp-u*u)
    else
    value=0.0_dp
    end if
    case('biweight')
    if(abs(u)<=1.0_dp)then
    value=15.0_dp/16.0_dp*(1.0_dp-u*u)**2
    else
    value=0.0_dp
    end if
    case('triweight')
    if(abs(u)<=1.0_dp)then
    value=35.0_dp/32.0_dp*(1.0_dp-u*u)**3
    else
    value=0.0_dp
    end if
    case default
    value=normal_pdf(u)
    end select
  end function kernel_pdf

  pure real(dp) function kernel_cdf(u,name) result(value)
    real(dp),intent(in)::u
    character(len=*),intent(in)::name
    if(trim(adjustl(name))=='normal')then
    value=normal_cdf(u)
    return
    end if
    if(u<=-1.0_dp)then
    value=0.0_dp
    return
    else if(u>=1.0_dp)then
    value=1.0_dp
    return
    end if
    select case(trim(adjustl(name)))
    case('box')
    value=0.5_dp*(u+1.0_dp)
    case('epanech','epanechnikov')
    value=0.5_dp+0.75_dp*(u-u**3/3.0_dp)
    case('biweight')
    value=0.5_dp+15.0_dp/16.0_dp*(u-2.0_dp*u**3/3.0_dp+u**5/5.0_dp)
    case('triweight')
    value=0.5_dp+35.0_dp/32.0_dp*(u-u**3+3.0_dp*u**5/5.0_dp-u**7/7.0_dp)
    case default
    value=normal_cdf(u)
    end select
  end function kernel_cdf

end module tsd_spd
