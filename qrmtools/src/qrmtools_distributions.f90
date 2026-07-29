! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_distributions
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_positive_inf
  use qrmtools_kinds, only : dp
  use qrmtools_stats, only : random_uniform, normal_pdf, normal_cdf, normal_quantile, &
    student_t_pdf, student_t_cdf, student_t_quantile
  implicit none
  private
  integer, parameter, public :: distribution_normal=1, distribution_student=2
  integer, parameter, public :: distribution_gev=3, distribution_gpd=4, distribution_pareto=5
  type, public :: distribution_component
    integer :: family=distribution_normal
    real(dp) :: parameters(4)=0.0_dp
  end type distribution_component
  public :: dgev,pgev,qgev,rgev,dgpd,pgpd,qgpd,rgpd
  public :: dpar,ppar,qpar,rpar,dgpdtail,pgpdtail,qgpdtail,rgpdtail
  public :: component_pdf,component_cdf,component_quantile
  public :: composite_pdf,composite_cdf,composite_quantile,composite_random
contains
  pure real(dp) function dgev(x,shape,loc,scale,log_density) result(value)
    real(dp), intent(in) :: x,shape
    real(dp), intent(in), optional :: loc,scale
    logical, intent(in), optional :: log_density
    real(dp) :: mu,sigma,y,z,lv
    logical :: lg
    mu=0.0_dp; sigma=1.0_dp; lg=.false.
    if(present(loc))mu=loc; if(present(scale))sigma=scale; if(present(log_density))lg=log_density
    if(sigma<=0.0_dp) then
      lv=ieee_value(1.0_dp,ieee_negative_inf)
    else
      y=(x-mu)/sigma
      if(abs(shape)<=1.0e-12_dp) then
        lv=-log(sigma)-y-exp(-y)
      else
        z=1.0_dp+shape*y
        if(z<=0.0_dp) then
          lv=ieee_value(1.0_dp,ieee_negative_inf)
        else
          lv=-log(sigma)+(-1.0_dp/shape-1.0_dp)*log(z)-z**(-1.0_dp/shape)
        end if
      end if
    end if
    if(lg) then; value=lv; else; value=exp(lv); end if
  end function dgev

  pure real(dp) function pgev(q,shape,loc,scale,lower_tail,log_probability) result(value)
    real(dp), intent(in) :: q,shape
    real(dp), intent(in), optional :: loc,scale
    logical, intent(in), optional :: lower_tail,log_probability
    real(dp) :: mu,sigma,y,z,cdf
    logical :: lower,lp
    mu=0.0_dp; sigma=1.0_dp; lower=.true.; lp=.false.
    if(present(loc))mu=loc; if(present(scale))sigma=scale
    if(present(lower_tail))lower=lower_tail; if(present(log_probability))lp=log_probability
    if(sigma<=0.0_dp) then; value=ieee_value(1.0_dp,ieee_negative_inf); return; end if
    y=(q-mu)/sigma
    if(abs(shape)<=1.0e-12_dp) then
      cdf=exp(-exp(-y))
    else
      z=1.0_dp+shape*y
      if(z<=0.0_dp) then
        if(shape>0.0_dp) then; cdf=0.0_dp; else; cdf=1.0_dp; end if
      else
        cdf=exp(-z**(-1.0_dp/shape))
      end if
    end if
    if(.not.lower)cdf=1.0_dp-cdf
    if(lp) then; value=log(max(cdf,tiny(1.0_dp))); else; value=cdf; end if
  end function pgev

  pure real(dp) function qgev(p,shape,loc,scale,lower_tail,log_probability) result(value)
    real(dp), intent(in) :: p,shape
    real(dp), intent(in), optional :: loc,scale
    logical, intent(in), optional :: lower_tail,log_probability
    real(dp) :: mu,sigma,prob,y
    logical :: lower,lp
    mu=0.0_dp; sigma=1.0_dp; lower=.true.; lp=.false.
    if(present(loc))mu=loc; if(present(scale))sigma=scale
    if(present(lower_tail))lower=lower_tail; if(present(log_probability))lp=log_probability
    prob=p; if(lp)prob=exp(min(p,0.0_dp)); if(.not.lower)prob=1.0_dp-prob
    prob=min(max(prob,0.0_dp),1.0_dp)
    if(prob<=0.0_dp) then
      if(shape<0.0_dp) then; value=mu-sigma/shape; else; value=ieee_value(1.0_dp,ieee_negative_inf); end if
      return
    else if(prob>=1.0_dp) then
      if(shape<0.0_dp) then; value=mu-sigma/shape; else; value=ieee_value(1.0_dp,ieee_positive_inf); end if
      return
    end if
    if(abs(shape)<=1.0e-12_dp) then
      y=-log(-log(prob))
    else
      y=((-log(prob))**(-shape)-1.0_dp)/shape
    end if
    value=mu+sigma*y
  end function qgev

  real(dp) function rgev(shape,loc,scale) result(value)
    real(dp), intent(in) :: shape
    real(dp), intent(in), optional :: loc,scale
    value=qgev(random_uniform(),shape,loc,scale)
  end function rgev

  pure real(dp) function dgpd(x,shape,scale,log_density) result(value)
    real(dp), intent(in) :: x,shape,scale
    logical, intent(in), optional :: log_density
    logical :: lg
    real(dp) :: z,lv
    lg=.false.; if(present(log_density))lg=log_density
    lv=ieee_value(1.0_dp,ieee_negative_inf)
    if(scale>0.0_dp .and. x>=0.0_dp) then
      if(abs(shape)<=1.0e-12_dp) then
        lv=-x/scale-log(scale)
      else
        z=1.0_dp+shape*x/scale
        if(z>0.0_dp) lv=-(1.0_dp/shape+1.0_dp)*log(z)-log(scale)
      end if
    end if
    if(lg) then; value=lv; else; value=exp(lv); end if
  end function dgpd

  pure real(dp) function pgpd(q,shape,scale,lower_tail,log_probability) result(value)
    real(dp), intent(in) :: q,shape,scale
    logical, intent(in), optional :: lower_tail,log_probability
    logical :: lower,lp
    real(dp) :: cdf,z
    lower=.true.; lp=.false.; if(present(lower_tail))lower=lower_tail
    if(present(log_probability))lp=log_probability
    if(q<=0.0_dp) then
      cdf=0.0_dp
    else if(shape<0.0_dp .and. q>=-scale/shape) then
      cdf=1.0_dp
    else if(abs(shape)<=1.0e-12_dp) then
      cdf=1.0_dp-exp(-q/scale)
    else
      z=1.0_dp+shape*q/scale; cdf=1.0_dp-z**(-1.0_dp/shape)
    end if
    if(.not.lower)cdf=1.0_dp-cdf
    if(lp) then; value=log(max(cdf,tiny(1.0_dp))); else; value=cdf; end if
  end function pgpd

  pure real(dp) function qgpd(p,shape,scale,lower_tail,log_probability) result(value)
    real(dp), intent(in) :: p,shape,scale
    logical, intent(in), optional :: lower_tail,log_probability
    logical :: lower,lp
    real(dp) :: prob
    lower=.true.; lp=.false.; if(present(lower_tail))lower=lower_tail
    if(present(log_probability))lp=log_probability
    prob=p; if(lp)prob=exp(min(p,0.0_dp)); if(.not.lower)prob=1.0_dp-prob
    prob=min(max(prob,0.0_dp),1.0_dp)
    if(prob<=0.0_dp) then; value=0.0_dp; return; end if
    if(prob>=1.0_dp) then
      if(shape<0.0_dp) then; value=-scale/shape; else; value=ieee_value(1.0_dp,ieee_positive_inf); end if
      return
    end if
    if(abs(shape)<=1.0e-12_dp) then
      value=-scale*log(1.0_dp-prob)
    else
      value=(scale/shape)*((1.0_dp-prob)**(-shape)-1.0_dp)
    end if
  end function qgpd

  real(dp) function rgpd(shape,scale) result(value)
    real(dp), intent(in) :: shape,scale
    value=qgpd(random_uniform(),shape,scale)
  end function rgpd

  pure real(dp) function dpar(x,shape,scale,log_density) result(value)
    real(dp), intent(in) :: x,shape
    real(dp), intent(in), optional :: scale
    logical, intent(in), optional :: log_density
    real(dp) :: k,lv
    logical :: lg
    k=1.0_dp; lg=.false.; if(present(scale))k=scale; if(present(log_density))lg=log_density
    if(x<0.0_dp .or. shape<=0.0_dp .or. k<=0.0_dp) then
      lv=ieee_value(1.0_dp,ieee_negative_inf)
    else
      lv=log(shape/k)+(shape+1.0_dp)*log(k/(k+x))
    end if
    if(lg) then; value=lv; else; value=exp(lv); end if
  end function dpar

  pure real(dp) function ppar(q,shape,scale,lower_tail,log_probability) result(value)
    real(dp), intent(in) :: q,shape
    real(dp), intent(in), optional :: scale
    logical, intent(in), optional :: lower_tail,log_probability
    real(dp) :: k,cdf
    logical :: lower,lp
    k=1.0_dp; lower=.true.; lp=.false.; if(present(scale))k=scale
    if(present(lower_tail))lower=lower_tail; if(present(log_probability))lp=log_probability
    if(q<=0.0_dp) then; cdf=0.0_dp; else; cdf=1.0_dp-(1.0_dp+q/k)**(-shape); end if
    if(.not.lower)cdf=1.0_dp-cdf
    if(lp) then; value=log(max(cdf,tiny(1.0_dp))); else; value=cdf; end if
  end function ppar

  pure real(dp) function qpar(p,shape,scale,lower_tail,log_probability) result(value)
    real(dp), intent(in) :: p,shape
    real(dp), intent(in), optional :: scale
    logical, intent(in), optional :: lower_tail,log_probability
    real(dp) :: k,prob
    logical :: lower,lp
    k=1.0_dp; lower=.true.; lp=.false.; if(present(scale))k=scale
    if(present(lower_tail))lower=lower_tail; if(present(log_probability))lp=log_probability
    prob=p; if(lp)prob=exp(min(p,0.0_dp)); if(.not.lower)prob=1.0_dp-prob
    if(prob<=0.0_dp) then; value=0.0_dp
    else if(prob>=1.0_dp) then; value=ieee_value(1.0_dp,ieee_positive_inf)
    else; value=k*((1.0_dp-prob)**(-1.0_dp/shape)-1.0_dp); end if
  end function qpar

  real(dp) function rpar(shape,scale) result(value)
    real(dp), intent(in) :: shape
    real(dp), intent(in), optional :: scale
    value=qpar(random_uniform(),shape,scale)
  end function rpar

  pure real(dp) function dgpdtail(x,threshold,p_exceed,shape,scale,log_density) result(value)
    real(dp), intent(in) :: x,threshold,p_exceed,shape,scale
    logical, intent(in), optional :: log_density
    logical :: lg
    real(dp) :: lv
    lg=.false.; if(present(log_density))lg=log_density
    lv=log(max(p_exceed,tiny(1.0_dp)))+dgpd(x-threshold,shape,scale,.true.)
    if(lg) then; value=lv; else; value=exp(lv); end if
  end function dgpdtail

  pure real(dp) function pgpdtail(q,threshold,p_exceed,shape,scale,lower_tail,log_probability) result(value)
    real(dp), intent(in) :: q,threshold,p_exceed,shape,scale
    logical, intent(in), optional :: lower_tail,log_probability
    logical :: lower,lp
    real(dp) :: surv,prob
    lower=.true.; lp=.false.; if(present(lower_tail))lower=lower_tail
    if(present(log_probability))lp=log_probability
    surv=p_exceed*pgpd(q-threshold,shape,scale,.false.)
    if(lower) then; prob=1.0_dp-surv; else; prob=surv; end if
    if(lp) then; value=log(max(prob,tiny(1.0_dp))); else; value=prob; end if
  end function pgpdtail

  pure real(dp) function qgpdtail(p,threshold,p_exceed,shape,scale,lower_tail,log_probability) result(value)
    real(dp), intent(in) :: p,threshold,p_exceed,shape,scale
    logical, intent(in), optional :: lower_tail,log_probability
    logical :: lower,lp
    real(dp) :: prob
    lower=.true.; lp=.false.; if(present(lower_tail))lower=lower_tail
    if(present(log_probability))lp=log_probability
    prob=p; if(lp)prob=exp(min(p,0.0_dp)); if(.not.lower)prob=1.0_dp-prob
    value=threshold+qgpd(1.0_dp-(1.0_dp-prob)/p_exceed,shape,scale)
  end function qgpdtail

  real(dp) function rgpdtail(threshold,p_exceed,shape,scale) result(value)
    real(dp), intent(in) :: threshold,p_exceed,shape,scale
    value=qgpdtail(random_uniform(),threshold,p_exceed,shape,scale)
  end function rgpdtail

  pure real(dp) function component_pdf(x,component) result(value)
    real(dp), intent(in) :: x
    type(distribution_component), intent(in) :: component
    select case(component%family)
    case(distribution_normal)
      value=normal_pdf((x-component%parameters(1))/component%parameters(2))/component%parameters(2)
    case(distribution_student)
      value=student_t_pdf((x-component%parameters(1))/component%parameters(2),component%parameters(3))/component%parameters(2)
    case(distribution_gev)
      value=dgev(x,component%parameters(1),component%parameters(2),component%parameters(3))
    case(distribution_gpd)
      value=dgpd(x,component%parameters(1),component%parameters(2))
    case(distribution_pareto)
      value=dpar(x,component%parameters(1),component%parameters(2))
    case default
      value=0.0_dp
    end select
  end function component_pdf

  pure real(dp) function component_cdf(x,component) result(value)
    real(dp), intent(in) :: x
    type(distribution_component), intent(in) :: component
    select case(component%family)
    case(distribution_normal)
      value=normal_cdf((x-component%parameters(1))/component%parameters(2))
    case(distribution_student)
      value=student_t_cdf((x-component%parameters(1))/component%parameters(2),component%parameters(3))
    case(distribution_gev)
      value=pgev(x,component%parameters(1),component%parameters(2),component%parameters(3))
    case(distribution_gpd)
      value=pgpd(x,component%parameters(1),component%parameters(2))
    case(distribution_pareto)
      value=ppar(x,component%parameters(1),component%parameters(2))
    case default
      value=0.0_dp
    end select
  end function component_cdf

  real(dp) function component_quantile(p,component) result(value)
    real(dp), intent(in) :: p
    type(distribution_component), intent(in) :: component
    select case(component%family)
    case(distribution_normal)
      value=component%parameters(1)+component%parameters(2)*normal_quantile(p)
    case(distribution_student)
      value=component%parameters(1)+component%parameters(2)*student_t_quantile(p,component%parameters(3))
    case(distribution_gev)
      value=qgev(p,component%parameters(1),component%parameters(2),component%parameters(3))
    case(distribution_gpd)
      value=qgpd(p,component%parameters(1),component%parameters(2))
    case(distribution_pareto)
      value=qpar(p,component%parameters(1),component%parameters(2))
    case default
      value=0.0_dp
    end select
  end function component_quantile

  real(dp) function composite_pdf(x,cuts,components,weights) result(value)
    real(dp), intent(in) :: x,cuts(:),weights(:)
    type(distribution_component), intent(in) :: components(:)
    integer :: i
    real(dp) :: low,high,mass
    i=bucket_index(x,cuts)
    low=0.0_dp
    if(i>1) low=component_cdf(cuts(i-1),components(i))
    high=1.0_dp
    if(i<=size(cuts)) high=component_cdf(cuts(i),components(i))
    mass=high-low
    value=weights(i)*component_pdf(x,components(i))/mass
  end function composite_pdf

  real(dp) function composite_cdf(x,cuts,components,weights) result(value)
    real(dp), intent(in) :: x,cuts(:),weights(:)
    type(distribution_component), intent(in) :: components(:)
    integer :: i
    real(dp) :: low,high,mass
    i=bucket_index(x,cuts)
    low=0.0_dp; if(i>1)low=component_cdf(cuts(i-1),components(i))
    high=1.0_dp; if(i<=size(cuts))high=component_cdf(cuts(i),components(i))
    mass=high-low
    value=sum(weights(1:i-1))+weights(i)*(component_cdf(x,components(i))-low)/mass
  end function composite_cdf

  real(dp) function composite_quantile(p,cuts,components,weights) result(value)
    real(dp), intent(in) :: p,cuts(:),weights(:)
    type(distribution_component), intent(in) :: components(:)
    integer :: i
    real(dp) :: cumulative,low,high,mass,prob
    cumulative=0.0_dp; i=1
    do while(i<size(weights) .and. p>cumulative+weights(i)); cumulative=cumulative+weights(i); i=i+1; end do
    low=0.0_dp; if(i>1)low=component_cdf(cuts(i-1),components(i))
    high=1.0_dp; if(i<=size(cuts))high=component_cdf(cuts(i),components(i))
    mass=high-low; prob=low+mass*(p-cumulative)/weights(i)
    value=component_quantile(prob,components(i))
  end function composite_quantile

  real(dp) function composite_random(cuts,components,weights) result(value)
    real(dp), intent(in) :: cuts(:),weights(:)
    type(distribution_component), intent(in) :: components(:)
    value=composite_quantile(random_uniform(),cuts,components,weights)
  end function composite_random

  pure integer function bucket_index(x,cuts) result(index_value)
    real(dp), intent(in) :: x,cuts(:)
    integer :: i
    index_value=size(cuts)+1
    do i=1,size(cuts)
      if(x<=cuts(i)) then; index_value=i; return; end if
    end do
  end function bucket_index
end module qrmtools_distributions
