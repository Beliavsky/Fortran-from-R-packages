! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_jumps
  use highfrequency_kinds, only: dp, pi
  use highfrequency_stats, only: normal_cdf, normal_quantile, mu_abs_normal
  use highfrequency_types, only: jump_test_result, iv_inference_result
  use highfrequency_realized, only: realized_variance, bipower_variation
  use highfrequency_realized, only: tripower_quarticity, realized_quarticity
  use highfrequency_data, only: make_returns
  implicit none
  private
  public :: bns_jump_test, aj_jump_test, iv_inference, abd_jump_test

contains

  function bns_jump_test(returns, ratio_form, log_transform, confidence) result(result)
    real(dp), intent(in) :: returns(:)
    logical, intent(in), optional :: ratio_form, log_transform
    real(dp), intent(in), optional :: confidence
    type(jump_test_result) :: result
    real(dp) :: rv, bv, tq, product, theta, alpha, z
    logical :: use_ratio, use_log
    integer :: n
    n=size(returns)
    if(n<4)return
    use_ratio=.false.;if(present(ratio_form))use_ratio=ratio_form
    use_log=.false.;if(present(log_transform))use_log=log_transform
    alpha=0.975_dp;if(present(confidence))alpha=confidence
    rv=realized_variance(returns)
    bv=bipower_variation(returns)
    tq=tripower_quarticity(returns)
    theta=pi*pi/4.0_dp+pi-5.0_dp
    if(rv<=0.0_dp .or. bv<=0.0_dp .or. tq<=0.0_dp)return
    if(use_ratio)then
      product=max(1.0_dp,tq/(bv*bv))
      z=sqrt(real(n,dp))*(1.0_dp-bv/rv)/sqrt(theta*product)
    else
      product=tq
      if(use_log)then
        z=sqrt(real(n,dp))*(log(rv)-log(bv))/sqrt(theta*product/(bv*bv))
      else
        z=sqrt(real(n,dp))*(rv-bv)/sqrt(theta*product)
      end if
    end if
    result%statistic=z
    result%critical_lower=normal_quantile(1.0_dp-alpha)
    result%critical_upper=normal_quantile(alpha)
    result%p_value=2.0_dp*normal_cdf(-abs(z))
    result%reject=abs(z)>result%critical_upper
  end function bns_jump_test

  function abd_jump_test(returns, confidence) result(result)
    real(dp), intent(in) :: returns(:)
    real(dp), intent(in), optional :: confidence
    type(jump_test_result) :: result
    real(dp) :: rv,bv,tq,alpha,z,denom
    integer :: n
    n=size(returns)
    if(n<4)return
    alpha=0.975_dp;if(present(confidence))alpha=confidence
    rv=realized_variance(returns)
    bv=bipower_variation(returns)
    tq=tripower_quarticity(returns)
    denom=sqrt(max(tiny(1.0_dp),(pi*pi/4.0_dp+pi-5.0_dp)*tq/real(n,dp)))
    z=(rv-bv)/denom
    result%statistic=z
    result%critical_lower=normal_quantile(1.0_dp-alpha)
    result%critical_upper=normal_quantile(alpha)
    result%p_value=2.0_dp*normal_cdf(-abs(z))
    result%reject=z>result%critical_upper
  end function abd_jump_test

  function aj_jump_test(prices, p, k, alpha_multiplier, confidence) result(result)
    real(dp), intent(in) :: prices(:)
    real(dp), intent(in), optional :: p, alpha_multiplier, confidence
    integer, intent(in), optional :: k
    type(jump_test_result) :: result
    real(dp), allocatable :: r(:), coarse(:), selected(:)
    real(dp) :: power, alpha_mult, conf, rv, cutoff, pv1, pv2, s, variance
    integer :: kk,n,i,m,ns
    n=size(prices)-1
    if(n<20 .or. any(prices<=0.0_dp))return
    power=4.0_dp;if(present(p))power=p
    kk=2;if(present(k))kk=k
    alpha_mult=4.0_dp;if(present(alpha_multiplier))alpha_mult=alpha_multiplier
    conf=0.975_dp;if(present(confidence))conf=confidence
    if(kk<2)return
    r=abs(make_returns(prices))
    rv=realized_variance(r)
    cutoff=alpha_mult*sqrt(max(rv,tiny(1.0_dp)))*real(n,dp)**(-0.47_dp)
    m=n/kk
    allocate(coarse(m))
    do i=1,m
      coarse(i)=abs(log(prices(1+i*kk)/prices(1+(i-1)*kk)))
    end do
    pv1=sum(r**power)
    pv2=sum(coarse**power)
    if(pv1<=0.0_dp)return
    s=pv2/pv1
    ns=count(r<cutoff)
    if(ns<5)return
    allocate(selected(ns))
    selected=pack(r,r<cutoff)
    variance=aj_asymptotic_variance(selected,power,kk,n)
    if(variance<=0.0_dp)return
    result%statistic=(s-real(kk,dp)**(0.5_dp*power-1.0_dp))/sqrt(variance)
    result%critical_lower=normal_quantile(1.0_dp-conf)
    result%critical_upper=normal_quantile(conf)
    result%p_value=2.0_dp*normal_cdf(-abs(result%statistic))
    result%reject=abs(result%statistic)>result%critical_upper
  end function aj_jump_test

  pure real(dp) function aj_asymptotic_variance(r,p,k,n) result(value)
    real(dp),intent(in)::r(:),p
    integer,intent(in)::k,n
    real(dp)::m2p,mp,scale
    mp=mu_abs_normal(p)
    m2p=mu_abs_normal(2.0_dp*p)
    if(mp<=0.0_dp .or. size(r)<2)then
      value=0.0_dp
      return
    end if
    scale=sum(r**(2.0_dp*p))/max(tiny(1.0_dp),sum(r**p)**2)
    value=real(n,dp)*scale*(m2p/(mp*mp)-1.0_dp)* &
      real(k,dp)**(p-2.0_dp)*(real(k,dp)-1.0_dp)
  end function aj_asymptotic_variance

  function iv_inference(returns, use_bipower, confidence) result(result)
    real(dp),intent(in)::returns(:)
    logical,intent(in),optional::use_bipower
    real(dp),intent(in),optional::confidence
    type(iv_inference_result)::result
    real(dp)::estimate,iq,theta,z,conf
    logical::bp
    integer::n
    n=size(returns)
    if(n<4)return
    bp=.false.;if(present(use_bipower))bp=use_bipower
    conf=0.95_dp;if(present(confidence))conf=confidence
    if(bp)then
      estimate=bipower_variation(returns)
      theta=pi*pi/4.0_dp+pi-5.0_dp
      iq=tripower_quarticity(returns)
    else
      estimate=realized_variance(returns)
      theta=2.0_dp
      iq=realized_quarticity(returns)
    end if
    result%estimate=estimate
    result%standard_error=sqrt(max(0.0_dp,theta*iq/real(n,dp)))
    z=normal_quantile(0.5_dp+0.5_dp*conf)
    result%lower=estimate-z*result%standard_error
    result%upper=estimate+z*result%standard_error
  end function iv_inference

end module highfrequency_jumps
