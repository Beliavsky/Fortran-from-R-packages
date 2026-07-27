! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_risk
  use fextremes_kinds, only: dp, huge_penalty
  use fextremes_stats, only: quantile_type1, mean_value, normal_quantile, chi_square1_quantile, nan_value
  use fextremes_distributions, only: gev_quantile, gev_logpdf
  use fextremes_fit, only: gev_fit_result, gpd_fit_result, fit_gpd
  use fextremes_optimize, only: nelder_mead_bounded
  implicit none
  private
  public :: risk_result, threshold_stability_result, return_level_result, tail_profile_result
  public :: sample_var, sample_cvar, gpd_risk_measures, gpd_threshold_stability
  public :: gev_return_level, gev_return_level_profile, gpd_profile_risk, gpd_tail_curve

  type :: risk_result
    real(dp), allocatable :: probability(:), value_at_risk(:), expected_shortfall(:)
    real(dp), allocatable :: var_se(:), es_se(:)
  end type risk_result
  type :: threshold_stability_result
    real(dp), allocatable :: threshold(:), xi(:), beta(:), modified_scale(:), xi_se(:), beta_se(:), high_quantile(:)
    logical, allocatable :: converged(:)
  end type threshold_stability_result
  type :: return_level_result
    real(dp) :: estimate=0.0_dp, lower=0.0_dp, upper=0.0_dp, standard_error=0.0_dp
    real(dp) :: blocks=0.0_dp, confidence=0.95_dp
  end type return_level_result
  type :: tail_profile_result
    real(dp) :: var_estimate=0.0_dp, var_lower=0.0_dp, var_upper=0.0_dp
    real(dp) :: es_estimate=0.0_dp, es_lower=0.0_dp, es_upper=0.0_dp
    real(dp) :: probability=0.99_dp, confidence=0.95_dp
  end type tail_profile_result

  real(dp), allocatable, save :: profile_data(:), tail_excess(:)
  real(dp), save :: profile_pp=0.0_dp, profile_level=0.0_dp
  real(dp), save :: tail_a=0.0_dp, tail_threshold=0.0_dp, tail_target=0.0_dp
  integer, save :: tail_mode=1
contains
  real(dp) function sample_var(x,alpha,upper_tail) result(v)
    real(dp),intent(in)::x(:),alpha
    logical,intent(in),optional::upper_tail
    logical::up
    up=.false.; if(present(upper_tail)) up=upper_tail
    if(up) then; v=quantile_type1(x,1.0_dp-alpha); else; v=quantile_type1(x,alpha); end if
  end function sample_var

  real(dp) function sample_cvar(x,alpha,upper_tail) result(v)
    real(dp),intent(in)::x(:),alpha
    logical,intent(in),optional::upper_tail
    logical::up; real(dp)::q
    up=.false.; if(present(upper_tail)) up=upper_tail
    q=sample_var(x,alpha,up)
    if(up) then
      if(count(x>=q)>0) then; v=mean_value(pack(x,x>=q)); else; v=q; end if
    else
      if(count(x<=q)>0) then; v=mean_value(pack(x,x<=q)); else; v=q; end if
    end if
  end function sample_cvar

  subroutine gpd_risk_measures(fit,probabilities,result)
    type(gpd_fit_result),intent(in)::fit
    real(dp),intent(in)::probabilities(:)
    type(risk_result),intent(out)::result
    real(dp)::a,g,gx,q,es,dqdx,dqdb,desdx,desdb,numer,la
    real(dp)::grad(2),cov(2,2)
    integer::i
    allocate(result%probability(size(probabilities)),result%value_at_risk(size(probabilities)), &
      result%expected_shortfall(size(probabilities)), &
      result%var_se(size(probabilities)),result%es_se(size(probabilities)))
    result%probability=probabilities; cov=fit%covariance
    do i=1,size(probabilities)
      a=(1.0_dp-probabilities(i))/fit%exceedance_probability
      if(a<=0.0_dp) then
        result%value_at_risk(i)=nan_value(); result%expected_shortfall(i)=nan_value(); cycle
      end if
      la=log(a)
      if(abs(fit%xi)<1.0e-7_dp) then
        g=-la; gx=0.5_dp*la*la
      else
        g=(a**(-fit%xi)-1.0_dp)/fit%xi
        gx=(-fit%xi*a**(-fit%xi)*la-(a**(-fit%xi)-1.0_dp))/fit%xi**2
      end if
      q=fit%threshold+fit%beta*g; dqdx=fit%beta*gx; dqdb=g
      result%value_at_risk(i)=q
      grad=[dqdx,dqdb]; result%var_se(i)=sqrt(max(0.0_dp,dot_product(grad,matmul(cov,grad))))
      if(fit%xi<1.0_dp) then
        numer=q+fit%beta-fit%xi*fit%threshold; es=numer/(1.0_dp-fit%xi)
        result%expected_shortfall(i)=es
        desdb=(dqdb+1.0_dp)/(1.0_dp-fit%xi)
        desdx=((dqdx-fit%threshold)*(1.0_dp-fit%xi)+numer)/(1.0_dp-fit%xi)**2
        grad=[desdx,desdb]; result%es_se(i)=sqrt(max(0.0_dp,dot_product(grad,matmul(cov,grad))))
      else
        result%expected_shortfall(i)=huge(1.0_dp); result%es_se(i)=huge(1.0_dp)
      end if
    end do
  end subroutine gpd_risk_measures

  subroutine gpd_threshold_stability(x,thresholds,probability,result)
    real(dp),intent(in)::x(:),thresholds(:),probability
    type(threshold_stability_result),intent(out)::result
    type(gpd_fit_result)::fit
    type(risk_result)::risk
    integer::i
    allocate(result%threshold(size(thresholds)),result%xi(size(thresholds)),result%beta(size(thresholds)), &
      result%modified_scale(size(thresholds)),result%xi_se(size(thresholds)),result%beta_se(size(thresholds)), &
      result%high_quantile(size(thresholds)),result%converged(size(thresholds)))
    result%threshold=thresholds
    do i=1,size(thresholds)
      call fit_gpd(x,thresholds(i),fit,'mle','observed')
      result%xi(i)=fit%xi; result%beta(i)=fit%beta; result%modified_scale(i)=fit%beta-fit%xi*thresholds(i)
      result%xi_se(i)=fit%se(1); result%beta_se(i)=fit%se(2); result%converged(i)=fit%converged
      call gpd_risk_measures(fit,[probability],risk); result%high_quantile(i)=risk%value_at_risk(1)
    end do
  end subroutine gpd_threshold_stability


  subroutine gpd_tail_curve(fit, xvalues, survival)
    type(gpd_fit_result),intent(in)::fit
    real(dp),intent(in)::xvalues(:)
    real(dp),intent(out)::survival(size(xvalues))
    real(dp)::z
    integer::i
    do i=1,size(xvalues)
      if (xvalues(i) < fit%threshold) then
        survival(i)=1.0_dp
      else if (abs(fit%xi)<1.0e-10_dp) then
        survival(i)=fit%exceedance_probability*exp(-(xvalues(i)-fit%threshold)/fit%beta)
      else
        z=1.0_dp+fit%xi*(xvalues(i)-fit%threshold)/fit%beta
        if (z<=0.0_dp) then; survival(i)=0.0_dp
        else; survival(i)=fit%exceedance_probability*z**(-1.0_dp/fit%xi); end if
      end if
    end do
  end subroutine gpd_tail_curve

  subroutine gpd_profile_risk(fit,probability,confidence,result,n_grid)
    type(gpd_fit_result),intent(in)::fit
    real(dp),intent(in)::probability,confidence
    type(tail_profile_result),intent(out)::result
    integer,intent(in),optional::n_grid
    type(risk_result)::base
    integer::ng,i,it,ev
    real(dp)::start(1),lower(1),upper(1),best(1),fb,cutoff,peak,lo,hi,target
    real(dp),allocatable::grid(:),prof(:)
    logical::conv
    call gpd_risk_measures(fit,[probability],base)
    result%probability=probability; result%confidence=confidence
    result%var_estimate=base%value_at_risk(1); result%es_estimate=base%expected_shortfall(1)
    ng=81; if(present(n_grid)) ng=max(31,n_grid)
    tail_excess=fit%exceedances-fit%threshold
    tail_a=(1.0_dp-probability)/fit%exceedance_probability
    tail_threshold=fit%threshold
    start=[fit%xi]; lower=[-2.0_dp]; upper=[5.0_dp]
    peak=-fit%nll; cutoff=peak-0.5_dp*chi_square1_quantile(confidence)
    allocate(grid(ng),prof(ng))
    tail_mode=1
    lo=max(fit%threshold+epsilon(1.0_dp),result%var_estimate-5.0_dp*max(base%var_se(1),0.2_dp*abs(result%var_estimate)))
    hi=result%var_estimate+5.0_dp*max(base%var_se(1),0.2_dp*abs(result%var_estimate))
    do i=1,ng
      target=lo+(hi-lo)*real(i-1,dp)/real(ng-1,dp); grid(i)=target; tail_target=target
      call nelder_mead_bounded(tail_profile_obj,start,lower,upper,best,fb,conv,it,ev,max_iter=500)
      prof(i)=-fb
    end do
    result%var_lower=result%var_estimate; result%var_upper=result%var_estimate
    if(any(prof>=cutoff)) then
      result%var_lower=minval(pack(grid,prof>=cutoff)); result%var_upper=maxval(pack(grid,prof>=cutoff))
    end if
    tail_mode=2; upper=[0.99_dp]
    lo=max(result%var_estimate,result%es_estimate-5.0_dp*max(base%es_se(1),0.2_dp*abs(result%es_estimate)))
    hi=result%es_estimate+5.0_dp*max(base%es_se(1),0.2_dp*abs(result%es_estimate))
    do i=1,ng
      target=lo+(hi-lo)*real(i-1,dp)/real(ng-1,dp); grid(i)=target; tail_target=target
      call nelder_mead_bounded(tail_profile_obj,start,lower,upper,best,fb,conv,it,ev,max_iter=500)
      prof(i)=-fb
    end do
    result%es_lower=result%es_estimate; result%es_upper=result%es_estimate
    if(any(prof>=cutoff)) then
      result%es_lower=minval(pack(grid,prof>=cutoff)); result%es_upper=maxval(pack(grid,prof>=cutoff))
    end if
    if(allocated(tail_excess)) deallocate(tail_excess)
  end subroutine gpd_profile_risk

  real(dp) function tail_profile_obj(p) result(v)
    real(dp),intent(in)::p(:)
    real(dp)::xi,beta,g,z
    integer::i
    xi=p(1)
    if(abs(xi)<1.0e-8_dp) then; g=-log(tail_a)
    else; g=(tail_a**(-xi)-1.0_dp)/xi; end if
    if(tail_mode==1) then
      if(abs(g)<1.0e-12_dp) then; v=huge_penalty; return; end if
      beta=(tail_target-tail_threshold)/g
    else
      beta=(1.0_dp-xi)*(tail_target-tail_threshold)/(g+1.0_dp)
    end if
    if(beta<=0.0_dp) then; v=huge_penalty; return; end if
    v=0.0_dp
    do i=1,size(tail_excess)
      if(abs(xi)<1.0e-10_dp) then
        v=v+log(beta)+tail_excess(i)/beta
      else
        z=1.0_dp+xi*tail_excess(i)/beta
        if(z<=0.0_dp) then; v=huge_penalty; return; end if
        v=v+log(beta)+(1.0_dp/xi+1.0_dp)*log(z)
      end if
    end do
  end function tail_profile_obj

  subroutine gev_return_level(fit,blocks,confidence,result)
    type(gev_fit_result),intent(in)::fit
    real(dp),intent(in)::blocks,confidence
    type(return_level_result),intent(out)::result
    real(dp)::p,z,grad(3),h,params(3)
    integer::j
    p=1.0_dp-1.0_dp/blocks
    result%estimate=gev_quantile(p,fit%xi,fit%mu,fit%beta)
    params=[fit%xi,fit%mu,fit%beta]
    result%blocks=blocks; result%confidence=confidence; grad=0.0_dp
    do j=1,3
      h=1.0e-5_dp*max(1.0_dp,abs(params(j)))
      select case(j)
      case(1); grad(j)=(gev_quantile(p,fit%xi+h,fit%mu,fit%beta)-gev_quantile(p,fit%xi-h,fit%mu,fit%beta))/(2.0_dp*h)
      case(2); grad(j)=1.0_dp
      case(3); grad(j)=(result%estimate-fit%mu)/fit%beta
      end select
    end do
    result%standard_error=sqrt(max(0.0_dp,dot_product(grad,matmul(fit%covariance,grad))))
    z=normal_quantile(0.5_dp*(1.0_dp+confidence))
    result%lower=result%estimate-z*result%standard_error
    result%upper=result%estimate+z*result%standard_error
  end subroutine gev_return_level

  subroutine gev_return_level_profile(data,fit,blocks,confidence,result,n_grid)
    real(dp),intent(in)::data(:)
    type(gev_fit_result),intent(in)::fit
    real(dp),intent(in)::blocks,confidence
    type(return_level_result),intent(out)::result
    integer,intent(in),optional::n_grid
    integer::ng,i,it,ev; real(dp)::lo_rl,hi_rl,rl,cutoff,bestll,start(2),lower(2),upper(2),best(2),fb
    logical::conv; real(dp),allocatable::levels(:),prof(:)
    call gev_return_level(fit,blocks,confidence,result); ng=121; if(present(n_grid)) ng=max(31,n_grid)
    lo_rl=min(minval(data),result%estimate-5.0_dp*max(result%standard_error,0.1_dp*abs(result%estimate)+1.0_dp))
    hi_rl=max(maxval(data),result%estimate+5.0_dp*max(result%standard_error,0.1_dp*abs(result%estimate)+1.0_dp))
    allocate(levels(ng),prof(ng)); profile_data=data; profile_pp=1.0_dp/blocks
    bestll=-fit%nll; cutoff=bestll-0.5_dp*chi_square1_quantile(confidence)
    start=[fit%xi,log(fit%beta)]
    lower=[-2.0_dp,log(1.0e-8_dp)]
    upper=[2.0_dp, &
      log(max(1.0_dp,100.0_dp*(maxval(data)-minval(data)+1.0_dp)))]
    do i=1,ng
      rl=lo_rl+(hi_rl-lo_rl)*real(i-1,dp)/real(ng-1,dp); levels(i)=rl; profile_level=rl
      call nelder_mead_bounded(profile_obj,start,lower,upper,best,fb,conv,it,ev,max_iter=800)
      prof(i)=-fb
    end do
    if (any(prof>=cutoff)) then
      result%lower=minval(pack(levels,prof>=cutoff))
      result%upper=maxval(pack(levels,prof>=cutoff))
    end if
    if(allocated(profile_data)) deallocate(profile_data)
  end subroutine gev_return_level_profile

  real(dp) function profile_obj(p) result(v)
    real(dp),intent(in)::p(:)
    real(dp)::xi,beta,mu,prob,t
    integer::i
    xi=p(1); beta=exp(p(2)); prob=1.0_dp-profile_pp
    if(abs(xi)<1.0e-8_dp) then; mu=profile_level+beta*log(-log(prob))
    else; mu=profile_level+beta*(1.0_dp-(-log(prob))**(-xi))/xi; end if
    v=0.0_dp
    do i=1,size(profile_data)
      t=gev_logpdf(profile_data(i),xi,mu,beta)
      if(t<=-0.5_dp*huge(1.0_dp)) then; v=huge_penalty; return; end if
      v=v-t
    end do
  end function profile_obj
end module fextremes_risk
