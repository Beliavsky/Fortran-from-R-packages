! SPDX-License-Identifier: GPL-2.0-only
module tsd_profile
  use ghyp_kinds, only : dp
  use ghyp_rng, only : rng_state
  use tsd_types, only : parameter_specification, distribution_fit, profile_summary, &
                        tsd_success, tsd_invalid_argument
  use tsd_distributions, only : rdist
  use tsd_fit, only : estimate_distribution
  use tsd_moments, only : dskewness, dkurtosis
  implicit none
  private
  public :: tsprofile
contains

  function tsprofile(specification,sizes,repetitions,rng,max_iterations) result(summary)
    type(parameter_specification),intent(in)::specification
    integer,intent(in)::sizes(:),repetitions
    type(rng_state),intent(inout)::rng
    integer,intent(in),optional::max_iterations
    type(profile_summary)::summary
    real(dp),allocatable::draw(:),errors(:,:,:),estimate(:),actual(:)
    type(distribution_fit)::fit
    integer::i,j,k,ns,maxit,nok
    ns=size(sizes);maxit=500;if(present(max_iterations))maxit=max_iterations
    if(ns==0.or.repetitions<1.or.any(sizes<5))then
      summary%status=tsd_invalid_argument;summary%message='positive repetitions and sample sizes >= 5 are required';return
    end if
    allocate(summary%sizes(ns));summary%sizes=sizes
    allocate(summary%rmse(7,ns),summary%mae(7,ns),summary%mape(7,ns));summary%rmse=0.0_dp;summary%mae=0.0_dp;summary%mape=0.0_dp
    allocate(errors(7,repetitions,ns));errors=0.0_dp
    actual=parameter_vector(specification)
    summary%actual=actual
    summary%distribution=specification%distribution
    nok=0
    do k=1,ns
      do j=1,repetitions
        draw=rdist(specification%distribution,sizes(k),rng,specification%parameters)
        fit=estimate_distribution(draw,specification,max_iterations=maxit,use_hessian=.false.)
        summary%attempted_fits=summary%attempted_fits+1
        if(fit%status==tsd_success)then
          estimate=fit_vector(fit)
          errors(:,j,k)=estimate-actual;nok=nok+1;summary%successful_fits=summary%successful_fits+1
        else
          errors(:,j,k)=0.0_dp
        end if
      end do
      if(nok>0)then
        do i=1,7
          summary%rmse(i,k)=sqrt(sum(errors(i,:,k)**2)/real(max(count(abs(errors(i,:,k))>0.0_dp),1),dp))
          summary%mae(i,k)=sum(abs(errors(i,:,k)))/real(max(count(abs(errors(i,:,k))>0.0_dp),1),dp)
          if(abs(actual(i))>1.0e-12_dp)summary%mape(i,k)=100.0_dp*summary%mae(i,k)/abs(actual(i))
        end do
      end if
      nok=0
    end do
    if(summary%successful_fits>0)then;summary%status=tsd_success;summary%message='simulation profile completed';else
      summary%status=tsd_invalid_argument;summary%message='all simulated fits failed';end if
  end function tsprofile

  function parameter_vector(spec) result(v)
    type(parameter_specification),intent(in)::spec
    real(dp),allocatable::v(:)
    allocate(v(7))
    v=[spec%parameters%mu,spec%parameters%sigma,spec%parameters%skew,spec%parameters%shape,spec%parameters%lambda, &
       dskewness(spec%distribution,spec%parameters%skew,spec%parameters%shape,spec%parameters%lambda), &
       dkurtosis(spec%distribution,spec%parameters%skew,spec%parameters%shape,spec%parameters%lambda)]
  end function parameter_vector

  function fit_vector(fit) result(v)
    type(distribution_fit),intent(in)::fit
    real(dp),allocatable::v(:)
    allocate(v(7))
    v=[fit%parameters%mu,fit%parameters%sigma,fit%parameters%skew,fit%parameters%shape,fit%parameters%lambda, &
       dskewness(fit%distribution,fit%parameters%skew,fit%parameters%shape,fit%parameters%lambda), &
       dkurtosis(fit%distribution,fit%parameters%skew,fit%parameters%shape,fit%parameters%lambda)]
  end function fit_vector

end module tsd_profile
