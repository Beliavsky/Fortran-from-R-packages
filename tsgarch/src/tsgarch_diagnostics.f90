! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_diagnostics
  use ghyp_kinds, only : dp
  use tsd_distributions, only : pdist
  use tsd_math, only : invert_matrix
  use tsgarch_types
  implicit none
  private
  public :: probability_integral_transform, half_life, effective_sample_omega
  public :: covariance_opg, covariance_sandwich, confidence_intervals
contains
  function probability_integral_transform(fit) result(pit)
    type(garch_fit),intent(in)::fit
    real(dp),allocatable::pit(:)
    integer::i
    allocate(pit(fit%filtered%nobs))
    pit=0.0_dp
    do i=1,size(pit)
      pit(i)=pdist(fit%spec%distribution,fit%filtered%standardized_residuals(i),fit%parameters%dist)
    end do
  end function probability_integral_transform

  real(dp) function half_life(fit) result(value)
    type(garch_fit),intent(in)::fit
    real(dp)::p
    p=fit%filtered%persistence
    if(p<=0.0_dp)then
    value=0.0_dp
    else if(p>=1.0_dp)then
    value=huge(1.0_dp)
    else
    value=log(0.5_dp)/log(p)
    end if
  end function half_life

  real(dp) function effective_sample_omega(fit) result(value)
    type(garch_fit),intent(in)::fit
    value=fit%filtered%effective_omega
  end function effective_sample_omega

  subroutine covariance_opg(fit,covariance,status)
    type(garch_fit),intent(in)::fit
    real(dp),allocatable,intent(out)::covariance(:,:)
    integer,intent(out)::status
    real(dp),allocatable::opg(:,:)
    logical::ok
    if(.not.allocated(fit%scores))then
    allocate(covariance(0,0))
    status=tsg_invalid_argument
    return
    end if
    opg=matmul(transpose(fit%scores),fit%scores)
    call invert_matrix(opg,covariance,ok)
    status=merge(tsg_success,tsg_singular,ok)
  end subroutine covariance_opg

  subroutine covariance_sandwich(fit,covariance,status)
    type(garch_fit),intent(in)::fit
    real(dp),allocatable,intent(out)::covariance(:,:)
    integer,intent(out)::status
    real(dp),allocatable::hinv(:,:),meat(:,:)
    logical::ok
    if(.not.allocated(fit%scores).or..not.allocated(fit%hessian))then
    allocate(covariance(0,0))
    status=tsg_invalid_argument
    return
    end if
    call invert_matrix(fit%hessian,hinv,ok)
    if(.not.ok)then
    allocate(covariance(0,0))
    status=tsg_singular
    return
    end if
    meat=matmul(transpose(fit%scores),fit%scores)
    covariance=matmul(hinv,matmul(meat,hinv))
    status=tsg_success
  end subroutine covariance_sandwich

  subroutine confidence_intervals(fit,level,lower,upper,status)
    type(garch_fit),intent(in)::fit
    real(dp),intent(in),optional::level
    real(dp),allocatable,intent(out)::lower(:),upper(:)
    integer,intent(out)::status
    real(dp)::lev,z
    lev=0.95_dp
    if(present(level))lev=level
    if(lev<=0.0_dp.or.lev>=1.0_dp.or..not.allocated(fit%standard_errors))then
      allocate(lower(0),upper(0))
      status=tsg_invalid_argument
      return
    end if
    z=normal_quantile_local(0.5_dp+0.5_dp*lev)
    lower=fit%packed_parameters-z*fit%standard_errors
    upper=fit%packed_parameters+z*fit%standard_errors
    status=tsg_success
  end subroutine confidence_intervals

  real(dp) function normal_quantile_local(p) result(x)
    real(dp),intent(in)::p
    real(dp)::lo,hi,mid,cdf
    integer::i
    lo=-10.0_dp
    hi=10.0_dp
    do i=1,100
      mid=0.5_dp*(lo+hi)
      cdf=0.5_dp*erfc(-mid/sqrt(2.0_dp))
      if(cdf<p)then
      lo=mid
      else
      hi=mid
      end if
    end do
    x=0.5_dp*(lo+hi)
  end function normal_quantile_local
end module tsgarch_diagnostics
