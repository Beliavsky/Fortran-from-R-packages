! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_forecast_module
  use ghyp_kinds, only : dp, i8
  use tsgarch_types
  use tsgarch_simulation_module, only : simulate_conditional
  implicit none
  private
  public :: forecast_garch, type7_quantile
contains
  function forecast_garch(history,fit,horizon,probabilities,paths,seed,vreg_history,vreg_future) result(out)
    real(dp),intent(in)::history(:)
    type(garch_fit),intent(in)::fit
    integer,intent(in)::horizon
    real(dp),intent(in),optional::probabilities(:)
    integer,intent(in),optional::paths
    integer(i8),intent(in),optional::seed
    real(dp),intent(in),optional::vreg_history(:,:),vreg_future(:,:)
    type(garch_forecast)::out
    type(garch_simulation)::sim
    real(dp),allocatable::probs(:)
    integer::h,j,np
    if(fit%status/=tsg_success.and.fit%filtered%status/=tsg_success)then
    out%message='fit is not valid'
    return
    end if
    if(present(probabilities))then
      allocate(probs(size(probabilities)))
      probs=probabilities
    else
      allocate(probs(5))
      probs=[0.01_dp,0.05_dp,0.50_dp,0.95_dp,0.99_dp]
    end if
    if(any(probs<=0.0_dp).or.any(probs>=1.0_dp))then
    out%message='probabilities must lie in (0,1)'
    return
    end if
    np=1000
    if(present(paths))np=paths
    if(present(vreg_history))then
      sim=simulate_conditional(history,fit%spec,fit%parameters,horizon,np,seed,vreg_history,vreg_future)
    else
      sim=simulate_conditional(history,fit%spec,fit%parameters,horizon,np,seed)
    end if
    if(sim%status/=tsg_success)then
    out%message=sim%message
    return
    end if
    allocate(out%mean(horizon),out%sigma(horizon),out%variance(horizon))
    allocate(out%probabilities(size(probs)),out%quantiles(horizon,size(probs)))
    out%probabilities=probs
    do h=1,horizon
      out%mean(h)=sum(sim%series(h,:))/real(size(sim%series,2),dp)
      out%variance(h)=sum((sim%series(h,:)-out%mean(h))**2)/real(max(1,size(sim%series,2)-1),dp)
      out%sigma(h)=sqrt(max(out%variance(h),0.0_dp))
      do j=1,size(probs)
      out%quantiles(h,j)=type7_quantile(sim%series(h,:),probs(j))
      end do
    end do
    out%status=tsg_success
    out%message='ok'
  end function forecast_garch

  real(dp) function type7_quantile(x,p) result(value)
    real(dp),intent(in)::x(:),p
    real(dp),allocatable::s(:)
    real(dp)::h,g
    integer::j,n
    n=size(x)
    if(n==0)then
    value=0.0_dp
    return
    end if
    allocate(s(n))
    s=x
    call sort_real(s)
    if(p<=0.0_dp)then
    value=s(1)
    return
    else if(p>=1.0_dp)then
    value=s(n)
    return
    end if
    h=1.0_dp+real(n-1,dp)*p
    j=int(floor(h))
    g=h-real(j,dp)
    if(j>=n)then
    value=s(n)
    else
    value=(1.0_dp-g)*s(j)+g*s(j+1)
    end if
  end function type7_quantile

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::key
    do i=2,size(x)
    key=x(i)
    j=i-1
    do while(j>=1)
    if(x(j)<=key)exit
    x(j+1)=x(j)
    j=j-1
    end do
    x(j+1)=key
    end do
  end subroutine sort_real
end module tsgarch_forecast_module
