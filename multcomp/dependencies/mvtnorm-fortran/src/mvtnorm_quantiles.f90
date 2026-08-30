! SPDX-License-Identifier: GPL-2.0-only
module mvtnorm_quantiles
  use mvtnorm_kinds, only : dp
  use mvtnorm_types, only : probability_control, probability_result, quantile_result
  use mvtnorm_probabilities, only : pmvnorm, pmvt
  implicit none
  private
  public :: qmvnorm, qmvt
contains

  function qmvnorm(p,mean,sigma,tail,control,lower_bound,upper_bound,ptol,maxiter) result(res)
    real(dp),intent(in)::p,mean(:),sigma(:,:)
    character(len=*),intent(in),optional::tail
    type(probability_control),intent(in),optional::control
    real(dp),intent(in),optional::lower_bound,upper_bound,ptol
    integer,intent(in),optional::maxiter
    type(quantile_result)::res
    character(len=16)::which
    which='lower'; if(present(tail)) which=adjustl(tail)
    res=quantile_core(p,mean,sigma,0.0_dp,which,control,lower_bound,upper_bound,ptol,maxiter)
  end function qmvnorm

  function qmvt(p,delta,sigma,df,tail,control,lower_bound,upper_bound,ptol,maxiter) result(res)
    real(dp),intent(in)::p,delta(:),sigma(:,:),df
    character(len=*),intent(in),optional::tail
    type(probability_control),intent(in),optional::control
    real(dp),intent(in),optional::lower_bound,upper_bound,ptol
    integer,intent(in),optional::maxiter
    type(quantile_result)::res
    character(len=16)::which
    which='lower'; if(present(tail)) which=adjustl(tail)
    res=quantile_core(p,delta,sigma,df,which,control,lower_bound,upper_bound,ptol,maxiter)
  end function qmvt

  function quantile_core(p,location,sigma,df,tail,control,lower_bound,upper_bound,ptol,maxiter) result(res)
    real(dp),intent(in)::p,location(:),sigma(:,:),df
    character(len=*),intent(in)::tail
    type(probability_control),intent(in),optional::control
    real(dp),intent(in),optional::lower_bound,upper_bound,ptol
    integer,intent(in),optional::maxiter
    type(quantile_result)::res
    real(dp)::lo,hi,mid,flo,fhi,fmid,tolerance,smax,center
    integer::iter,itmax,n
    n=size(location); tolerance=1.0e-4_dp; if(present(ptol)) tolerance=ptol
    itmax=100; if(present(maxiter)) itmax=maxiter
    if(p<=0.0_dp .or. p>=1.0_dp) then; res%message='p must be strictly between zero and one'; return; end if
    smax=sqrt(maxval([(sigma(iter,iter),iter=1,n)])); center=maxval(abs(location))
    lo=-center-10.0_dp*smax; hi=center+10.0_dp*smax
    if(index(tail,'both')>0) lo=0.0_dp
    if(present(lower_bound)) lo=lower_bound
    if(present(upper_bound)) hi=upper_bound
    flo=objective(lo,p,location,sigma,df,tail,control)
    fhi=objective(hi,p,location,sigma,df,tail,control)
    do while(flo*fhi>0.0_dp .and. abs(lo)<1.0e8_dp .and. abs(hi)<1.0e8_dp)
      if(abs(flo)<abs(fhi)) then; hi=hi+2.0_dp*(hi-lo); fhi=objective(hi,p,location,sigma,df,tail,control)
      else; lo=lo-2.0_dp*(hi-lo); if(index(tail,'both')>0) lo=0.0_dp; flo=objective(lo,p,location,sigma,df,tail,control); end if
    end do
    if(flo*fhi>0.0_dp) then; res%message='failed to bracket quantile'; return; end if
    mid=0.5_dp*(lo+hi)
    fmid=objective(mid,p,location,sigma,df,tail,control)
    do iter=1,itmax
      mid=0.5_dp*(lo+hi); fmid=objective(mid,p,location,sigma,df,tail,control)
      if(abs(fmid)<=tolerance .or. abs(hi-lo)<=tolerance*max(1.0_dp,abs(mid))) exit
      if(flo*fmid<=0.0_dp) then; hi=mid; fhi=fmid; else; lo=mid; flo=fmid; end if
    end do
    res%quantile=mid; res%probability=p+fmid; res%error=fmid; res%iterations=iter
    res%converged=abs(fmid)<=tolerance .or. abs(hi-lo)<=tolerance*max(1.0_dp,abs(mid))
    if(res%converged) then; res%message='normal completion'; else; res%message='maximum iterations reached'; end if
  end function quantile_core

  real(dp) function objective(q,p,location,sigma,df,tail,control) result(f)
    real(dp),intent(in)::q,p,location(:),sigma(:,:),df
    character(len=*),intent(in)::tail
    type(probability_control),intent(in),optional::control
    real(dp),allocatable::lo(:),up(:)
    type(probability_result)::pr
    integer::n
    n=size(location); allocate(lo(n),up(n))
    if(index(tail,'both')>0) then
      lo=-abs(q); up=abs(q)
    else if(index(tail,'upper')>0) then
      lo=q; up=1.0e100_dp
    else
      lo=-1.0e100_dp; up=q
    end if
    if(df>0.0_dp) then
      pr=pmvt(lo,up,location,sigma,df,control)
    else
      pr=pmvnorm(lo,up,location,sigma,control)
    end if
    if(index(tail,'upper')>0) then
      f=(1.0_dp-pr%value)-(1.0_dp-p)
    else
      f=pr%value-p
    end if
  end function objective
end module mvtnorm_quantiles
