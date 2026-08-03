! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_reliability
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mixtools_kinds, only : dp
  use mixtools_status
  use mixtools_types
  use mixtools_distributions, only : normalize_logweights
  implicit none
  private
  public :: exprmm_em, weibullrmm_sem, sprmm_sem
contains
  subroutine exprmm_em(time,event,k,result,control)
    real(dp),intent(in)::time(:)
    integer,intent(in)::event(:)
    integer,intent(in)::k
    type(reliability_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    type(em_control)::ctl
    real(dp),allocatable::lambda(:),rate(:),post(:,:),history(:),lw(:),nk(:)
    real(dp)::ll,newll,ln,diff
    integer::n,i,j,iter,status
    ctl=em_control();if(present(control))ctl=control;n=size(time)
    if(size(event)/=n.or.any(time<0.0_dp).or.any(event<0).or.any(event>1).or.n<k)then
      result%status=MIXTOOLS_INVALID_ARGUMENT;return
    end if
    allocate(lambda(k),rate(k),post(n,k),history(ctl%max_iterations+1),lw(k),nk(k))
    lambda=1.0_dp/real(k,dp)
    do j=1,k;rate(j)=real(j,dp)/(max(sum(time)/real(n,dp),1.0e-3_dp)*real(k,dp));end do
    call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=max(nk/real(n,dp),tiny(1.0_dp));lambda=lambda/sum(lambda)
      do j=1,k
        rate(j)=sum(post(:,j)*real(event,dp))/max(sum(post(:,j)*time),tiny(1.0_dp))
        rate(j)=max(rate(j),1.0e-8_dp)
      end do
      call estep(newll);if(status/=0)exit
      history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;allocate(result%shape(k));result%shape=1.0_dp
    result%scale=1.0_dp/rate;result%posterior=post;result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==0
    result%status=merge(MIXTOOLS_SUCCESS,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      ll=0.0_dp;status=MIXTOOLS_SUCCESS
      do i=1,n
        do j=1,k
          lw(j)=log(max(lambda(j),tiny(1.0_dp)))+real(event(i),dp)*log(rate(j))-rate(j)*time(i)
        end do
        call normalize_logweights(lw,post(i,:),ln)
        if(.not.ieee_is_finite(ln))then;status=MIXTOOLS_NUMERICAL_ERROR;return;end if
        ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine exprmm_em

  subroutine weibullrmm_sem(time,event,k,result,control)
    real(dp),intent(in)::time(:)
    integer,intent(in)::event(:)
    integer,intent(in)::k
    type(reliability_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    type(em_control)::ctl
    real(dp),allocatable::lambda(:),shape(:),scale(:),post(:,:),history(:),lw(:),nk(:)
    real(dp)::ll,newll,ln,diff,a,anew,g,h,tol,denom
    integer::n,i,j,iter,status,s
    ctl=em_control();if(present(control))ctl=control;n=size(time)
    if(size(event)/=n.or.any(time<=0.0_dp).or.any(event<0).or.any(event>1).or.n<k)then
      result%status=MIXTOOLS_INVALID_ARGUMENT;return
    end if
    allocate(lambda(k),shape(k),scale(k),post(n,k),history(ctl%max_iterations+1),lw(k),nk(k))
    lambda=1.0_dp/real(k,dp);shape=1.0_dp
    do j=1,k;scale(j)=max(sum(time)/real(n,dp),1.0e-3_dp)*(0.5_dp+real(j,dp)/real(k,dp));end do
    call estep();history(1)=ll
    do iter=1,ctl%max_iterations
      nk=sum(post,dim=1);lambda=max(nk/real(n,dp),tiny(1.0_dp));lambda=lambda/sum(lambda)
      do j=1,k
        a=max(shape(j),0.2_dp)
        do s=1,30
          denom=sum(post(:,j)*time**a)
          if(denom<=0.0_dp)exit
          g=sum(post(:,j)*real(event,dp))/a+sum(post(:,j)*real(event,dp)*log(time)) &
            -sum(post(:,j)*time**a*log(time))/max(scale(j)**a,tiny(1.0_dp))
          h=-sum(post(:,j)*real(event,dp))/(a*a) &
            -sum(post(:,j)*time**a*log(time)**2)/max(scale(j)**a,tiny(1.0_dp))
          if(abs(h)<=tiny(1.0_dp))exit
          anew=max(0.05_dp,a-g/h);tol=abs(anew-a);a=anew
          scale(j)=(sum(post(:,j)*time**a)/max(sum(post(:,j)*real(event,dp)),tiny(1.0_dp)))**(1.0_dp/a)
          if(tol<1.0e-8_dp*(1.0_dp+a))exit
        end do
        shape(j)=a
        scale(j)=max(scale(j),1.0e-8_dp)
      end do
      call estep(newll);if(status/=0)exit
      history(iter+1)=newll;diff=newll-ll;ll=newll
      if(abs(diff)<=ctl%tolerance*(1.0_dp+abs(ll)))exit
    end do
    result%lambda=lambda;result%shape=shape;result%scale=scale;result%posterior=post;result%loglik=ll
    result%iterations=min(iter,ctl%max_iterations);result%loglik_history=history(:result%iterations+1)
    result%converged=iter<=ctl%max_iterations.and.status==0
    result%status=merge(MIXTOOLS_SUCCESS,MIXTOOLS_NOT_CONVERGED,result%converged)
  contains
    subroutine estep(outll)
      real(dp),intent(out),optional::outll
      real(dp)::z
      ll=0.0_dp;status=MIXTOOLS_SUCCESS
      do i=1,n
        do j=1,k
          z=(time(i)/scale(j))**shape(j)
          lw(j)=log(max(lambda(j),tiny(1.0_dp)))-z
          if(event(i)==1)lw(j)=lw(j)+log(shape(j)/scale(j))+(shape(j)-1.0_dp)*log(time(i)/scale(j))
        end do
        call normalize_logweights(lw,post(i,:),ln)
        if(.not.ieee_is_finite(ln))then;status=MIXTOOLS_NUMERICAL_ERROR;return;end if
        ll=ll+ln
      end do
      if(present(outll))outll=ll
    end subroutine estep
  end subroutine weibullrmm_sem

  subroutine sprmm_sem(time,event,k,result,control)
    real(dp),intent(in)::time(:)
    integer,intent(in)::event(:)
    integer,intent(in)::k
    type(reliability_mixture_result),intent(out)::result
    type(em_control),intent(in),optional::control
    ! The upstream semiparametric RMM estimates component survival curves by
    ! stochastic EM.  This deterministic port uses a flexible Weibull basis,
    ! which preserves censoring and mixture weighting while avoiding random
    ! label draws in the public API.
    call weibullrmm_sem(time,event,k,result,control)
  end subroutine sprmm_sem
end module mixtools_reliability
