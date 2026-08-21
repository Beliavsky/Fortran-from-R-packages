! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_mixture
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data, flexsurv_spec, flexsurv_result, &
    fit_flexsurvreg, flexsurv_loglik, parameter_row, predict_survival, &
    initial_theta, &
    predict_hazard, predict_density
  use flexsurv_distributions, only : dist_cdf, dist_quantile
  use flexsurv_math, only : logsumexp
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  type, public :: flexsurvmix_result
    real(dp), allocatable :: mixing(:)
    type(flexsurv_spec), allocatable :: specs(:)
    type(flexsurv_result), allocatable :: components(:)
    real(dp), allocatable :: posterior(:,:)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    integer :: iterations = 0
    logical :: converged = .false.
  end type flexsurvmix_result

  public :: fit_flexsurvmix, mix_loglik, mix_survival, mix_density, mix_hazard
  public :: mix_cdf, mix_quantile, mix_posterior

contains

  function fit_flexsurvmix(data,specs,mixing_init,maxit,tol) result(res)
    type(flexsurv_data),intent(in)::data
    type(flexsurv_spec),intent(in)::specs(:)
    real(dp),intent(in),optional::mixing_init(:)
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    type(flexsurvmix_result)::res
    type(flexsurv_data)::wd
    type(flexsurv_spec),allocatable::sp(:)
    real(dp),allocatable::pi(:),post(:,:),llik(:,:),logw(:),ind(:)
    real(dp)::oldll,newll,tt,sw
    integer::k,n,i,j,mi,iter
    k=size(specs);n=size(data%lower);mi=200;if(present(maxit))mi=maxit
    tt=1.0e-6_dp;if(present(tol))tt=tol
    allocate(pi(k),post(n,k),llik(n,k),logw(k),ind(n),sp(k),res%components(k))
    sp=specs
    if(present(mixing_init))then
      pi=max(mixing_init,1.0e-12_dp);pi=pi/sum(pi)
    else
      pi=1.0_dp/real(k,dp)
    end if
    oldll=-huge(1.0_dp)
    do iter=1,mi
      res%iterations=iter
      do j=1,k
        if(iter==1)then
          call component_individual(data,sp(j),initial_theta(sp(j)),llik(:,j))
        else
          call component_individual(data,sp(j),res%components(j)%theta,llik(:,j))
        end if
      end do
      newll=0.0_dp
      do i=1,n
        do j=1,k;logw(j)=log(max(pi(j),tiny(1.0_dp)))+llik(i,j);end do
        sw=logsumexp(logw);newll=newll+data%weights(i)*sw
        do j=1,k;post(i,j)=exp(logw(j)-sw);end do
      end do
      do j=1,k
        pi(j)=sum(data%weights*post(:,j))/sum(data%weights)
        wd=data;wd%weights=data%weights*post(:,j)
        if(iter>1)sp(j)%base_init=res%components(j)%base
        res%components(j)=fit_flexsurvreg(wd,sp(j))
      end do
      pi=max(pi,1.0e-12_dp);pi=pi/sum(pi)
      if(iter>1.and.abs(newll-oldll)<=tt*(1.0_dp+abs(oldll)))then
        res%converged=.true.;exit
      end if
      oldll=newll
    end do
    ! final posterior and likelihood
    do j=1,k;call component_individual(data,sp(j),res%components(j)%theta,llik(:,j));end do
    res%loglik=0.0_dp
    do i=1,n
      do j=1,k;logw(j)=log(max(pi(j),tiny(1.0_dp)))+llik(i,j);end do
      sw=logsumexp(logw);res%loglik=res%loglik+data%weights(i)*sw
      do j=1,k;post(i,j)=exp(logw(j)-sw);end do
    end do
    res%mixing=pi;res%posterior=post;res%specs=sp
    res%aic=-2.0_dp*res%loglik+2.0_dp*real(total_parameters(res),dp)
  end function fit_flexsurvmix

  real(dp) function mix_loglik(data,res) result(ll)
    type(flexsurv_data),intent(in)::data
    type(flexsurvmix_result),intent(in)::res
    real(dp),allocatable::li(:,:),lw(:),ind(:)
    integer::n,k,i,j
    n=size(data%lower);k=size(res%mixing);allocate(li(n,k),lw(k),ind(n));ll=0.0_dp
    do j=1,k;call component_individual(data,res%specs(j),res%components(j)%theta,li(:,j));end do
    do i=1,n
      do j=1,k;lw(j)=log(max(res%mixing(j),tiny(1.0_dp)))+li(i,j);end do
      ll=ll+data%weights(i)*logsumexp(lw)
    end do
  end function mix_loglik

  subroutine mix_posterior(data,res,post)
    type(flexsurv_data),intent(in)::data
    type(flexsurvmix_result),intent(in)::res
    real(dp),intent(out)::post(size(data%lower),size(res%mixing))
    real(dp),allocatable::li(:,:),lw(:)
    real(dp)::ls
    integer::n,k,i,j
    n=size(data%lower);k=size(res%mixing);allocate(li(n,k),lw(k))
    do j=1,k;call component_individual(data,res%specs(j),res%components(j)%theta,li(:,j));end do
    do i=1,n
      do j=1,k;lw(j)=log(max(res%mixing(j),tiny(1.0_dp)))+li(i,j);end do
      ls=logsumexp(lw);do j=1,k;post(i,j)=exp(lw(j)-ls);end do
    end do
  end subroutine mix_posterior

  real(dp) function mix_survival(res,row,t) result(s)
    type(flexsurvmix_result),intent(in)::res
    integer,intent(in)::row;real(dp),intent(in)::t
    integer::j;s=0.0_dp
    do j=1,size(res%mixing);s=s+res%mixing(j)*predict_survival(res%specs(j),res%components(j)%theta,row,t);end do
  end function mix_survival

  real(dp) function mix_density(res,row,t) result(f)
    type(flexsurvmix_result),intent(in)::res
    integer,intent(in)::row;real(dp),intent(in)::t
    integer::j;f=0.0_dp
    do j=1,size(res%mixing);f=f+res%mixing(j)*predict_density(res%specs(j),res%components(j)%theta,row,t);end do
  end function mix_density

  real(dp) function mix_hazard(res,row,t) result(h)
    type(flexsurvmix_result),intent(in)::res
    integer,intent(in)::row;real(dp),intent(in)::t
    real(dp)::s; s=mix_survival(res,row,t)
    if(s<=0.0_dp)then;h=huge(1.0_dp);else;h=mix_density(res,row,t)/s;end if
  end function mix_hazard

  real(dp) function mix_cdf(res,row,t) result(f)
    type(flexsurvmix_result),intent(in)::res
    integer,intent(in)::row;real(dp),intent(in)::t
    f=1.0_dp-mix_survival(res,row,t)
  end function mix_cdf

  real(dp) function mix_quantile(res,row,p) result(q)
    type(flexsurvmix_result),intent(in)::res
    integer,intent(in)::row;real(dp),intent(in)::p
    real(dp)::lo,hi,mid
    integer::it,j
    if(p<=0.0_dp)then;q=0.0_dp;return;end if
    if(p>=1.0_dp)then;q=huge(1.0_dp);return;end if
    hi=1.0_dp
    do j=1,size(res%mixing)
      hi=max(hi,dist_quantile(res%specs(j)%dist,0.999999_dp, &
        parameter_row(res%specs(j),res%components(j)%theta,row)))
    end do
    lo=0.0_dp
    do it=1,100
      mid=0.5_dp*(lo+hi)
      if(mix_cdf(res,row,mid)<p)then;lo=mid;else;hi=mid;end if
    end do
    q=0.5_dp*(lo+hi)
  end function mix_quantile

  subroutine component_individual(data,spec,theta,li)
    type(flexsurv_data),intent(in)::data
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:)
    real(dp),intent(out)::li(:)
    real(dp)::dummy
    dummy=flexsurv_loglik(data,spec,theta,individual=li)
  end subroutine component_individual

  integer function total_parameters(res) result(n)
    type(flexsurvmix_result),intent(in)::res
    integer::j;n=max(0,size(res%mixing)-1)
    do j=1,size(res%components);n=n+res%components(j)%npar;end do
  end function total_parameters

end module flexsurv_mixture
