! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_custom
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data, fs_status_event, fs_status_left, &
    fs_status_interval, bfgs_minimize
  use flexsurv_math, only : logdiffexp, near_positive_definite
  use numderiv, only : hessian
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  abstract interface
    real(dp) function custom_scalar_fn(t,par) result(v)
      import dp
      real(dp),intent(in)::t,par(:)
    end function custom_scalar_fn
  end interface

  type, public :: custom_distribution
    procedure(custom_scalar_fn), pointer, nopass :: logpdf => null()
    procedure(custom_scalar_fn), pointer, nopass :: cdf => null()
    procedure(custom_scalar_fn), pointer, nopass :: hazard => null()
  end type custom_distribution

  type, public :: custom_fit_result
    real(dp),allocatable::parameters(:),theta(:),covariance(:,:),gradient(:)
    real(dp)::loglik=-huge(1.0_dp),aic=huge(1.0_dp)
    integer::iterations=0,status=1
    logical::converged=.false.
  end type custom_fit_result

  public :: fit_custom_survival, custom_loglik, custom_survival

contains

  function fit_custom_survival(data,dist,par_init,positive,maxit,tol) result(res)
    type(flexsurv_data),intent(in)::data
    type(custom_distribution),intent(in)::dist
    real(dp),intent(in)::par_init(:)
    logical,intent(in),optional::positive(:)
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    type(custom_fit_result)::res
    logical,allocatable::pos(:)
    real(dp),allocatable::z0(:),zhat(:),g(:),hh(:,:),hpd(:,:),inv(:,:)
    real(dp)::f,tt
    integer::i,mi,st,n
    n=size(par_init);allocate(pos(n));pos=.false.;if(present(positive))pos=positive
    allocate(z0(n));do i=1,n;if(pos(i))then;z0(i)=log(max(par_init(i),tiny(1.0_dp)));else;z0(i)=par_init(i);end if;end do
    mi=500;if(present(maxit))mi=maxit;tt=1.0e-8_dp;if(present(tol))tt=tol
    call bfgs_minimize(objective,z0,zhat,f,g,res%iterations,res%status,mi,tt)
    allocate(res%parameters(n));do i=1,n;if(pos(i))then;res%parameters(i)=exp(zhat(i));else;res%parameters(i)=zhat(i);end if;end do
    res%theta=zhat;res%gradient=g;res%loglik=-f;res%aic=-2.0_dp*res%loglik+2.0_dp*real(n,dp);res%converged=res%status==0
    allocate(hh(n,n),hpd(n,n));call hessian(objective,zhat,hh);call near_positive_definite(hh,hpd,1.0e-9_dp)
    call invert_local(hpd,inv,st);if(st==0)res%covariance=inv
  contains
    real(dp) function objective(z) result(v)
      real(dp),intent(in)::z(:)
      real(dp),allocatable::p(:)
      integer::j
      allocate(p(n));do j=1,n;if(pos(j))then;p(j)=exp(min(z(j),700.0_dp));else;p(j)=z(j);end if;end do
      v=-custom_loglik(data,dist,p)
      if(.not.ieee_is_finite(v))v=1.0e100_dp
    end function objective
  end function fit_custom_survival

  real(dp) function custom_loglik(data,dist,par) result(ll)
    type(flexsurv_data),intent(in)::data
    type(custom_distribution),intent(in)::dist
    real(dp),intent(in)::par(:)
    real(dp)::li,fu,fl,pobs,hx
    integer::i
    if(.not.associated(dist%logpdf).or..not.associated(dist%cdf))then;ll=-huge(1.0_dp);return;end if
    ll=0.0_dp
    do i=1,size(data%lower)
      pobs=dist%cdf(data%rtrunc(i),par)-dist%cdf(data%start(i),par)
      if(pobs<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
      select case(data%status(i))
      case(fs_status_event)
        li=dist%logpdf(data%lower(i),par)
        if(data%bhazard(i)>0.0_dp.and.associated(dist%hazard))then
          hx=dist%hazard(data%lower(i),par);li=li+log(1.0_dp+data%bhazard(i)/max(hx,tiny(1.0_dp)))
        end if
      case(fs_status_left)
        li=log(max(dist%cdf(data%upper(i),par),tiny(1.0_dp)))
      case(fs_status_interval)
        fu=log(max(dist%cdf(data%upper(i),par),tiny(1.0_dp)))
        fl=log(max(dist%cdf(data%lower(i),par),tiny(1.0_dp)));li=logdiffexp(fu,fl)
      case default
        li=log(max(1.0_dp-dist%cdf(data%lower(i),par),tiny(1.0_dp)))
        if(data%bhazard(i)>0.0_dp)li=li+log(max(data%bcondsurv(i),tiny(1.0_dp)))
      end select
      li=li-log(pobs);ll=ll+data%weights(i)*li
    end do
  end function custom_loglik

  real(dp) function custom_survival(dist,t,par) result(s)
    type(custom_distribution),intent(in)::dist
    real(dp),intent(in)::t,par(:)
    s=1.0_dp-dist%cdf(t,par)
  end function custom_survival

  subroutine invert_local(a,ainv,status)
    real(dp),intent(in)::a(:,:);real(dp),allocatable,intent(out)::ainv(:,:);integer,intent(out)::status
    real(dp),allocatable::aug(:,:),tmp(:);real(dp)::piv,fac;integer::n,i,j,k,ip
    n=size(a,1);allocate(aug(n,2*n),tmp(2*n));aug=0.0_dp;aug(:,1:n)=a
    do i=1,n;aug(i,n+i)=1.0_dp;end do;status=0
    do i=1,n
      ip=i;do k=i+1,n;if(abs(aug(k,i))>abs(aug(ip,i)))ip=k;end do
      if(abs(aug(ip,i))<1e-12_dp)then;status=1;allocate(ainv(n,n));ainv=0.0_dp;return;end if
      if(ip/=i)then;tmp=aug(i,:);aug(i,:)=aug(ip,:);aug(ip,:)=tmp;end if
      piv=aug(i,i);aug(i,:)=aug(i,:)/piv
      do j=1,n;if(j/=i)then;fac=aug(j,i);aug(j,:)=aug(j,:)-fac*aug(i,:);end if;end do
    end do
    allocate(ainv(n,n));ainv=aug(:,n+1:)
  end subroutine invert_local
end module flexsurv_custom
