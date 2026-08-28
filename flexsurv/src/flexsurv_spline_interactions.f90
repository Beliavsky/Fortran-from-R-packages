! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_spline_interactions
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data, parameter_regression, fs_status_event, &
    fs_status_left, fs_status_interval, bfgs_minimize
  use flexsurv_spline, only : survspline_model, survspline_logpdf, survspline_survival, &
    survspline_cdf, survspline_hazard
  use flexsurv_spline_fit, only : flexsurvspline_result
  use flexsurv_math, only : logdiffexp, near_positive_definite
  use numderiv, only : hessian
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  type, public :: spline_gamma_effect
    real(dp), allocatable :: beta(:)
  end type spline_gamma_effect

  type, public :: flexsurvspline_interaction_result
    type(survspline_model) :: model
    type(spline_gamma_effect), allocatable :: effect(:)
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: gradient(:)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    integer :: npar = 0
    integer :: iterations = 0
    integer :: status = 1
    logical :: converged = .false.
  end type flexsurvspline_interaction_result

  public :: fit_flexsurvspline_interactions
  public :: flexsurvspline_interaction_loglik
  public :: spline_interaction_gamma
  public :: spline_interaction_survival, spline_interaction_hazard
  public :: spline_interaction_to_basic

contains

  function fit_flexsurvspline_interactions(data,knots,scale,timescale,reg,gamma_init, &
      beta_init,maxit,tol,basis) result(res)
    type(flexsurv_data),intent(in)::data
    real(dp),intent(in)::knots(:)
    integer,intent(in)::scale,timescale
    type(parameter_regression),intent(in)::reg(:)
    real(dp),intent(in),optional::gamma_init(:),beta_init(:),tol
    integer,intent(in),optional::maxit,basis
    type(flexsurvspline_interaction_result)::res
    real(dp),allocatable::z0(:),zhat(:),g(:),hh(:,:),hpd(:,:),inv(:,:)
    integer::ng,np,mi,ist,j,pos,pj
    real(dp)::f,tt

    ng=size(knots)
    if(size(reg)/=ng)then
      res%status=20
      return
    end if
    np=ng
    do j=1,ng
      if(allocated(reg(j)%x))np=np+size(reg(j)%x,2)
    end do
    if(.not.valid_rows(reg,size(data%lower)))then
      res%status=21
      return
    end if
    allocate(z0(np))
    z0=0.0_dp
    if(present(gamma_init))then
      if(size(gamma_init)==ng)z0(1:ng)=gamma_init
    else
      z0(1)=0.0_dp
      if(ng>=2)z0(2)=1.0_dp
    end if
    if(present(beta_init))then
      if(size(beta_init)==np-ng)z0(ng+1:)=beta_init
    end if
    mi=600
    if(present(maxit))mi=maxit
    tt=1.0e-5_dp
    if(present(tol))tt=tol
    call bfgs_minimize(objective,z0,zhat,f,g,res%iterations,res%status,mi,tt)
    res%theta=zhat
    res%gradient=g
    res%npar=np
    res%loglik=-f
    res%aic=-2.0_dp*res%loglik+2.0_dp*real(np,dp)
    res%converged=res%status==0
    allocate(res%model%knots(ng),res%model%gamma(ng),res%effect(ng))
    res%model%knots=knots
    res%model%gamma=zhat(1:ng)
    res%model%scale=scale
    res%model%timescale=timescale
    if(present(basis))res%model%basis=basis
    pos=ng
    do j=1,ng
      pj=0
      if(allocated(reg(j)%x))pj=size(reg(j)%x,2)
      allocate(res%effect(j)%beta(pj))
      if(pj>0)then
        res%effect(j)%beta=zhat(pos+1:pos+pj)
        pos=pos+pj
      end if
    end do
    allocate(hh(np,np))
    call hessian(objective,zhat,hh)
    allocate(hpd(np,np))
    call near_positive_definite(hh,hpd,1.0e-9_dp)
    call invert_local(hpd,inv,ist)
    if(ist==0)res%covariance=inv
  contains
    real(dp) function objective(z) result(v)
      real(dp),intent(in)::z(:)
      v=-flexsurvspline_interaction_loglik(data,knots,scale,timescale,reg,z,basis)
      if(.not.ieee_is_finite(v))v=1.0e100_dp
    end function objective
  end function fit_flexsurvspline_interactions

  real(dp) function flexsurvspline_interaction_loglik(data,knots,scale,timescale,reg,z,basis) result(ll)
    type(flexsurv_data),intent(in)::data
    real(dp),intent(in)::knots(:),z(:)
    integer,intent(in)::scale,timescale
    type(parameter_regression),intent(in)::reg(:)
    integer,intent(in),optional::basis
    type(survspline_model)::m
    real(dp),allocatable::gam(:)
    real(dp)::li,fu,fl,pobs,hx
    integer::i,ng
    ng=size(knots)
    allocate(m%knots(ng),m%gamma(ng))
    m%knots=knots
    m%scale=scale
    m%timescale=timescale
    if(present(basis))m%basis=basis
    ll=0.0_dp
    do i=1,size(data%lower)
      gam=spline_interaction_gamma(reg,z,i,ng)
      m%gamma=gam
      pobs=survspline_cdf(m,data%rtrunc(i))-survspline_cdf(m,data%start(i))
      if(pobs<=0.0_dp.or..not.ieee_is_finite(pobs))then
        ll=-huge(1.0_dp)
        return
      end if
      select case(data%status(i))
      case(fs_status_event)
        li=survspline_logpdf(m,data%lower(i))
        if(data%bhazard(i)>0.0_dp)then
          hx=survspline_hazard(m,data%lower(i))
          if(hx<=0.0_dp)then
          ll=-huge(1.0_dp)
          return
          end if
          li=li+log(1.0_dp+data%bhazard(i)/hx)
        end if
      case(fs_status_left)
        li=log(max(survspline_cdf(m,data%upper(i)),tiny(1.0_dp)))
      case(fs_status_interval)
        fu=log(max(survspline_cdf(m,data%upper(i)),tiny(1.0_dp)))
        fl=log(max(survspline_cdf(m,data%lower(i)),tiny(1.0_dp)))
        li=logdiffexp(fu,fl)
      case default
        li=log(max(survspline_survival(m,data%lower(i)),tiny(1.0_dp)))
        if(data%bhazard(i)>0.0_dp)li=li+log(max(data%bcondsurv(i),tiny(1.0_dp)))
      end select
      li=li-log(pobs)
      if(.not.ieee_is_finite(li))then
      ll=-huge(1.0_dp)
      return
      end if
      ll=ll+data%weights(i)*li
    end do
  end function flexsurvspline_interaction_loglik

  function spline_interaction_gamma(reg,z,row,ng) result(gamma)
    type(parameter_regression),intent(in)::reg(:)
    real(dp),intent(in)::z(:)
    integer,intent(in)::row,ng
    real(dp)::gamma(ng)
    integer::j,pj,pos
    gamma=z(1:ng)
    pos=ng
    do j=1,ng
      pj=0
      if(allocated(reg(j)%x))pj=size(reg(j)%x,2)
      if(pj>0)then
        gamma(j)=gamma(j)+dot_product(reg(j)%x(row,:),z(pos+1:pos+pj))
        pos=pos+pj
      end if
    end do
  end function spline_interaction_gamma

  real(dp) function spline_interaction_survival(res,row,t,xreg) result(s)
    type(flexsurvspline_interaction_result),intent(in)::res
    integer,intent(in)::row
    real(dp),intent(in)::t
    type(parameter_regression),intent(in)::xreg(:)
    type(survspline_model)::m
    m=res%model
    m%gamma=gamma_from_result(res,row,xreg)
    s=survspline_survival(m,t)
  end function spline_interaction_survival

  real(dp) function spline_interaction_hazard(res,row,t,xreg) result(h)
    type(flexsurvspline_interaction_result),intent(in)::res
    integer,intent(in)::row
    real(dp),intent(in)::t
    type(parameter_regression),intent(in)::xreg(:)
    type(survspline_model)::m
    m=res%model
    m%gamma=gamma_from_result(res,row,xreg)
    h=survspline_hazard(m,t)
  end function spline_interaction_hazard

  function spline_interaction_to_basic(res,row,xreg) result(out)
    type(flexsurvspline_interaction_result),intent(in)::res
    integer,intent(in)::row
    type(parameter_regression),intent(in)::xreg(:)
    type(flexsurvspline_result)::out
    out%model=res%model
    out%model%gamma=gamma_from_result(res,row,xreg)
    allocate(out%beta(0))
    out%loglik=res%loglik
    out%aic=res%aic
    out%iterations=res%iterations
    out%status=res%status
    out%converged=res%converged
  end function spline_interaction_to_basic

  function gamma_from_result(res,row,xreg) result(gamma)
    type(flexsurvspline_interaction_result),intent(in)::res
    integer,intent(in)::row
    type(parameter_regression),intent(in)::xreg(:)
    real(dp)::gamma(size(res%model%gamma))
    integer::j,pj
    gamma=res%model%gamma
    do j=1,size(gamma)
      pj=size(res%effect(j)%beta)
      if(pj>0)gamma(j)=gamma(j)+dot_product(xreg(j)%x(row,:),res%effect(j)%beta)
    end do
  end function gamma_from_result

  logical function valid_rows(reg,n) result(ok)
    type(parameter_regression),intent(in)::reg(:)
    integer,intent(in)::n
    integer::j
    ok=.true.
    do j=1,size(reg)
      if(allocated(reg(j)%x))then
        if(size(reg(j)%x,1)/=n)then
        ok=.false.
        return
        end if
      end if
    end do
  end function valid_rows

  subroutine invert_local(a,ainv,status)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::ainv(:,:)
    integer,intent(out)::status
    real(dp),allocatable::aug(:,:),tmp(:)
    real(dp)::piv,fac
    integer::n,i,j,k,ip
    n=size(a,1)
    allocate(aug(n,2*n),tmp(2*n))
    aug=0.0_dp
    aug(:,1:n)=a
    do i=1,n
    aug(i,n+i)=1.0_dp
    end do
    status=0
    do i=1,n
      ip=i
      do k=i+1,n
        if(abs(aug(k,i))>abs(aug(ip,i)))ip=k
      end do
      if(abs(aug(ip,i))<1.0e-14_dp)then
        status=1
        allocate(ainv(n,n))
        ainv=0.0_dp
        return
      end if
      if(ip/=i)then
      tmp=aug(i,:)
      aug(i,:)=aug(ip,:)
      aug(ip,:)=tmp
      end if
      piv=aug(i,i)
      aug(i,:)=aug(i,:)/piv
      do j=1,n
        if(j/=i)then
        fac=aug(j,i)
        aug(j,:)=aug(j,:)-fac*aug(i,:)
        end if
      end do
    end do
    allocate(ainv(n,n))
    ainv=aug(:,n+1:)
  end subroutine invert_local

end module flexsurv_spline_interactions
