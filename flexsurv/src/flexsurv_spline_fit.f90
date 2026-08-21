! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_spline_fit
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data, fs_status_event, fs_status_left, &
    fs_status_interval, bfgs_minimize
  use flexsurv_spline, only : survspline_model, survspline_logpdf, &
    survspline_survival, survspline_cdf
  use flexsurv_math, only : logdiffexp, near_positive_definite
  use numderiv, only : hessian
  use quadprog, only : solve_qp, qp_result
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  type, public :: flexsurvspline_result
    type(survspline_model) :: model
    real(dp), allocatable :: beta(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: gradient(:)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: status = 1
    logical :: converged = .false.
  end type flexsurvspline_result

  public :: fit_flexsurvspline, flexsurvspline_loglik
  public :: spline_predict_survival, spline_predict_hazard
  public :: spline_qp_initial

contains

  function fit_flexsurvspline(data, knots, scale, timescale, gamma_init, x, &
      beta_init, maxit, tol, basis) result(res)
    type(flexsurv_data), intent(in) :: data
    real(dp), intent(in) :: knots(:)
    integer, intent(in) :: scale, timescale
    real(dp), intent(in), optional :: gamma_init(:), x(:,:), beta_init(:)
    integer, intent(in), optional :: maxit,basis
    real(dp), intent(in), optional :: tol
    type(flexsurvspline_result) :: res
    real(dp), allocatable :: z0(:), zhat(:), g(:), hh(:,:), hpd(:,:), inv(:,:)
    integer :: ng, p, mi, ist
    real(dp) :: f, tt
    ng=size(knots); p=0
    if(present(x))p=size(x,2)
    allocate(z0(ng+p)); z0=0.0_dp
    if(present(gamma_init))then
      if(size(gamma_init)==ng)z0(1:ng)=gamma_init
    else
      call heuristic_gamma(data,knots,z0(1:ng))
    end if
    if(p>0.and.present(beta_init))then
      if(size(beta_init)==p)z0(ng+1:)=beta_init
    end if
    mi=500;if(present(maxit))mi=maxit
    tt=1.0e-7_dp;if(present(tol))tt=tol
    call bfgs_minimize(objective,z0,zhat,f,g,res%iterations,res%status,mi,tt)
    allocate(res%model%knots(ng),res%model%gamma(ng))
    res%model%knots=knots;res%model%gamma=zhat(1:ng)
    res%model%scale=scale;res%model%timescale=timescale
    if(present(basis))res%model%basis=basis
    allocate(res%beta(p));if(p>0)res%beta=zhat(ng+1:)
    res%gradient=g;res%loglik=-f;res%aic=-2.0_dp*res%loglik+2.0_dp*size(zhat)
    res%converged=res%status==0
    allocate(hh(size(zhat),size(zhat)))
    call hessian(objective,zhat,hh)
    allocate(hpd(size(hh,1),size(hh,2))); call near_positive_definite(hh,hpd,1.0e-9_dp)
    call invert_local(hpd,inv,ist)
    if(ist==0)res%covariance=inv
  contains
    real(dp) function objective(z) result(v)
      real(dp),intent(in)::z(:)
      v=-flexsurvspline_loglik(data,knots,scale,timescale,z,x,basis)
      if(.not.ieee_is_finite(v))v=1.0e100_dp
    end function objective
  end function fit_flexsurvspline

  real(dp) function flexsurvspline_loglik(data,knots,scale,timescale,z,x,basis) result(ll)
    type(flexsurv_data),intent(in)::data
    real(dp),intent(in)::knots(:),z(:)
    integer,intent(in)::scale,timescale
    real(dp),intent(in),optional::x(:,:)
    integer,intent(in),optional::basis
    type(survspline_model)::m
    integer::i,ng,p
    real(dp)::shift,li,fu,fl,pobs
    ng=size(knots);p=size(z)-ng
    allocate(m%knots(ng),m%gamma(ng));m%knots=knots;m%scale=scale;m%timescale=timescale
    if(present(basis))m%basis=basis
    ll=0.0_dp
    do i=1,size(data%lower)
      shift=0.0_dp
      if(p>0.and.present(x))shift=dot_product(x(i,:),z(ng+1:))
      m%gamma=z(1:ng);m%gamma(1)=m%gamma(1)+shift
      pobs=survspline_cdf(m,data%rtrunc(i))-survspline_cdf(m,data%start(i))
      if(pobs<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
      select case(data%status(i))
      case(fs_status_event)
        li=survspline_logpdf(m,data%lower(i))
        if(data%bhazard(i)>0.0_dp)then
          li=li+log(1.0_dp+data%bhazard(i)/max(spline_haz(m,data%lower(i)),tiny(1.0_dp)))
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
      if(.not.ieee_is_finite(li))then;ll=-huge(1.0_dp);return;end if
      ll=ll+data%weights(i)*li
    end do
  end function flexsurvspline_loglik

  real(dp) function spline_predict_survival(res,row,t,x) result(s)
    type(flexsurvspline_result),intent(in)::res
    integer,intent(in)::row
    real(dp),intent(in)::t
    real(dp),intent(in),optional::x(:,:)
    type(survspline_model)::m
    real(dp)::shift
    m=res%model;shift=0.0_dp
    if(size(res%beta)>0.and.present(x))shift=dot_product(x(row,:),res%beta)
    m%gamma(1)=m%gamma(1)+shift;s=survspline_survival(m,t)
  end function spline_predict_survival

  real(dp) function spline_predict_hazard(res,row,t,x) result(h)
    type(flexsurvspline_result),intent(in)::res
    integer,intent(in)::row
    real(dp),intent(in)::t
    real(dp),intent(in),optional::x(:,:)
    type(survspline_model)::m
    real(dp)::shift
    m=res%model;shift=0.0_dp
    if(size(res%beta)>0.and.present(x))shift=dot_product(x(row,:),res%beta)
    m%gamma(1)=m%gamma(1)+shift;h=spline_haz(m,t)
  end function spline_predict_hazard

  real(dp) function spline_haz(m,t) result(h)
    use flexsurv_spline, only : survspline_hazard
    type(survspline_model),intent(in)::m
    real(dp),intent(in)::t
    h=survspline_hazard(m,t)
  end function spline_haz

  subroutine spline_qp_initial(y,xmat,dxmat,coef,status,eps)
    real(dp),intent(in)::y(:),xmat(:,:),dxmat(:,:)
    real(dp),allocatable,intent(out)::coef(:)
    integer,intent(out)::status
    real(dp),intent(in),optional::eps
    type(qp_result)::qr
    real(dp),allocatable::dmat(:,:),dvec(:),amat(:,:),bvec(:)
    real(dp)::ee
    integer::p
    p=size(xmat,2);ee=1.0e-9_dp;if(present(eps))ee=eps
    allocate(dmat(p,p),dvec(p),amat(p,size(dxmat,1)),bvec(size(dxmat,1)))
    dmat=matmul(transpose(xmat),xmat)
    dmat=dmat+1.0e-10_dp*identity_matrix(p)
    dvec=matmul(transpose(xmat),y);amat=transpose(dxmat);bvec=ee
    qr=solve_qp(dmat,dvec,amat,bvec)
    status=qr%status;coef=qr%solution
  end subroutine spline_qp_initial

  pure function identity_matrix(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0.0_dp;do i=1,n;a(i,i)=1.0_dp;end do
  end function identity_matrix

  subroutine heuristic_gamma(data,knots,g)
    type(flexsurv_data),intent(in)::data
    real(dp),intent(in)::knots(:)
    real(dp),intent(out)::g(:)
    real(dp)::med
    integer::n
    g=0.0_dp;n=size(data%lower)
    med=sum(data%lower)/real(max(n,1),dp)
    if(size(g)>=2)then
      g(1)=0.0_dp;g(2)=1.0_dp
      if(med>0.0_dp.and.knots(2)/=knots(1))g(1)=-log(med)
    end if
  end subroutine heuristic_gamma

  subroutine invert_local(a,ainv,status)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::ainv(:,:)
    integer,intent(out)::status
    real(dp),allocatable::aug(:,:),tmp(:)
    real(dp)::piv,fac
    integer::n,i,j,k,ip
    n=size(a,1);allocate(aug(n,2*n),tmp(2*n));aug=0.0_dp;aug(:,1:n)=a
    do i=1,n;aug(i,n+i)=1.0_dp;end do
    status=0
    do i=1,n
      ip=i
      do k=i+1,n;if(abs(aug(k,i))>abs(aug(ip,i)))ip=k;end do
      if(abs(aug(ip,i))<1.0e-14_dp)then;status=1;allocate(ainv(n,n));ainv=0.0_dp;return;end if
      if(ip/=i)then;tmp=aug(i,:);aug(i,:)=aug(ip,:);aug(ip,:)=tmp;end if
      piv=aug(i,i);aug(i,:)=aug(i,:)/piv
      do j=1,n;if(j/=i)then;fac=aug(j,i);aug(j,:)=aug(j,:)-fac*aug(i,:);end if;end do
    end do
    allocate(ainv(n,n));ainv=aug(:,n+1:)
  end subroutine invert_local

end module flexsurv_spline_fit
