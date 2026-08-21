! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_fit
  use flexsurv_kinds, only : dp
  use flexsurv_distributions, only : dist_npar, dist_valid, dist_logpdf, &
    dist_cdf, dist_survival, dist_hazard, dist_name, &
    dist_exponential, dist_weibull, dist_weibull_ph, dist_gamma, &
    dist_lognormal, dist_gompertz, dist_loglogistic, dist_gengamma, dist_genf
  use flexsurv_math, only : logdiffexp, near_positive_definite
  use numderiv, only : hessian
  use survival_aft, only : survreg_fit
  use survival_types, only : aft_result
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
    ieee_positive_inf, ieee_quiet_nan
  implicit none
  private

  integer, parameter, public :: fs_status_censored = 0
  integer, parameter, public :: fs_status_event = 1
  integer, parameter, public :: fs_status_left = 2
  integer, parameter, public :: fs_status_interval = 3
  integer, parameter, public :: fs_optim_bfgs = 1
  integer, parameter, public :: fs_optim_nelder_mead = 2
  integer, parameter, public :: fs_optim_auto = 3

  type, public :: flexsurv_data
    real(dp), allocatable :: start(:)
    real(dp), allocatable :: lower(:)
    real(dp), allocatable :: upper(:)
    integer, allocatable :: status(:)
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: rtrunc(:)
    real(dp), allocatable :: bhazard(:)
    real(dp), allocatable :: bcondsurv(:)
  end type flexsurv_data

  type, public :: parameter_regression
    real(dp), allocatable :: x(:,:)
  end type parameter_regression

  type, public :: flexsurv_spec
    integer :: dist = dist_weibull
    real(dp), allocatable :: base_init(:)
    type(parameter_regression), allocatable :: reg(:)
  end type flexsurv_spec

  type, public :: flexsurv_result
    integer :: dist = 0
    real(dp), allocatable :: theta(:)
    real(dp), allocatable :: base(:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: covariance_natural(:,:)
    real(dp), allocatable :: gradient(:)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: nobs = 0
    integer :: npar = 0
    integer :: iterations = 0
    integer :: status = 1
    logical :: converged = .false.
    character(len=160) :: message = ''
  end type flexsurv_result

  public :: prepare_survival_data, validate_data
  public :: initialize_spec, parameter_count, location_parameter
  public :: fit_flexsurvreg, flexsurv_loglik
  public :: parameter_row, parameter_row_tcov, baseline_parameters
  public :: predict_survival, predict_hazard, predict_cumhaz, predict_density
  public :: predict_quantile
  public :: aicc, model_bic
  public :: bfgs_minimize, nelder_mead_minimize, objective_fn, initial_theta
  public :: source_initial_theta

  abstract interface
    function objective_fn(x) result(v)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: v
    end function objective_fn
  end interface

contains

  subroutine prepare_survival_data(data, time, status, start, upper, weights, &
      rtrunc, bhazard, bcondsurv)
    type(flexsurv_data), intent(out) :: data
    real(dp), intent(in) :: time(:)
    integer, intent(in), optional :: status(:)
    real(dp), intent(in), optional :: start(:), upper(:), weights(:), rtrunc(:)
    real(dp), intent(in), optional :: bhazard(:), bcondsurv(:)
    integer :: n
    n=size(time)
    allocate(data%start(n),data%lower(n),data%upper(n),data%status(n), &
      data%weights(n),data%rtrunc(n),data%bhazard(n),data%bcondsurv(n))
    data%start=0.0_dp;data%lower=time;data%upper=ieee_value(0.0_dp,ieee_positive_inf)
    data%status=fs_status_censored;data%weights=1.0_dp
    data%rtrunc=ieee_value(0.0_dp,ieee_positive_inf)
    data%bhazard=0.0_dp;data%bcondsurv=1.0_dp
    if(present(status))data%status=status
    if(present(start))data%start=start
    if(present(upper))data%upper=upper
    if(present(weights))data%weights=weights
    if(present(rtrunc))data%rtrunc=rtrunc
    if(present(bhazard))data%bhazard=bhazard
    if(present(bcondsurv))data%bcondsurv=bcondsurv
    where(data%status==fs_status_event)data%upper=data%lower
    where(data%status==fs_status_left)data%lower=0.0_dp
  end subroutine prepare_survival_data

  logical function validate_data(data,message) result(ok)
    type(flexsurv_data),intent(in)::data
    character(len=*),intent(out),optional::message
    integer::n
    ok=.false.;if(present(message))message=''
    if(.not.allocated(data%lower).or..not.allocated(data%status))then
      if(present(message))message='lower and status are required';return
    end if
    n=size(data%lower)
    if(n==0)then;if(present(message))message='empty data';return;end if
    if(.not.allocated(data%upper).or.size(data%upper)/=n)then
      if(present(message))message='upper has wrong size';return
    end if
    if(.not.allocated(data%start).or.size(data%start)/=n)then
      if(present(message))message='start has wrong size';return
    end if
    if(any(data%start<0.0_dp).or.any(data%lower<0.0_dp))then
      if(present(message))message='times must be nonnegative';return
    end if
    if(any(data%upper<data%lower))then
      if(present(message))message='upper must be >= lower';return
    end if
    if(any(data%start>data%lower))then
      if(present(message))message='start must be <= lower';return
    end if
    ok=.true.
  end function validate_data

  subroutine initialize_spec(spec,dist,nobs,base_init)
    type(flexsurv_spec),intent(out)::spec
    integer,intent(in)::dist,nobs
    real(dp),intent(in),optional::base_init(:)
    integer::nb,i
    spec%dist=dist;nb=dist_npar(dist)
    allocate(spec%base_init(nb),spec%reg(nb))
    spec%base_init=default_base(dist)
    if(present(base_init))then
      if(size(base_init)==nb)spec%base_init=base_init
    end if
    do i=1,nb;allocate(spec%reg(i)%x(nobs,0));end do
  end subroutine initialize_spec

  integer function parameter_count(spec) result(n)
    type(flexsurv_spec),intent(in)::spec
    integer::i
    n=size(spec%base_init)
    if(allocated(spec%reg))then
      do i=1,size(spec%reg)
        if(allocated(spec%reg(i)%x))n=n+size(spec%reg(i)%x,2)
      end do
    end if
  end function parameter_count

  integer function location_parameter(dist) result(i)
    integer,intent(in)::dist
    select case(dist)
    case(dist_exponential);i=1
    case(dist_weibull,dist_weibull_ph);i=2
    case(dist_gamma);i=2
    case(dist_lognormal);i=1
    case(dist_gompertz);i=2
    case(dist_loglogistic);i=2
    case(dist_gengamma,dist_genf);i=1
    case default;i=1
    end select
  end function location_parameter

  function fit_flexsurvreg(data,spec,fixed,control_maxit,tol,theta_init, &
      source_start,optim_method,fallback) result(res)
    type(flexsurv_data),intent(in)::data
    type(flexsurv_spec),intent(in)::spec
    logical,intent(in),optional::fixed(:)
    integer,intent(in),optional::control_maxit,optim_method
    real(dp),intent(in),optional::tol,theta_init(:)
    logical,intent(in),optional::source_start,fallback
    type(flexsurv_result)::res
    real(dp),allocatable::theta0(:),free0(:),freehat(:),full(:),hess(:,:), &
      hess_pd(:,:),covfree(:,:),jac(:,:),covnat(:,:),g(:)
    logical,allocatable::fix(:)
    integer,allocatable::freeidx(:)
    integer::npar,nfree,i,j,stat,maxit,nobs,omethod,it2,st2
    real(dp)::ftol,fval,fval2
    logical::do_source,do_fallback
    real(dp),allocatable::x2(:),g2(:)
    character(len=:),allocatable::ndmsg
    character(len=160)::msg

    nobs=size(data%lower);npar=parameter_count(spec)
    res%dist=spec%dist;res%nobs=nobs;res%npar=npar
    if(.not.validate_data(data,msg))then;res%message=msg;return;end if
    if(.not.validate_spec_rows(spec,nobs,msg))then;res%message=msg;return;end if
    do_source=.false.;if(present(source_start))do_source=source_start
    if(do_source)then
      theta0=source_initial_theta(data,spec)
    else
      theta0=initial_theta(spec)
    end if
    if(present(theta_init))then
      if(size(theta_init)==npar)theta0=theta_init
    end if
    allocate(fix(npar));fix=.false.;if(present(fixed))then;if(size(fixed)==npar)fix=fixed;end if
    nfree=count(.not.fix);allocate(freeidx(nfree));j=0
    do i=1,npar;if(.not.fix(i))then;j=j+1;freeidx(j)=i;end if;end do
    if(nfree==0)then
      full=theta0;fval=objective_full(full);res%theta=full;res%base=baseline_parameters(spec,full)
      res%loglik=-fval;res%aic=-2.0_dp*res%loglik;res%bic=-2.0_dp*res%loglik
      res%converged=.true.;res%status=0;res%message='all parameters fixed';return
    end if
    allocate(free0(nfree));free0=theta0(freeidx)
    maxit=500;if(present(control_maxit))maxit=control_maxit
    ftol=1.0e-8_dp;if(present(tol))ftol=tol
    omethod=fs_optim_bfgs;if(present(optim_method))omethod=optim_method
    do_fallback=.true.;if(present(fallback))do_fallback=fallback
    select case(omethod)
    case(fs_optim_nelder_mead)
      call nelder_mead_minimize(objective_free,free0,freehat,fval,res%iterations,stat,maxit,ftol)
      allocate(g(size(freehat)));call finite_gradient(objective_free,freehat,g)
    case default
      call bfgs_minimize(objective_free,free0,freehat,fval,g,res%iterations,stat,maxit,ftol)
      if((stat/=0.or..not.ieee_is_finite(fval)).and.do_fallback)then
        call nelder_mead_minimize(objective_free,free0,x2,fval2,it2,st2,maxit,ftol)
        if(ieee_is_finite(fval2).and.(.not.ieee_is_finite(fval).or.fval2<fval))then
          freehat=x2;fval=fval2;res%iterations=res%iterations+it2;stat=st2
          allocate(g2(size(freehat)));call finite_gradient(objective_free,freehat,g2);g=g2
        end if
      end if
    end select
    full=theta0;full(freeidx)=freehat
    res%theta=full;allocate(res%gradient(npar));res%gradient=0.0_dp;res%gradient(freeidx)=g
    res%base=baseline_parameters(spec,full);res%loglik=-fval
    res%aic=-2.0_dp*res%loglik+2.0_dp*real(nfree,dp)
    res%bic=-2.0_dp*res%loglik+log(real(max(1,nobs),dp))*real(nfree,dp)
    res%converged=(stat==0);res%status=stat
    if(res%converged)then;res%message='success';else;res%message='optimizer did not converge';end if

    call hessian(objective_free,freehat,hess,status=stat,message=ndmsg)
    allocate(res%covariance(npar,npar));res%covariance=0.0_dp
    if(stat==0.and.all(ieee_is_finite(hess)))then
      allocate(hess_pd(nfree,nfree));call near_positive_definite(hess,hess_pd,1.0e-9_dp)
      call invert_matrix(hess_pd,covfree,stat)
      if(stat==0)then
        do i=1,nfree;do j=1,nfree
          res%covariance(freeidx(i),freeidx(j))=covfree(i,j)
        end do;end do
      end if
    end if
    jac=natural_jacobian(spec,full)
    allocate(covnat(size(jac,1),size(jac,1)));covnat=matmul(jac,matmul(res%covariance,transpose(jac)))
    res%covariance_natural=covnat
  contains
    real(dp) function objective_free(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp),allocatable::z(:)
      z=theta0;z(freeidx)=x;v=objective_full(z)
    end function objective_free
    real(dp) function objective_full(x) result(v)
      real(dp),intent(in)::x(:)
      v=-flexsurv_loglik(data,spec,x)
      if(.not.ieee_is_finite(v))v=1.0e100_dp
    end function objective_full
  end function fit_flexsurvreg

  real(dp) function flexsurv_loglik(data,spec,theta,individual) result(ll)
    type(flexsurv_data),intent(in)::data
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:)
    real(dp),intent(out),optional::individual(:)
    integer::i,n
    real(dp),allocatable::par(:)
    real(dp)::li,logpobs,logfu,logfl,pobs,hx,lo,up,bg
    n=size(data%lower);ll=0.0_dp
    if(present(individual))individual=0.0_dp
    do i=1,n
      par=parameter_row(spec,theta,i)
      if(.not.dist_valid(spec%dist,par))then;ll=-huge(1.0_dp);return;end if
      pobs=dist_cdf(spec%dist,data%rtrunc(i),par)-dist_cdf(spec%dist,data%start(i),par)
      if(pobs<=0.0_dp.or..not.ieee_is_finite(pobs))then;ll=-huge(1.0_dp);return;end if
      logpobs=log(pobs)
      select case(data%status(i))
      case(fs_status_event)
        li=dist_logpdf(spec%dist,data%lower(i),par)
        if(data%bhazard(i)>0.0_dp)then
          hx=dist_hazard(spec%dist,data%lower(i),par)
          if(hx<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
          li=li+log(1.0_dp+data%bhazard(i)/hx)
        end if
      case default
        lo=max(0.0_dp,data%lower(i));up=data%upper(i)
        if(ieee_is_finite(up))then;logfu=log(max(dist_cdf(spec%dist,up,par),tiny(1.0_dp)))
        else;logfu=0.0_dp;end if
        if(lo>0.0_dp)then;logfl=log(max(dist_cdf(spec%dist,lo,par),tiny(1.0_dp)))
        else;logfl=-ieee_value(0.0_dp,ieee_positive_inf);end if
        if(data%bhazard(i)>0.0_dp)then
          bg=max(0.0_dp,min(1.0_dp,data%bcondsurv(i)))
          li=log(max((exp(logfu)-1.0_dp)*bg+1.0_dp-exp(logfl),tiny(1.0_dp)))
        else
          li=logdiffexp(logfu,logfl)
        end if
      end select
      li=li-logpobs
      if(.not.ieee_is_finite(li))then;ll=-huge(1.0_dp);return;end if
      ll=ll+data%weights(i)*li
      if(present(individual))individual(i)=li
    end do
  end function flexsurv_loglik

  function parameter_row(spec,theta,row) result(par)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:)
    integer,intent(in)::row
    real(dp),allocatable::par(:)
    real(dp),allocatable::z(:)
    integer::nb,i,j,k,p
    nb=size(spec%base_init);allocate(z(nb),par(nb));z=theta(1:nb);k=nb
    do i=1,nb
      if(allocated(spec%reg(i)%x))then
        p=size(spec%reg(i)%x,2)
        if(p>0)then
          do j=1,p;z(i)=z(i)+spec%reg(i)%x(row,j)*theta(k+j);end do
          k=k+p
        end if
      end if
      par(i)=inverse_transform(spec%dist,i,z(i))
    end do
  end function parameter_row

  function parameter_row_tcov(spec,theta,row,current_time,param_idx,col_idx,rate) result(par)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),current_time
    integer,intent(in)::row,param_idx(:),col_idx(:)
    real(dp),intent(in)::rate(:)
    real(dp),allocatable::par(:)
    real(dp),allocatable::z(:)
    real(dp)::xij
    integer::nb,i,j,k,p,r
    nb=size(spec%base_init);allocate(z(nb),par(nb));z=theta(1:nb);k=nb
    do i=1,nb
      if(allocated(spec%reg(i)%x))then
        p=size(spec%reg(i)%x,2)
        if(p>0)then
          do j=1,p
            xij=spec%reg(i)%x(row,j)
            do r=1,min(size(param_idx),min(size(col_idx),size(rate)))
              if(param_idx(r)==i.and.col_idx(r)==j) xij=xij+rate(r)*current_time
            end do
            z(i)=z(i)+xij*theta(k+j)
          end do
          k=k+p
        end if
      end if
      par(i)=inverse_transform(spec%dist,i,z(i))
    end do
  end function parameter_row_tcov

  function baseline_parameters(spec,theta) result(par)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:)
    real(dp),allocatable::par(:)
    integer::i,nb
    nb=size(spec%base_init);allocate(par(nb))
    do i=1,nb;par(i)=inverse_transform(spec%dist,i,theta(i));end do
  end function baseline_parameters

  real(dp) function predict_survival(spec,theta,row,t) result(v)
    type(flexsurv_spec),intent(in)::spec;real(dp),intent(in)::theta(:),t;integer,intent(in)::row
    real(dp),allocatable::p(:);p=parameter_row(spec,theta,row);v=dist_survival(spec%dist,t,p)
  end function predict_survival
  real(dp) function predict_hazard(spec,theta,row,t) result(v)
    type(flexsurv_spec),intent(in)::spec;real(dp),intent(in)::theta(:),t;integer,intent(in)::row
    real(dp),allocatable::p(:);p=parameter_row(spec,theta,row);v=dist_hazard(spec%dist,t,p)
  end function predict_hazard
  real(dp) function predict_cumhaz(spec,theta,row,t) result(v)
    type(flexsurv_spec),intent(in)::spec;real(dp),intent(in)::theta(:),t;integer,intent(in)::row
    real(dp),allocatable::p(:);p=parameter_row(spec,theta,row);v=-log(dist_survival(spec%dist,t,p))
  end function predict_cumhaz
  real(dp) function predict_density(spec,theta,row,t) result(v)
    type(flexsurv_spec),intent(in)::spec;real(dp),intent(in)::theta(:),t;integer,intent(in)::row
    real(dp),allocatable::p(:);p=parameter_row(spec,theta,row);v=exp(dist_logpdf(spec%dist,t,p))
  end function predict_density
  real(dp) function predict_quantile(spec,theta,row,prob) result(v)
    use flexsurv_distributions, only : dist_quantile
    type(flexsurv_spec),intent(in)::spec;real(dp),intent(in)::theta(:),prob;integer,intent(in)::row
    real(dp),allocatable::p(:);p=parameter_row(spec,theta,row);v=dist_quantile(spec%dist,prob,p)
  end function predict_quantile

  pure real(dp) function aicc(loglik,k,n) result(v)
    real(dp),intent(in)::loglik;integer,intent(in)::k,n
    v=-2.0_dp*loglik+2.0_dp*real(k,dp)
    if(n>k+1)v=v+2.0_dp*real(k*(k+1),dp)/real(n-k-1,dp)
  end function aicc
  pure real(dp) function model_bic(loglik,k,n) result(v)
    real(dp),intent(in)::loglik;integer,intent(in)::k,n
    v=-2.0_dp*loglik+log(real(n,dp))*real(k,dp)
  end function model_bic

  function source_initial_theta(data,spec) result(theta)
    type(flexsurv_data),intent(in)::data
    type(flexsurv_spec),intent(in)::spec
    real(dp),allocatable::theta(:),base(:),tt(:),xloc(:,:)
    real(dp)::sw,m,v,sd,shape,scale,q25,med
    integer::i,n,nb,loc,p,k,j
    logical::aft_ok
    type(aft_result)::aft
    n=size(data%lower);nb=size(spec%base_init);allocate(base(nb));base=spec%base_init
    allocate(tt(n));tt=data%lower
    do i=1,n
      if(data%status(i)==fs_status_interval.and.ieee_is_finite(data%upper(i))) &
        tt(i)=0.5_dp*(data%lower(i)+data%upper(i))
    end do
    sw=sum(data%weights)
    if(sw>0.0_dp)tt=tt*data%weights*real(n,dp)/sw
    aft_ok=all(data%start==0.0_dp).and.all(data%status==fs_status_event.or. &
      data%status==fs_status_censored).and.all(tt>0.0_dp)
    loc=location_parameter(spec%dist);p=0
    if(allocated(spec%reg(loc)%x))p=size(spec%reg(loc)%x,2)
    if(aft_ok.and.any(spec%dist==[dist_exponential,dist_weibull,dist_weibull_ph, &
        dist_lognormal,dist_loglogistic]))then
      allocate(xloc(n,p+1));xloc(:,1)=1.0_dp
      if(p>0)xloc(:,2:)=spec%reg(loc)%x
      select case(spec%dist)
      case(dist_exponential)
        call survreg_fit(tt,data%status,xloc,'exponential',aft,weights=data%weights)
        if(aft%converged)base(1)=exp(-aft%coef(1))
      case(dist_weibull,dist_weibull_ph)
        call survreg_fit(tt,data%status,xloc,'weibull',aft,weights=data%weights)
        if(aft%converged)then
          shape=1.0_dp/max(aft%scale,tiny(1.0_dp));scale=exp(aft%coef(1));base(1)=shape
          if(spec%dist==dist_weibull)then;base(2)=scale;else;base(2)=scale**(-shape);end if
        end if
      case(dist_lognormal)
        call survreg_fit(tt,data%status,xloc,'lognormal',aft,weights=data%weights)
        if(aft%converged)then;base(1)=aft%coef(1);base(2)=aft%scale;end if
      case(dist_loglogistic)
        call survreg_fit(tt,data%status,xloc,'loglogistic',aft,weights=data%weights)
        if(aft%converged)then;base(1)=1.0_dp/max(aft%scale,tiny(1.0_dp));base(2)=exp(aft%coef(1));end if
      end select
    else
      call positive_moments(tt,m,v,sd,med,q25)
      select case(spec%dist)
      case(dist_exponential);base(1)=1.0_dp/max(m,tiny(1.0_dp))
      case(dist_weibull,dist_weibull_ph)
        call log_moments(tt,m,v,sd);shape=max(1.0e-6_dp,1.64_dp/max(v,1.0e-8_dp))
        scale=exp(m+0.572_dp);base(1)=shape
        if(spec%dist==dist_weibull)then;base(2)=scale;else;base(2)=scale**(-shape);end if
      case(dist_lognormal)
        call log_moments(tt,m,v,sd);base(1)=m;base(2)=max(sd,1.0e-6_dp)
      case(dist_loglogistic)
        call positive_moments(tt,m,v,sd,med,q25);scale=max(med,1.0e-8_dp)
        if(q25>0.0_dp.and.q25/=scale)then
          shape=1.0_dp/(log(q25/scale)/log(3.0_dp));if(shape<0.0_dp)shape=1.0_dp
        else;shape=1.0_dp;end if
        base(1)=shape;base(2)=scale
      case(dist_gamma)
        call positive_moments(tt,m,v,sd,med,q25)
        base(1)=max(m*m/max(v,1.0e-8_dp),1.0e-6_dp);base(2)=max(m/max(v,1.0e-8_dp),1.0e-6_dp)
      case(dist_gompertz)
        call positive_moments(tt,m,v,sd,med,q25);base(1)=0.001_dp;base(2)=1.0_dp/max(m,tiny(1.0_dp))
      case(dist_gengamma)
        call log_moments(tt,m,v,sd);base(1)=m;base(2)=max(sd,1.0e-6_dp);base(3)=0.0_dp
      case(dist_genf)
        call log_moments(tt,m,v,sd);base(1)=m;base(2)=max(sd,1.0e-6_dp);base(3)=0.0_dp;base(4)=1.0_dp
      end select
    end if
    theta=initial_theta(spec)
    do i=1,nb;theta(i)=forward_transform(spec%dist,i,base(i));end do
    if(aft_ok.and.aft%converged.and.p>0)then
      k=nb
      do i=1,nb
        if(allocated(spec%reg(i)%x))then
          j=size(spec%reg(i)%x,2)
          if(i==loc.and.j==p)then
            select case(spec%dist)
            case(dist_exponential);theta(k+1:k+p)=-aft%coef(2:p+1)
            case(dist_weibull);theta(k+1:k+p)=aft%coef(2:p+1)
            case(dist_weibull_ph);theta(k+1:k+p)=-aft%coef(2:p+1)*base(1)
            case(dist_lognormal,dist_loglogistic);theta(k+1:k+p)=aft%coef(2:p+1)
            end select
          end if
          k=k+j
        end if
      end do
    end if
  end function source_initial_theta

  subroutine positive_moments(t,m,v,sd,med,q25)
    real(dp),intent(in)::t(:)
    real(dp),intent(out)::m,v,sd,med,q25
    real(dp),allocatable::z(:)
    integer::n
    z=pack(t,t>0.0_dp);n=size(z)
    if(n==0)then;m=1.0_dp;v=1.0_dp;sd=1.0_dp;med=1.0_dp;q25=0.5_dp;return;end if
    call sort_real(z);m=sum(z)/real(n,dp)
    if(n>1)then;v=sum((z-m)**2)/real(n-1,dp);else;v=max(m*m,1.0_dp);end if
    sd=sqrt(max(v,0.0_dp));med=sample_quantile_sorted(z,0.5_dp);q25=sample_quantile_sorted(z,0.25_dp)
  end subroutine positive_moments

  subroutine log_moments(t,m,v,sd)
    real(dp),intent(in)::t(:)
    real(dp),intent(out)::m,v,sd
    real(dp),allocatable::z(:)
    integer::n
    z=log(pack(t,t>0.0_dp));n=size(z)
    if(n==0)then;m=0.0_dp;v=1.0_dp;sd=1.0_dp;return;end if
    m=sum(z)/real(n,dp)
    if(n>1)then;v=sum((z-m)**2)/real(n-1,dp);else;v=1.0_dp;end if
    sd=sqrt(max(v,0.0_dp))
  end subroutine log_moments

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    real(dp)::key
    integer::i,j
    do i=2,size(x)
      key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  real(dp) function sample_quantile_sorted(x,p) result(q)
    real(dp),intent(in)::x(:),p
    real(dp)::h,a
    integer::lo,hi,n
    n=size(x);if(n==1)then;q=x(1);return;end if
    h=1.0_dp+(real(n-1,dp))*p;lo=floor(h);hi=ceiling(h);a=h-real(lo,dp)
    q=(1.0_dp-a)*x(max(1,lo))+a*x(min(n,hi))
  end function sample_quantile_sorted

  function initial_theta(spec) result(theta)
    type(flexsurv_spec),intent(in)::spec
    real(dp),allocatable::theta(:)
    integer::nb,n,i,k,p
    nb=size(spec%base_init);n=parameter_count(spec);allocate(theta(n));theta=0.0_dp
    do i=1,nb;theta(i)=forward_transform(spec%dist,i,spec%base_init(i));end do
    k=nb
    do i=1,nb
      if(allocated(spec%reg(i)%x))then;p=size(spec%reg(i)%x,2);if(p>0)theta(k+1:k+p)=0.0_dp;k=k+p;end if
    end do
  end function initial_theta

  pure real(dp) function forward_transform(dist,idx,x) result(z)
    integer,intent(in)::dist,idx;real(dp),intent(in)::x
    logical::positive
    positive=positive_parameter(dist,idx)
    if(positive)then;z=log(max(x,1.0e-12_dp));else;z=x;end if
  end function forward_transform
  pure real(dp) function inverse_transform(dist,idx,z) result(x)
    integer,intent(in)::dist,idx;real(dp),intent(in)::z
    if(positive_parameter(dist,idx))then
      if(z>700.0_dp)then;x=huge(1.0_dp);else if(z<-700.0_dp)then;x=tiny(1.0_dp);else;x=exp(z);end if
    else;x=z;end if
  end function inverse_transform
  pure logical function positive_parameter(dist,idx) result(pos)
    integer,intent(in)::dist,idx
    pos=.false.
    select case(dist)
    case(dist_exponential);pos=(idx==1)
    case(dist_weibull,dist_weibull_ph);pos=.true.
    case(dist_gamma);pos=.true.
    case(dist_lognormal);pos=(idx==2)
    case(dist_gompertz);pos=(idx==2)
    case(dist_loglogistic);pos=.true.
    case(dist_gengamma);pos=(idx==2)
    case(dist_genf);pos=(idx==2.or.idx==4)
    end select
  end function positive_parameter

  function default_base(dist) result(p)
    integer,intent(in)::dist
    real(dp),allocatable::p(:)
    allocate(p(dist_npar(dist)))
    select case(dist)
    case(dist_exponential);p=[1.0_dp]
    case(dist_weibull);p=[1.0_dp,1.0_dp]
    case(dist_weibull_ph);p=[1.0_dp,1.0_dp]
    case(dist_gamma);p=[1.0_dp,1.0_dp]
    case(dist_lognormal);p=[0.0_dp,1.0_dp]
    case(dist_gompertz);p=[0.001_dp,1.0_dp]
    case(dist_loglogistic);p=[1.0_dp,1.0_dp]
    case(dist_gengamma);p=[0.0_dp,1.0_dp,0.0_dp]
    case(dist_genf);p=[0.0_dp,1.0_dp,0.0_dp,1.0_dp]
    end select
  end function default_base

  logical function validate_spec_rows(spec,n,msg) result(ok)
    type(flexsurv_spec),intent(in)::spec;integer,intent(in)::n;character(len=*),intent(out)::msg
    integer::i
    ok=.false.;msg=''
    if(size(spec%base_init)/=dist_npar(spec%dist))then;msg='base_init size does not match distribution';return;end if
    if(.not.dist_valid(spec%dist,spec%base_init))then;msg='invalid base_init';return;end if
    if(.not.allocated(spec%reg))then;msg='regression array is not allocated';return;end if
    if(size(spec%reg)/=size(spec%base_init))then;msg='regression array size mismatch';return;end if
    do i=1,size(spec%reg)
      if(allocated(spec%reg(i)%x))then
        if(size(spec%reg(i)%x,1)/=n)then;msg='covariate matrix row count mismatch';return;end if
      end if
    end do
    ok=.true.
  end function validate_spec_rows

  function natural_jacobian(spec,theta) result(jac)
    type(flexsurv_spec),intent(in)::spec;real(dp),intent(in)::theta(:)
    real(dp),allocatable::jac(:,:)
    real(dp),allocatable::p(:)
    integer::i,n
    n=size(theta);allocate(jac(n,n));jac=0.0_dp
    p=baseline_parameters(spec,theta)
    do i=1,size(p)
      if(positive_parameter(spec%dist,i))then;jac(i,i)=p(i);else;jac(i,i)=1.0_dp;end if
    end do
    do i=size(p)+1,n;jac(i,i)=1.0_dp;end do
  end function natural_jacobian

  subroutine bfgs_minimize(fn,x0,x,fval,gradv,iters,status,maxit,tol)
    procedure(objective_fn)::fn
    real(dp),intent(in)::x0(:)
    real(dp),allocatable,intent(out)::x(:),gradv(:)
    real(dp),intent(out)::fval
    integer,intent(out)::iters,status
    integer,intent(in)::maxit
    real(dp),intent(in)::tol
    real(dp),allocatable::h(:,:),g(:),gnew(:),p(:),s(:),y(:),xnew(:),hy(:)
    real(dp)::fnew,alpha,rho,ys,gg
    integer::n,i
    n=size(x0);allocate(x(n),gradv(n),h(n,n),g(n),gnew(n),p(n),s(n),y(n),xnew(n),hy(n))
    x=x0;h=0.0_dp;do i=1,n;h(i,i)=1.0_dp;end do
    fval=fn(x);call finite_gradient(fn,x,g)
    status=1
    do iters=1,maxit
      gg=maxval(abs(g));if(gg<tol)then;status=0;exit;end if
      p=-matmul(h,g)
      if(dot_product(p,g)>=0.0_dp)p=-g
      alpha=1.0_dp
      do i=1,40
        xnew=x+alpha*p;fnew=fn(xnew)
        if(ieee_is_finite(fnew).and.fnew<=fval+1.0e-4_dp*alpha*dot_product(g,p))exit
        alpha=0.5_dp*alpha
      end do
      if(alpha<1.0e-12_dp)exit
      call finite_gradient(fn,xnew,gnew)
      s=xnew-x;y=gnew-g;ys=dot_product(y,s)
      if(ys>1.0e-12_dp*sqrt(max(dot_product(y,y)*dot_product(s,s),tiny(1.0_dp))))then
        rho=1.0_dp/ys;hy=matmul(h,y)
        h=h+(1.0_dp+rho*dot_product(y,hy))*rho*outer(s,s)-rho*(outer(hy,s)+outer(s,hy))
      else
        h=0.0_dp;do i=1,n;h(i,i)=1.0_dp;end do
      end if
      x=xnew;g=gnew;fval=fnew
    end do
    gradv=g
    if(maxval(abs(g))<10.0_dp*tol)status=0
  end subroutine bfgs_minimize

  subroutine nelder_mead_minimize(fn,x0,x,fval,iters,status,maxit,tol)
    procedure(objective_fn)::fn
    real(dp),intent(in)::x0(:)
    real(dp),allocatable,intent(out)::x(:)
    real(dp),intent(out)::fval
    integer,intent(out)::iters,status
    integer,intent(in)::maxit
    real(dp),intent(in)::tol
    real(dp),allocatable::simp(:,:),fv(:),cent(:),xr(:),xe(:),xc(:),tmp(:)
    real(dp)::fr,fe,fc,spread,step
    integer::n,i,j,best,worst,second,ii
    n=size(x0);allocate(simp(n,n+1),fv(n+1),cent(n),xr(n),xe(n),xc(n),tmp(n))
    simp(:,1)=x0
    do j=1,n
      simp(:,j+1)=x0;step=0.05_dp*max(1.0_dp,abs(x0(j)));simp(j,j+1)=simp(j,j+1)+step
    end do
    do j=1,n+1;fv(j)=fn(simp(:,j));end do
    status=1
    do iters=1,maxit
      call order_simplex(fv,best,worst,second)
      spread=maxval(abs(fv-fv(best)))
      if(spread<tol*(1.0_dp+abs(fv(best))))then;status=0;exit;end if
      cent=0.0_dp
      do j=1,n+1;if(j/=worst)cent=cent+simp(:,j);end do
      cent=cent/real(n,dp);xr=cent+(cent-simp(:,worst));fr=fn(xr)
      if(fr<fv(best))then
        xe=cent+2.0_dp*(xr-cent);fe=fn(xe)
        if(fe<fr)then;simp(:,worst)=xe;fv(worst)=fe;else;simp(:,worst)=xr;fv(worst)=fr;end if
      else if(fr<fv(second))then
        simp(:,worst)=xr;fv(worst)=fr
      else
        if(fr<fv(worst))then;xc=cent+0.5_dp*(xr-cent);else;xc=cent+0.5_dp*(simp(:,worst)-cent);end if
        fc=fn(xc)
        if(fc<min(fr,fv(worst)))then
          simp(:,worst)=xc;fv(worst)=fc
        else
          do j=1,n+1
            if(j/=best)then;simp(:,j)=simp(:,best)+0.5_dp*(simp(:,j)-simp(:,best));fv(j)=fn(simp(:,j));end if
          end do
        end if
      end if
    end do
    call order_simplex(fv,best,worst,second);allocate(x(n));x=simp(:,best);fval=fv(best)
  contains
    subroutine order_simplex(v,ib,iw,is)
      real(dp),intent(in)::v(:);integer,intent(out)::ib,iw,is
      integer::a
      ib=1;iw=1
      do a=2,size(v);if(v(a)<v(ib))ib=a;if(v(a)>v(iw))iw=a;end do
      is=ib
      do a=1,size(v)
        if(a/=iw)then
          if(is==iw.or.v(a)>v(is))is=a
        end if
      end do
    end subroutine order_simplex
  end subroutine nelder_mead_minimize

  subroutine finite_gradient(fn,x,g)
    procedure(objective_fn)::fn
    real(dp),intent(in)::x(:);real(dp),intent(out)::g(size(x))
    real(dp),allocatable::xp(:),xm(:)
    real(dp)::h,fp,fm
    integer::i
    allocate(xp(size(x)),xm(size(x)))
    do i=1,size(x)
      h=1.0e-5_dp*max(1.0_dp,abs(x(i)));xp=x;xm=x;xp(i)=xp(i)+h;xm(i)=xm(i)-h
      fp=fn(xp);fm=fn(xm);g(i)=(fp-fm)/(2.0_dp*h)
    end do
  end subroutine finite_gradient

  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:);real(dp)::c(size(a),size(b));integer::i,j
    do i=1,size(a);do j=1,size(b);c(i,j)=a(i)*b(j);end do;end do
  end function outer

  subroutine invert_matrix(a,ainv,status)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::ainv(:,:)
    integer,intent(out)::status
    real(dp),allocatable::aug(:,:),tmp(:)
    real(dp)::pivot,factor
    integer::n,i,j,k,ip
    n=size(a,1);allocate(aug(n,2*n),tmp(2*n));aug=0.0_dp;aug(:,1:n)=a
    do i=1,n;aug(i,n+i)=1.0_dp;end do
    status=0
    do i=1,n
      ip=i
      do k=i+1,n;if(abs(aug(k,i))>abs(aug(ip,i)))ip=k;end do
      if(abs(aug(ip,i))<1.0e-14_dp)then;status=1;allocate(ainv(n,n));ainv=0.0_dp;return;end if
      if(ip/=i)then;tmp=aug(i,:);aug(i,:)=aug(ip,:);aug(ip,:)=tmp;end if
      pivot=aug(i,i);aug(i,:)=aug(i,:)/pivot
      do j=1,n
        if(j/=i)then;factor=aug(j,i);aug(j,:)=aug(j,:)-factor*aug(i,:);end if
      end do
    end do
    allocate(ainv(n,n));ainv=aug(:,n+1:2*n)
  end subroutine invert_matrix

end module flexsurv_fit
