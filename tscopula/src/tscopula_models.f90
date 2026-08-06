! SPDX-License-Identifier: GPL-3.0-only
module tscopula_models
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tscopula_kinds, only : dp, pi
  use tscopula_status, only : tsc_error
  use tscopula_math, only : uniform_random, normal_pdf, empirical_quantile_type7, safe_standard_errors
  use tscopula_margins, only : margin_spec, margin_fit_result, pmarg, qmarg, dmarg, fit_margin
  use tscopula_vtransforms, only : vtransform_spec, vlinear, vtrans, vinverse, vdownprob, vgradient, stochinverse
  use tscopula_paircopula, only : pair_copula, pair_hinv2, pair_cdf
  use tscopula_timeseries, only : arma_copula, sarma_copula, arma_fit_result, sarma2arma, &
    sim_arma_copula, sim_sarma_copula, arma_objective, predict_arma_cdf, &
    predict_arma_quantile, predict_arma_density, resid_arma_copula, fit_arma_copula, kendall_arma
  use tscopula_dvine, only : dvine_copula, dvine_fit_result, sim_dvine, dvine_loglik, &
    predict_dvine_cdf, predict_dvine_quantile, predict_dvine_density, resid_dvine, &
    fit_dvine, kendall_dvine
  implicit none
  private

  integer,parameter,public :: copula_white=0,copula_arma=1,copula_sarma=2,copula_dvine_kind=3

  type,public :: tscopula_spec
    integer :: kind=copula_white
    type(arma_copula) :: arma
    type(sarma_copula) :: sarma
    type(dvine_copula) :: dvine
  end type tscopula_spec

  type,public :: vtscopula_spec
    type(tscopula_spec) :: base
    type(vtransform_spec) :: transform
    logical :: has_wcopula=.false.
    type(pair_copula) :: wcopula
  end type vtscopula_spec

  type,public :: tscm_spec
    type(vtscopula_spec) :: copula
    type(margin_spec) :: margin
  end type tscm_spec

  type,public :: tscm_fit_result
    type(tscm_spec) :: model
    real(dp) :: log_likelihood=-huge(1.0_dp)
    real(dp) :: aic=huge(1.0_dp),bic=huge(1.0_dp)
    integer :: convergence=1
  end type tscm_fit_result

  type,public :: empirical_distribution
    real(dp),allocatable :: data(:)
    real(dp) :: bandwidth=0.0_dp
  end type empirical_distribution

  public :: swncopula, tscopula_from_arma, tscopula_from_sarma, tscopula_from_dvine
  public :: vtscopula, setwcopula, tscm, sim_tscopula, sim_vtscopula, sim_tscm
  public :: tscopula_loglik, tscopula_objective, tscopula_residuals
  public :: predict_tscopula_cdf, predict_tscopula_quantile, predict_tscopula_density
  public :: kendall_tscopula, vtscopula_objective, wobjective, predict_vtscopula_cdf
  public :: predict_vtscopula_quantile, predict_vtscopula_density, profilefulcrum
  public :: tscm_loglik, predict_tscm_cdf, predict_tscm_quantile, predict_tscm_density
  public :: fit_tscm_steps, fit_full, fit_edf, pedf, predict_empirical, aicc, safe_ses

contains

  function swncopula() result(model)
    type(tscopula_spec)::model;model%kind=copula_white
  end function swncopula
  function tscopula_from_arma(arma) result(model)
    type(arma_copula),intent(in)::arma;type(tscopula_spec)::model;model%kind=copula_arma;model%arma=arma
  end function tscopula_from_arma
  function tscopula_from_sarma(sarma) result(model)
    type(sarma_copula),intent(in)::sarma;type(tscopula_spec)::model;model%kind=copula_sarma;model%sarma=sarma
  end function tscopula_from_sarma
  function tscopula_from_dvine(dvine) result(model)
    type(dvine_copula),intent(in)::dvine;type(tscopula_spec)::model;model%kind=copula_dvine_kind;model%dvine=dvine
  end function tscopula_from_dvine

  function vtscopula(base,transform,wcopula) result(model)
    type(tscopula_spec),intent(in)::base;type(vtransform_spec),intent(in),optional::transform
    type(pair_copula),intent(in),optional::wcopula;type(vtscopula_spec)::model
    model%base=base;model%transform=vlinear();if(present(transform))model%transform=transform
    if(present(wcopula))then;model%has_wcopula=.true.;model%wcopula=wcopula;end if
  end function vtscopula

  function setwcopula(model) result(w)
    type(vtscopula_spec),intent(in)::model;type(pair_copula)::w
    w=model%wcopula
  end function setwcopula

  function tscm(copula,margin) result(model)
    type(vtscopula_spec),intent(in)::copula;type(margin_spec),intent(in)::margin;type(tscm_spec)::model
    model%copula=copula;model%margin=margin
  end function tscm

  function sim_tscopula(model,n) result(u)
    type(tscopula_spec),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::u(:);integer::i
    select case(model%kind)
    case(copula_white);allocate(u(n));do i=1,n;u(i)=uniform_random();end do
    case(copula_arma);u=sim_arma_copula(model%arma,n)
    case(copula_sarma);u=sim_sarma_copula(model%sarma,n)
    case(copula_dvine_kind);u=sim_dvine(model%dvine,n)
    case default;allocate(u(n));u=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function sim_tscopula

  function sim_vtscopula(model,n) result(u)
    type(vtscopula_spec),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::u(:),v(:),w(:);integer::i
    v=sim_tscopula(model%base,n)
    if(.not.model%has_wcopula)then;u=stochinverse(model%transform,v);return;end if
    allocate(w(n));w(1)=uniform_random();do i=2,n;w(i)=pair_hinv2(model%wcopula,uniform_random(),w(i-1));end do
    u=stochinverse(model%transform,v,w)
  end function sim_vtscopula

  function sim_tscm(model,n) result(x)
    type(tscm_spec),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:),u(:)
    u=sim_vtscopula(model%copula,n);allocate(x(n));x=qmarg(model%margin,u)
  end function sim_tscm

  real(dp) function tscopula_loglik(model,u) result(value)
    type(tscopula_spec),intent(in)::model;real(dp),intent(in)::u(:);type(arma_copula)::a
    select case(model%kind)
    case(copula_white);value=0.0_dp
    case(copula_arma);value=-arma_objective(model%arma,u)
    case(copula_sarma);a=sarma2arma(model%sarma);value=-arma_objective(a,u)
    case(copula_dvine_kind);value=dvine_loglik(model%dvine,u)
    case default;value=-huge(1.0_dp)
    end select
  end function tscopula_loglik

  real(dp) function tscopula_objective(model,u) result(value)
    type(tscopula_spec),intent(in)::model;real(dp),intent(in)::u(:);value=-tscopula_loglik(model,u)
  end function tscopula_objective

  function tscopula_residuals(model,u) result(r)
    type(tscopula_spec),intent(in)::model;real(dp),intent(in)::u(:);real(dp),allocatable::r(:);type(arma_copula)::a
    select case(model%kind)
    case(copula_white);allocate(r(size(u)));r=u
    case(copula_arma);r=resid_arma_copula(model%arma,u)
    case(copula_sarma);a=sarma2arma(model%sarma);r=resid_arma_copula(a,u)
    case(copula_dvine_kind);r=resid_dvine(model%dvine,u)
    case default;allocate(r(size(u)));r=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function tscopula_residuals

  real(dp) function predict_tscopula_cdf(model,history,x) result(value)
    type(tscopula_spec),intent(in)::model;real(dp),intent(in)::history(:),x;type(arma_copula)::a
    select case(model%kind)
    case(copula_white);value=x
    case(copula_arma);value=predict_arma_cdf(model%arma,history,x)
    case(copula_sarma);a=sarma2arma(model%sarma);value=predict_arma_cdf(a,history,x)
    case(copula_dvine_kind);value=predict_dvine_cdf(model%dvine,history,x)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function predict_tscopula_cdf

  real(dp) function predict_tscopula_quantile(model,history,p) result(value)
    type(tscopula_spec),intent(in)::model;real(dp),intent(in)::history(:),p;type(arma_copula)::a
    select case(model%kind)
    case(copula_white);value=p
    case(copula_arma);value=predict_arma_quantile(model%arma,history,p)
    case(copula_sarma);a=sarma2arma(model%sarma);value=predict_arma_quantile(a,history,p)
    case(copula_dvine_kind);value=predict_dvine_quantile(model%dvine,history,p)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function predict_tscopula_quantile

  real(dp) function predict_tscopula_density(model,history,x) result(value)
    type(tscopula_spec),intent(in)::model;real(dp),intent(in)::history(:),x;type(arma_copula)::a
    select case(model%kind)
    case(copula_white);value=1.0_dp
    case(copula_arma);value=predict_arma_density(model%arma,history,x)
    case(copula_sarma);a=sarma2arma(model%sarma);value=predict_arma_density(a,history,x)
    case(copula_dvine_kind);value=predict_dvine_density(model%dvine,history,x)
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function predict_tscopula_density

  function kendall_tscopula(model,maxlag) result(tau)
    type(tscopula_spec),intent(in)::model;integer,intent(in)::maxlag;real(dp),allocatable::tau(:),tmp(:);type(arma_copula)::a
    allocate(tau(maxlag));tau=0.0_dp
    select case(model%kind)
    case(copula_arma);tau=kendall_arma(model%arma,maxlag)
    case(copula_sarma);a=sarma2arma(model%sarma);tau=kendall_arma(a,maxlag)
    case(copula_dvine_kind);tmp=kendall_dvine(model%dvine);tau(1:min(maxlag,size(tmp)))=tmp(1:min(maxlag,size(tmp)))
    end select
  end function kendall_tscopula

  real(dp) function vtscopula_objective(model,u) result(value)
    type(vtscopula_spec),intent(in)::model;real(dp),intent(in)::u(:);real(dp),allocatable::v(:)
    allocate(v(size(u)));v=vtrans(model%transform,u);value=tscopula_objective(model%base,v)
    if(model%has_wcopula)value=value+wobjective(model%wcopula,model%transform,u)
  end function vtscopula_objective

  real(dp) function wobjective(wcopula,transform,u) result(value)
    type(pair_copula),intent(in)::wcopula;type(vtransform_spec),intent(in)::transform;real(dp),intent(in)::u(:)
    real(dp),allocatable::delta(:);real(dp)::cval,prob,d1,d2;integer::i
    value=0.0_dp;if(size(u)<2)return;allocate(delta(size(u)))
    do i=1,size(u)
      if(transform%family==3)then;delta(i)=transform%delta
      else;delta(i)=-1.0_dp/vgradient(transform,u(i));if(u(i)>transform%delta)delta(i)=delta(i)+1.0_dp;end if
      delta(i)=min(max(delta(i),epsilon(1.0_dp)),1.0_dp-epsilon(1.0_dp))
    end do
    do i=1,size(u)-1
      d1=delta(i);d2=delta(i+1);cval=pair_cdf(wcopula,d1,d2)
      if(u(i)<=transform%delta.and.u(i+1)<=transform%delta)then;prob=cval/(d1*d2)
      else if(u(i)<=transform%delta.and.u(i+1)>transform%delta)then;prob=(d1-cval)/(d1*(1.0_dp-d2))
      else if(u(i)>transform%delta.and.u(i+1)<=transform%delta)then;prob=(d2-cval)/((1.0_dp-d1)*d2)
      else;prob=(1.0_dp-d1-d2+cval)/((1.0_dp-d1)*(1.0_dp-d2));end if
      value=value-log(max(prob,tiny(1.0_dp)))
    end do
  end function wobjective

  real(dp) function predict_vtscopula_density(model,history,x) result(value)
    type(vtscopula_spec),intent(in)::model;real(dp),intent(in)::history(:),x;real(dp),allocatable::vh(:)
    allocate(vh(size(history)));vh=vtrans(model%transform,history)
    value=predict_tscopula_density(model%base,vh,vtrans(model%transform,x))
  end function predict_vtscopula_density

  real(dp) function predict_vtscopula_cdf(model,history,x) result(value)
    type(vtscopula_spec),intent(in)::model;real(dp),intent(in)::history(:),x;real(dp),allocatable::vh(:);real(dp)::d,mult,vx
    allocate(vh(size(history)));vh=vtrans(model%transform,history);d=model%transform%delta
    if(model%transform%family==3.and..not.model%has_wcopula)then
      vx=vtrans(model%transform,x);if(x<=d)then;mult=-d;else;mult=1.0_dp-d;end if
      value=d+mult*predict_tscopula_cdf(model%base,vh,vx)
    else
      value=numeric_cdf(model,history,x)
    end if
    value=min(max(value,0.0_dp),1.0_dp)
  end function predict_vtscopula_cdf

  real(dp) function numeric_cdf(model,history,x) result(value)
    type(vtscopula_spec),intent(in)::model;real(dp),intent(in)::history(:),x
    integer,parameter::nq=400;real(dp)::h,s,t;integer::i
    if(x<=0.0_dp)then;value=0.0_dp;return;else if(x>=1.0_dp)then;value=1.0_dp;return;end if
    h=x/real(nq,dp);s=0.5_dp*(predict_vtscopula_density(model,history,0.0_dp)+predict_vtscopula_density(model,history,x))
    do i=1,nq-1;t=real(i,dp)*h;s=s+predict_vtscopula_density(model,history,t);end do;value=s*h
  end function numeric_cdf

  real(dp) function predict_vtscopula_quantile(model,history,p) result(value)
    type(vtscopula_spec),intent(in)::model;real(dp),intent(in)::history(:),p;real(dp)::lo,hi,mid;integer::i
    lo=0.0_dp;hi=1.0_dp;do i=1,60;mid=0.5_dp*(lo+hi);if(predict_vtscopula_cdf(model,history,mid)<p)then;lo=mid;else;hi=mid;end if;end do;value=0.5_dp*(lo+hi)
  end function predict_vtscopula_quantile

  function profilefulcrum(model,u,grid) result(objective)
    type(vtscopula_spec),intent(in)::model;real(dp),intent(in)::u(:),grid(:);real(dp),allocatable::objective(:);type(vtscopula_spec)::work;integer::i
    allocate(objective(size(grid)));work=model
    do i=1,size(grid);work%transform%delta=grid(i);objective(i)=vtscopula_objective(work,u);end do
  end function profilefulcrum

  real(dp) function tscm_loglik(model,x) result(value)
    type(tscm_spec),intent(in)::model;real(dp),intent(in)::x(:);real(dp),allocatable::u(:);integer::i
    allocate(u(size(x)));u=pmarg(model%margin,x);value=-vtscopula_objective(model%copula,u)
    do i=1,size(x);value=value+dmarg(model%margin,x(i),.true.);end do
  end function tscm_loglik

  real(dp) function predict_tscm_cdf(model,history,x) result(value)
    type(tscm_spec),intent(in)::model;real(dp),intent(in)::history(:),x;real(dp),allocatable::u(:)
    allocate(u(size(history)));u=pmarg(model%margin,history);value=predict_vtscopula_cdf(model%copula,u,pmarg(model%margin,x))
  end function predict_tscm_cdf
  real(dp) function predict_tscm_quantile(model,history,p) result(value)
    type(tscm_spec),intent(in)::model;real(dp),intent(in)::history(:),p;real(dp),allocatable::u(:);real(dp)::up
    allocate(u(size(history)));u=pmarg(model%margin,history);up=predict_vtscopula_quantile(model%copula,u,p);value=qmarg(model%margin,up)
  end function predict_tscm_quantile
  real(dp) function predict_tscm_density(model,history,x) result(value)
    type(tscm_spec),intent(in)::model;real(dp),intent(in)::history(:),x;real(dp),allocatable::u(:);real(dp)::ux
    allocate(u(size(history)));u=pmarg(model%margin,history);ux=pmarg(model%margin,x)
    value=predict_vtscopula_density(model%copula,u,ux)*dmarg(model%margin,x)
  end function predict_tscm_density

  function fit_tscm_steps(template,x,max_iter) result(fit)
    type(tscm_spec),intent(in)::template;real(dp),intent(in)::x(:);integer,intent(in),optional::max_iter
    type(tscm_fit_result)::fit;type(margin_fit_result)::mf;type(tsc_error)::err;real(dp),allocatable::u(:),v(:)
    type(arma_fit_result)::af;type(dvine_fit_result)::df;character(len=12),allocatable::families(:);integer::j,k,npar
    fit%model=template;call fit_margin(template%margin,x,mf,err,max_iter=max_iter)
    if(.not.err%ok())then;fit%convergence=2;return;end if
    fit%model%margin=mf%margin;allocate(u(size(x)));u=pmarg(fit%model%margin,x);allocate(v(size(x)));v=vtrans(fit%model%copula%transform,u)
    select case(fit%model%copula%base%kind)
    case(copula_arma)
      af=fit_arma_copula(v,size(fit%model%copula%base%arma%ar),size(fit%model%copula%base%arma%ma),max_iter=max_iter)
      fit%model%copula%base%arma=af%model;fit%convergence=af%convergence
    case(copula_dvine_kind)
      k=size(fit%model%copula%base%dvine%pairs);allocate(families(k));do j=1,k;families(j)='gauss';end do
      df=fit_dvine(v,families);fit%model%copula%base%dvine=df%model;fit%convergence=df%convergence
    case default;fit%convergence=0
    end select
    fit%log_likelihood=tscm_loglik(fit%model,x);npar=margin_parameter_count(fit%model%margin)+copula_parameter_count(fit%model%copula%base)
    fit%aic=-2.0_dp*fit%log_likelihood+2.0_dp*real(npar,dp);fit%bic=-2.0_dp*fit%log_likelihood+log(real(size(x),dp))*real(npar,dp)
  end function fit_tscm_steps

  function fit_full(template,x,max_iter) result(fit)
    type(tscm_spec),intent(in)::template;real(dp),intent(in)::x(:);integer,intent(in),optional::max_iter;type(tscm_fit_result)::fit
    if(present(max_iter))then;fit=fit_tscm_steps(template,x,max_iter);else;fit=fit_tscm_steps(template,x);end if
  end function fit_full

  integer function margin_parameter_count(m) result(n)
    type(margin_spec),intent(in)::m
    select case(m%family);case(0);n=0;case(1,3);n=2;case(2,4);n=1;case(5);n=3;case(6,8);n=3;case(7,10);n=4;case(9);n=2;case default;n=0;end select
  end function margin_parameter_count
  integer function copula_parameter_count(c) result(n)
    type(tscopula_spec),intent(in)::c;integer::j
    select case(c%kind);case(copula_white);n=0;case(copula_arma);n=size(c%arma%ar)+size(c%arma%ma);case(copula_sarma);n=size(c%sarma%ar)+size(c%sarma%ma)+size(c%sarma%sar)+size(c%sarma%sma)
    case(copula_dvine_kind);n=0;do j=1,size(c%dvine%pairs);n=n+merge(2,1,c%dvine%pairs(j)%family==2.or.c%dvine%pairs(j)%family==7);end do;case default;n=0;end select
  end function copula_parameter_count

  function fit_edf(x,bandwidth) result(model)
    real(dp),intent(in)::x(:);real(dp),intent(in),optional::bandwidth;type(empirical_distribution)::model
    real(dp)::sd,m;allocate(model%data(size(x)));model%data=x;m=sum(x)/real(size(x),dp);sd=sqrt(sum((x-m)**2)/max(real(size(x)-1,dp),1.0_dp))
    model%bandwidth=1.06_dp*sd*real(size(x),dp)**(-0.2_dp);if(present(bandwidth))model%bandwidth=bandwidth
  end function fit_edf

  real(dp) function pedf(model,x) result(value)
    type(empirical_distribution),intent(in)::model;real(dp),intent(in)::x;integer::i
    value=0.0_dp;do i=1,size(model%data);if(model%data(i)<=x)value=value+1.0_dp;end do;value=value/real(size(model%data),dp)
  end function pedf

  real(dp) function predict_empirical(model,x,type) result(value)
    type(empirical_distribution),intent(in)::model;real(dp),intent(in)::x;character(len=*),intent(in),optional::type
    character(len=8)::kind;integer::i;real(dp)::h
    kind='cdf';if(present(type))kind=type
    if(trim(kind)=='quantile')then;value=empirical_quantile_type7(model%data,x)
    else if(trim(kind)=='density')then;h=max(model%bandwidth,sqrt(epsilon(1.0_dp)));value=0.0_dp;do i=1,size(model%data);value=value+normal_pdf((x-model%data(i))/h)/h;end do;value=value/real(size(model%data),dp)
    else;value=pedf(model,x);end if
  end function predict_empirical

  real(dp) function aicc(loglik,k,n) result(value)
    real(dp),intent(in)::loglik;integer,intent(in)::k,n
    value=-2.0_dp*loglik+2.0_dp*real(k,dp)
    if(n>k+1)value=value+2.0_dp*real(k*(k+1),dp)/real(n-k-1,dp)
  end function aicc

  subroutine safe_ses(hessian,se)
    real(dp),intent(in)::hessian(:,:);real(dp),allocatable,intent(out)::se(:);call safe_standard_errors(hessian,se)
  end subroutine safe_ses

end module tscopula_models
