! SPDX-License-Identifier: GPL-3.0-only
module tscopula_dvine
  use tscopula_kinds, only : dp
  use tscopula_math, only : uniform_random, integrate_simpson
  use tscopula_paircopula, only : pair_copula, paircop, pair_density, pair_log_density, &
    pair_h1, pair_h2, pair_kendall, kendall_to_parameter, family_from_name
  use tscopula_timeseries, only : kpacf_arma, kpacf_sarma4, kpacf_sarma12, kpacf_fbn, kpacf_arfima
  implicit none
  private

  type, public :: dvine_copula
    type(pair_copula), allocatable :: pairs(:)
  end type dvine_copula

  type, public :: dvine2_copula
    type(dvine_copula) :: vine
    character(len=16) :: kpacf_model = 'arma'
  end type dvine2_copula

  type, public :: dvine3_copula
    type(dvine_copula) :: vine
  end type dvine3_copula

  type, public :: dvine_fit_result
    type(dvine_copula) :: model
    real(dp) :: objective = huge(1.0_dp)
    real(dp), allocatable :: kendall_parameters(:)
    integer :: convergence = 0
  end type dvine_fit_result

  type :: cond_context
    type(dvine_copula) :: model
    real(dp), allocatable :: history(:)
  end type cond_context

  public :: dvinecopula, dvinecopula2, dvinecopula3, mklist_dvine, mklist_dvine2, mklist_dvine3
  public :: dvine_loglik, dvine_objective, rblatt, irblatt, rblattdens
  public :: simdvine, sim_dvine, resid_dvine, predict_dvine_cdf
  public :: predict_dvine_quantile, predict_dvine_density, kendall_dvine
  public :: fit_dvine, arma2dvine, sarma2dvine

contains

  function dvinecopula(pairs) result(model)
    type(pair_copula),intent(in)::pairs(:);type(dvine_copula)::model
    allocate(model%pairs(size(pairs)));model%pairs=pairs
  end function dvinecopula

  function dvinecopula2(family,maxlag,kpacf_model,ar,ma,sar,sma,d,hurst,rotation,df) result(model)
    character(len=*),intent(in)::family;integer,intent(in)::maxlag
    character(len=*),intent(in),optional::kpacf_model
    real(dp),intent(in),optional::ar(:),ma(:),sar(:),sma(:),d,hurst,df
    integer,intent(in),optional::rotation
    type(dvine2_copula)::model
    real(dp),allocatable::a(:),b(:),sa(:),sm(:),tau(:);real(dp)::dd,hh,nu;integer::rot,j
    character(len=16)::kind
    allocate(a(0),b(0),sa(0),sm(0));if(present(ar))a=ar;if(present(ma))b=ma;if(present(sar))sa=sar;if(present(sma))sm=sma
    kind='arma';if(present(kpacf_model))kind=lower(kpacf_model);dd=0.2_dp;if(present(d))dd=d;hh=0.7_dp;if(present(hurst))hh=hurst
    select case(trim(kind))
    case('sarma4');tau=kpacf_sarma4(a,b,sa,sm,maxlag)
    case('sarma12');tau=kpacf_sarma12(a,b,sa,sm,maxlag)
    case('arfima');tau=kpacf_arfima(a,b,dd,maxlag)
    case('fbn');tau=kpacf_fbn(hh,maxlag)
    case default;tau=kpacf_arma(a,b,maxlag)
    end select
    model%kpacf_model=kind;allocate(model%vine%pairs(maxlag));rot=0;if(present(rotation))rot=rotation;nu=6.0_dp;if(present(df))nu=df
    do j=1,maxlag
      model%vine%pairs(j)=paircop(family,kendall_to_parameter(family_from_name(family),abs(tau(j))),nu,merge(rot,90, tau(j)>=0.0_dp))
      if(tau(j)<0.0_dp)model%vine%pairs(j)%rotation=merge(90,270,rot==0)
    end do
  end function dvinecopula2

  function dvinecopula3(base,locations,replacements) result(model)
    type(dvine2_copula),intent(in)::base;integer,intent(in)::locations(:)
    type(pair_copula),intent(in)::replacements(:);type(dvine3_copula)::model;integer::j
    allocate(model%vine%pairs(size(base%vine%pairs)));model%vine%pairs=base%vine%pairs
    do j=1,min(size(locations),size(replacements))
      if(locations(j)>=1.and.locations(j)<=size(model%vine%pairs))model%vine%pairs(locations(j))=replacements(j)
    end do
  end function dvinecopula3

  pure function lower(s) result(out)
    character(len=*),intent(in)::s;character(len=len(s))::out;integer::i,c
    do i=1,len(s);c=iachar(s(i:i));if(c>=65.and.c<=90)then;out(i:i)=achar(c+32);else;out(i:i)=s(i:i);end if;end do
  end function lower

  function mklist_dvine(model) result(pairs)
    type(dvine_copula),intent(in)::model;type(pair_copula),allocatable::pairs(:)
    allocate(pairs(size(model%pairs)));pairs=model%pairs
  end function mklist_dvine
  function mklist_dvine2(model) result(pairs)
    type(dvine2_copula),intent(in)::model;type(pair_copula),allocatable::pairs(:)
    allocate(pairs(size(model%vine%pairs)));pairs=model%vine%pairs
  end function mklist_dvine2
  function mklist_dvine3(model) result(pairs)
    type(dvine3_copula),intent(in)::model;type(pair_copula),allocatable::pairs(:)
    allocate(pairs(size(model%vine%pairs)));pairs=model%vine%pairs
  end function mklist_dvine3

  real(dp) function dvine_loglik(model,u) result(value)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::u(:)
    real(dp),allocatable::left(:),right(:),nl(:),nr(:);integer::j,n,m,i
    value=0.0_dp;n=size(u);if(n<2.or.size(model%pairs)==0)return
    allocate(left(n-1),right(n-1));left=u(1:n-1);right=u(2:n)
    do j=1,min(size(model%pairs),n-1)
      m=size(left);do i=1,m;value=value+pair_log_density(model%pairs(j),left(i),right(i));end do
      if(j<min(size(model%pairs),n-1).and.m>1)then
        allocate(nl(m-1),nr(m-1))
        do i=1,m-1;nl(i)=pair_h2(model%pairs(j),left(i),right(i));nr(i)=pair_h1(model%pairs(j),left(i+1),right(i+1));end do
        call move_alloc(nl,left);call move_alloc(nr,right)
      end if
    end do
  end function dvine_loglik

  real(dp) function dvine_objective(model,u) result(value)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::u(:);value=-dvine_loglik(model,u)
  end function dvine_objective

  real(dp) function rblattdens(model,history,x) result(value)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::history(:),x
    real(dp),allocatable::data(:),left(:),right(:),nl(:),nr(:);integer::k,j,m,i,n
    k=min(size(model%pairs),size(history));if(k==0)then;value=1.0_dp;return;end if
    allocate(data(k+1));n=size(history);data(1:k)=history(n-k+1:n);data(k+1)=x
    allocate(left(k),right(k));left=data(1:k);right=data(2:k+1);value=1.0_dp
    do j=1,k
      m=size(left);value=value*pair_density(model%pairs(j),left(m),right(m))
      if(j<k.and.m>1)then
        allocate(nl(m-1),nr(m-1));do i=1,m-1;nl(i)=pair_h2(model%pairs(j),left(i),right(i));nr(i)=pair_h1(model%pairs(j),left(i+1),right(i+1));end do
        call move_alloc(nl,left);call move_alloc(nr,right)
      end if
    end do
  end function rblattdens

  real(dp) function rblatt(model,history,x) result(value)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::history(:),x;type(cond_context)::ctx
    if(x<=0.0_dp)then;value=0.0_dp;return;else if(x>=1.0_dp)then;value=1.0_dp;return;end if
    ctx%model=model;allocate(ctx%history(size(history)));ctx%history=history
    value=integrate_simpson(cond_density_integrand,0.0_dp,x,ctx,2.0e-7_dp);value=min(max(value,0.0_dp),1.0_dp)
  end function rblatt

  real(dp) function cond_density_integrand(x,context_any) result(value)
    real(dp),intent(in)::x;class(*),intent(inout)::context_any
    select type(ctx=>context_any);type is(cond_context);value=rblattdens(ctx%model,ctx%history,x);class default;value=1.0_dp;end select
  end function cond_density_integrand

  real(dp) function irblatt(model,history,p) result(x)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::history(:),p;real(dp)::lo,hi,mid;integer::iter
    if(p<=0.0_dp)then;x=0.0_dp;return;else if(p>=1.0_dp)then;x=1.0_dp;return;end if
    lo=0.0_dp;hi=1.0_dp;do iter=1,60;mid=0.5_dp*(lo+hi);if(rblatt(model,history,mid)<p)then;lo=mid;else;hi=mid;end if;end do;x=0.5_dp*(lo+hi)
  end function irblatt

  function simdvine(pairs,n,innov,start) result(u)
    type(pair_copula),intent(in)::pairs(:);integer,intent(in)::n;real(dp),intent(in),optional::innov(:),start(:)
    real(dp),allocatable::u(:),e(:);type(dvine_copula)::model;integer::t,k,ns
    model=dvinecopula(pairs);allocate(u(n),e(n));if(present(innov))then;e=innov;else;do t=1,n;e(t)=uniform_random();end do;end if
    k=size(pairs);ns=0;if(present(start))ns=min(size(start),min(k,n));if(ns>0)u(1:ns)=start(1:ns)
    do t=ns+1,n
      if(t==1)then;u(t)=e(t);else;u(t)=irblatt(model,u(1:t-1),e(t));end if
    end do
  end function simdvine

  function sim_dvine(model,n,innov,start) result(u)
    type(dvine_copula),intent(in)::model;integer,intent(in)::n;real(dp),intent(in),optional::innov(:),start(:);real(dp),allocatable::u(:)
    if(present(innov).and.present(start))then;u=simdvine(model%pairs,n,innov,start)
    else if(present(innov))then;u=simdvine(model%pairs,n,innov=innov)
    else if(present(start))then;u=simdvine(model%pairs,n,start=start)
    else;u=simdvine(model%pairs,n);end if
  end function sim_dvine

  function resid_dvine(model,u) result(residuals)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::u(:);real(dp),allocatable::residuals(:);integer::t
    allocate(residuals(size(u)));if(size(u)==0)return;residuals(1)=u(1)
    do t=2,size(u);residuals(t)=rblatt(model,u(1:t-1),u(t));end do
  end function resid_dvine

  real(dp) function predict_dvine_cdf(model,history,x) result(value)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::history(:),x;value=rblatt(model,history,x)
  end function predict_dvine_cdf
  real(dp) function predict_dvine_quantile(model,history,p) result(value)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::history(:),p;value=irblatt(model,history,p)
  end function predict_dvine_quantile
  real(dp) function predict_dvine_density(model,history,x) result(value)
    type(dvine_copula),intent(in)::model;real(dp),intent(in)::history(:),x;value=rblattdens(model,history,x)
  end function predict_dvine_density

  function kendall_dvine(model) result(tau)
    type(dvine_copula),intent(in)::model;real(dp),allocatable::tau(:);integer::j
    allocate(tau(size(model%pairs)));do j=1,size(tau);tau(j)=pair_kendall(model%pairs(j));end do
  end function kendall_dvine

  function fit_dvine(u,families,rotations,df) result(fit)
    real(dp),intent(in)::u(:);character(len=*),intent(in)::families(:);integer,intent(in),optional::rotations(:);real(dp),intent(in),optional::df
    type(dvine_fit_result)::fit;real(dp),allocatable::left(:),right(:),nl(:),nr(:);real(dp)::tau,nu;integer::j,i,m,rot
    allocate(fit%model%pairs(size(families)),fit%kendall_parameters(size(families)));left=u(1:size(u)-1);right=u(2:size(u));nu=6.0_dp;if(present(df))nu=df
    do j=1,size(families)
      tau=sample_kendall(left,right);rot=0;if(present(rotations))rot=rotations(j)
      if(tau<0.0_dp.and.(rot==0.or.rot==180))rot=90
      fit%model%pairs(j)=paircop(families(j),kendall_to_parameter(family_from_name(families(j)),abs(tau)),nu,rot)
      fit%kendall_parameters(j)=tau;m=size(left)
      if(j<size(families).and.m>1)then
        allocate(nl(m-1),nr(m-1));do i=1,m-1;nl(i)=pair_h2(fit%model%pairs(j),left(i),right(i));nr(i)=pair_h1(fit%model%pairs(j),left(i+1),right(i+1));end do
        call move_alloc(nl,left);call move_alloc(nr,right)
      end if
    end do
    fit%objective=dvine_objective(fit%model,u);fit%convergence=0
  end function fit_dvine

  real(dp) function sample_kendall(x,y) result(tau)
    real(dp),intent(in)::x(:),y(:);integer::i,j,n;real(dp)::s
    n=min(size(x),size(y));s=0.0_dp
    do i=1,n-1
      do j=i+1,n
        if(abs((x(j)-x(i))*(y(j)-y(i)))>epsilon(1.0_dp))s=s+sign(1.0_dp,(x(j)-x(i))*(y(j)-y(i)))
      end do
    end do
    tau=2.0_dp*s/max(real(n*(n-1),dp),1.0_dp)
  end function sample_kendall

  function arma2dvine(family,ar,ma,maxlag,rotation,df) result(model)
    character(len=*),intent(in)::family;real(dp),intent(in)::ar(:),ma(:);integer,intent(in)::maxlag;integer,intent(in),optional::rotation;real(dp),intent(in),optional::df
    type(dvine2_copula)::model
    if(present(rotation).and.present(df))then;model=dvinecopula2(family,maxlag,'arma',ar=ar,ma=ma,rotation=rotation,df=df)
    else if(present(rotation))then;model=dvinecopula2(family,maxlag,'arma',ar=ar,ma=ma,rotation=rotation)
    else if(present(df))then;model=dvinecopula2(family,maxlag,'arma',ar=ar,ma=ma,df=df)
    else;model=dvinecopula2(family,maxlag,'arma',ar=ar,ma=ma);end if
  end function arma2dvine

  function sarma2dvine(family,ar,ma,sar,sma,period,maxlag,rotation,df) result(model)
    character(len=*),intent(in)::family;real(dp),intent(in)::ar(:),ma(:),sar(:),sma(:);integer,intent(in)::period,maxlag;integer,intent(in),optional::rotation;real(dp),intent(in),optional::df
    type(dvine2_copula)::model;character(len=8)::kind
    if(period==4)then;kind='sarma4';else;kind='sarma12';end if
    if(present(rotation).and.present(df))then;model=dvinecopula2(family,maxlag,kind,ar,ma,sar,sma,rotation=rotation,df=df)
    else if(present(rotation))then;model=dvinecopula2(family,maxlag,kind,ar,ma,sar,sma,rotation=rotation)
    else if(present(df))then;model=dvinecopula2(family,maxlag,kind,ar,ma,sar,sma,df=df)
    else;model=dvinecopula2(family,maxlag,kind,ar,ma,sar,sma);end if
  end function sarma2dvine

end module tscopula_dvine
