! SPDX-License-Identifier: GPL-3.0-only
module tscopula_compat
  use tscopula_kinds, only : dp
  use tscopula_status, only : tsc_error
  use tscopula_margins, only : margin_spec, margin_fit_result, fit_margin
  use tscopula_vtransforms, only : vtransform_spec
  use tscopula_paircopula, only : pair_copula
  use tscopula_timeseries, only : arma_copula, sarma_copula, arma_fit_result, &
    sim_arma_copula, sim_sarma_copula, fit_arma_copula, kendall_arma
  use tscopula_dvine, only : dvine_copula, dvine2_copula, dvine3_copula, dvine_fit_result, &
    sim_dvine, fit_dvine, kendall_dvine, arma2dvine
  use tscopula_models, only : tscopula_spec, vtscopula_spec, tscm_spec, tscm_fit_result, &
    sim_tscopula, sim_vtscopula, sim_tscm, kendall_tscopula, fit_tscm_steps, &
    tscopula_from_arma, vtscopula, predict_vtscopula_cdf, predict_vtscopula_quantile, &
    predict_vtscopula_density, copula_arma, copula_dvine_kind
  implicit none
  private

  interface sim
    module procedure sim_arma_wrap, sim_sarma_wrap, sim_dvine_wrap
    module procedure sim_dvine2_wrap, sim_dvine3_wrap, sim_tscopula_wrap
    module procedure sim_vtscopula_wrap, sim_tscm_wrap
  end interface sim
  interface kendall
    module procedure kendall_arma_wrap, kendall_dvine_wrap, kendall_dvine2_wrap
    module procedure kendall_dvine3_wrap, kendall_tscopula_wrap
  end interface kendall
  interface fit
    module procedure fit_margin_wrap, fit_arma_wrap, fit_dvine_wrap
    module procedure fit_vtscopula_wrap, fit_tscm_wrap
  end interface fit

  public :: sim, kendall, fit, armafit2dvine, vtparlist
  public :: dcondvtarma, pcondvtarma, qcondvtarma

contains

  function sim_arma_wrap(model,n) result(x)
    type(arma_copula),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:);x=sim_arma_copula(model,n)
  end function sim_arma_wrap
  function sim_sarma_wrap(model,n) result(x)
    type(sarma_copula),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:);x=sim_sarma_copula(model,n)
  end function sim_sarma_wrap
  function sim_dvine_wrap(model,n) result(x)
    type(dvine_copula),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:);x=sim_dvine(model,n)
  end function sim_dvine_wrap
  function sim_dvine2_wrap(model,n) result(x)
    type(dvine2_copula),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:);x=sim_dvine(model%vine,n)
  end function sim_dvine2_wrap
  function sim_dvine3_wrap(model,n) result(x)
    type(dvine3_copula),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:);x=sim_dvine(model%vine,n)
  end function sim_dvine3_wrap
  function sim_tscopula_wrap(model,n) result(x)
    type(tscopula_spec),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:);x=sim_tscopula(model,n)
  end function sim_tscopula_wrap
  function sim_vtscopula_wrap(model,n) result(x)
    type(vtscopula_spec),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:);x=sim_vtscopula(model,n)
  end function sim_vtscopula_wrap
  function sim_tscm_wrap(model,n) result(x)
    type(tscm_spec),intent(in)::model;integer,intent(in)::n;real(dp),allocatable::x(:);x=sim_tscm(model,n)
  end function sim_tscm_wrap

  function kendall_arma_wrap(model,maxlag) result(tau)
    type(arma_copula),intent(in)::model;integer,intent(in)::maxlag;real(dp),allocatable::tau(:);tau=kendall_arma(model,maxlag)
  end function kendall_arma_wrap
  function kendall_dvine_wrap(model,maxlag) result(tau)
    type(dvine_copula),intent(in)::model;integer,intent(in),optional::maxlag;real(dp),allocatable::tau(:),alltau(:);integer::m
    alltau=kendall_dvine(model);m=size(alltau);if(present(maxlag))m=min(m,maxlag);allocate(tau(m));tau=alltau(1:m)
  end function kendall_dvine_wrap
  function kendall_dvine2_wrap(model,maxlag) result(tau)
    type(dvine2_copula),intent(in)::model;integer,intent(in),optional::maxlag;real(dp),allocatable::tau(:)
    if(present(maxlag))then;tau=kendall_dvine_wrap(model%vine,maxlag);else;tau=kendall_dvine_wrap(model%vine);end if
  end function kendall_dvine2_wrap
  function kendall_dvine3_wrap(model,maxlag) result(tau)
    type(dvine3_copula),intent(in)::model;integer,intent(in),optional::maxlag;real(dp),allocatable::tau(:)
    if(present(maxlag))then;tau=kendall_dvine_wrap(model%vine,maxlag);else;tau=kendall_dvine_wrap(model%vine);end if
  end function kendall_dvine3_wrap
  function kendall_tscopula_wrap(model,maxlag) result(tau)
    type(tscopula_spec),intent(in)::model;integer,intent(in)::maxlag;real(dp),allocatable::tau(:);tau=kendall_tscopula(model,maxlag)
  end function kendall_tscopula_wrap

  function fit_margin_wrap(template,data,max_iter) result(result)
    type(margin_spec),intent(in)::template;real(dp),intent(in)::data(:);integer,intent(in),optional::max_iter
    type(margin_fit_result)::result;type(tsc_error)::error
    if(present(max_iter))then;call fit_margin(template,data,result,error,max_iter=max_iter);else;call fit_margin(template,data,result,error);end if
    if(.not.error%ok())result%convergence=2
  end function fit_margin_wrap
  function fit_arma_wrap(template,data,max_iter) result(result)
    type(arma_copula),intent(in)::template;real(dp),intent(in)::data(:);integer,intent(in),optional::max_iter;type(arma_fit_result)::result
    if(present(max_iter))then;result=fit_arma_copula(data,size(template%ar),size(template%ma),max_iter=max_iter)
    else;result=fit_arma_copula(data,size(template%ar),size(template%ma));end if
  end function fit_arma_wrap
  function fit_dvine_wrap(template,data) result(result)
    type(dvine_copula),intent(in)::template;real(dp),intent(in)::data(:);type(dvine_fit_result)::result
    character(len=12),allocatable::families(:);integer,allocatable::rot(:);integer::j
    allocate(families(size(template%pairs)),rot(size(template%pairs)))
    do j=1,size(template%pairs);families(j)=family_label(template%pairs(j));rot(j)=template%pairs(j)%rotation;end do
    result=fit_dvine(data,families,rot)
  end function fit_dvine_wrap
  function fit_vtscopula_wrap(template,data,max_iter) result(result)
    type(vtscopula_spec),intent(in)::template;real(dp),intent(in)::data(:);integer,intent(in),optional::max_iter;type(vtscopula_spec)::result
    type(arma_fit_result)::af;type(dvine_fit_result)::df;real(dp),allocatable::v(:);character(len=12),allocatable::families(:);integer::j
    result=template;allocate(v(size(data)));v=tscopula_vvalues(template,data)
    select case(template%base%kind)
    case(copula_arma)
      if(present(max_iter))then;af=fit_arma_copula(v,size(template%base%arma%ar),size(template%base%arma%ma),max_iter=max_iter)
      else;af=fit_arma_copula(v,size(template%base%arma%ar),size(template%base%arma%ma));end if
      result%base%arma=af%model
    case(copula_dvine_kind)
      allocate(families(size(template%base%dvine%pairs)));do j=1,size(families);families(j)=family_label(template%base%dvine%pairs(j));end do
      df=fit_dvine(v,families);result%base%dvine=df%model
    end select
  end function fit_vtscopula_wrap
  function fit_tscm_wrap(template,data,max_iter) result(result)
    type(tscm_spec),intent(in)::template;real(dp),intent(in)::data(:);integer,intent(in),optional::max_iter;type(tscm_fit_result)::result
    if(present(max_iter))then;result=fit_tscm_steps(template,data,max_iter);else;result=fit_tscm_steps(template,data);end if
  end function fit_tscm_wrap

  function armafit2dvine(family,fit_result,maxlag,rotation,df) result(model)
    character(len=*),intent(in)::family;type(arma_fit_result),intent(in)::fit_result;integer,intent(in)::maxlag
    integer,intent(in),optional::rotation;real(dp),intent(in),optional::df;type(dvine2_copula)::model
    if(present(rotation).and.present(df))then;model=arma2dvine(family,fit_result%model%ar,fit_result%model%ma,maxlag,rotation,df)
    else if(present(rotation))then;model=arma2dvine(family,fit_result%model%ar,fit_result%model%ma,maxlag,rotation)
    else if(present(df))then;model=arma2dvine(family,fit_result%model%ar,fit_result%model%ma,maxlag,df=df)
    else;model=arma2dvine(family,fit_result%model%ar,fit_result%model%ma,maxlag);end if
  end function armafit2dvine

  function vtparlist(model) result(par)
    type(vtscopula_spec),intent(in)::model;real(dp),allocatable::par(:);integer::n
    n=3+merge(2,0,model%has_wcopula);allocate(par(n));par(1:3)=[model%transform%delta,model%transform%kappa,model%transform%xi]
    if(model%has_wcopula)par(4:5)=[model%wcopula%par1,model%wcopula%par2]
  end function vtparlist

  real(dp) function dcondvtarma(model,transform,history,x) result(value)
    type(arma_copula),intent(in)::model;type(vtransform_spec),intent(in)::transform;real(dp),intent(in)::history(:),x
    value=predict_vtscopula_density(vtscopula(tscopula_from_arma(model),transform),history,x)
  end function dcondvtarma
  real(dp) function pcondvtarma(model,transform,history,x) result(value)
    type(arma_copula),intent(in)::model;type(vtransform_spec),intent(in)::transform;real(dp),intent(in)::history(:),x
    value=predict_vtscopula_cdf(vtscopula(tscopula_from_arma(model),transform),history,x)
  end function pcondvtarma
  real(dp) function qcondvtarma(model,transform,history,p) result(value)
    type(arma_copula),intent(in)::model;type(vtransform_spec),intent(in)::transform;real(dp),intent(in)::history(:),p
    value=predict_vtscopula_quantile(vtscopula(tscopula_from_arma(model),transform),history,p)
  end function qcondvtarma

  function tscopula_vvalues(model,u) result(v)
    use tscopula_vtransforms, only : vtrans
    type(vtscopula_spec),intent(in)::model;real(dp),intent(in)::u(:);real(dp),allocatable::v(:);allocate(v(size(u)));v=vtrans(model%transform,u)
  end function tscopula_vvalues

  function family_label(cop) result(label)
    use tscopula_paircopula, only : family_name
    type(pair_copula),intent(in)::cop;character(len=12)::label;character(len=:),allocatable::tmp
    tmp=family_name(cop%family);label=' ';label(1:min(len_trim(tmp),len(label)))=tmp(1:min(len_trim(tmp),len(label)))
  end function family_label

end module tscopula_compat
