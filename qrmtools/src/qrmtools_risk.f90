! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_risk
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
  use qrmtools_kinds, only : dp
  use qrmtools_types, only : fit_result
  use qrmtools_stats, only : quantile_type1, normal_quantile, normal_pdf, student_t_quantile, &
    student_t_pdf, mean_value
  use qrmtools_distributions, only : qgpd,qpar,qgpdtail
  use qrmtools_optimization, only : nelder_mead
  implicit none
  private
  public :: var_np,var_t,var_t01,var_gpd,var_par,var_gpdtail
  public :: es_np,es_t,es_t01,es_gpd,es_par,es_gpdtail,rvar_np
  public :: geometric_objective,geometric_var,geometric_expectile
  real(dp), allocatable, save :: geom_data(:,:),geom_level(:)
  logical, save :: geom_expectile=.false.
contains
  real(dp) function var_np(x,level) result(value)
    real(dp), intent(in) :: x(:),level
    value=quantile_type1(x,level)
  end function var_np

  real(dp) function var_t(level,loc,scale,df) result(value)
    real(dp), intent(in) :: level
    real(dp), intent(in), optional :: loc,scale,df
    real(dp) :: mu,sigma,nu
    mu=0.0_dp; sigma=1.0_dp; nu=huge(1.0_dp)
    if(present(loc))mu=loc; if(present(scale))sigma=scale; if(present(df))nu=df
    if(nu>1.0e15_dp) then; value=mu+sigma*normal_quantile(level)
    else; value=mu+sigma*student_t_quantile(level,nu); end if
  end function var_t

  real(dp) function var_t01(level,df) result(value)
    real(dp), intent(in) :: level,df
    if(df<=2.0_dp) then; value=ieee_value(1.0_dp,ieee_quiet_nan)
    else if(df>1.0e15_dp) then; value=normal_quantile(level)
    else; value=sqrt((df-2.0_dp)/df)*student_t_quantile(level,df); end if
  end function var_t01

  pure real(dp) function var_gpd(level,shape,scale) result(value)
    real(dp), intent(in) :: level,shape,scale
    value=qgpd(level,shape,scale)
  end function var_gpd

  pure real(dp) function var_par(level,shape,scale) result(value)
    real(dp), intent(in) :: level,shape
    real(dp), intent(in), optional :: scale
    value=qpar(level,shape,scale)
  end function var_par

  pure real(dp) function var_gpdtail(level,threshold,p_exceed,shape,scale) result(value)
    real(dp), intent(in) :: level,threshold,p_exceed,shape,scale
    value=qgpdtail(level,threshold,p_exceed,shape,scale)
  end function var_gpdtail

  real(dp) function es_np(x,level,include_var) result(value)
    real(dp), intent(in) :: x(:),level
    logical, intent(in), optional :: include_var
    real(dp) :: v
    logical :: inclusive
    integer :: n
    inclusive=.false.; if(present(include_var))inclusive=include_var
    v=var_np(x,level)
    if(inclusive) then
      n=count(x>=v); if(n>0)value=sum(x,mask=x>=v)/real(n,dp)
    else
      n=count(x>v); if(n>0)value=sum(x,mask=x>v)/real(n,dp)
    end if
    if(n==0)value=ieee_value(1.0_dp,ieee_quiet_nan)
  end function es_np

  real(dp) function es_t(level,loc,scale,df) result(value)
    real(dp), intent(in) :: level
    real(dp), intent(in), optional :: loc,scale,df
    real(dp) :: mu,sigma,nu,q
    mu=0.0_dp; sigma=1.0_dp; nu=huge(1.0_dp)
    if(present(loc))mu=loc; if(present(scale))sigma=scale; if(present(df))nu=df
    if(level>=1.0_dp) then; value=ieee_value(1.0_dp,ieee_positive_inf); return; end if
    if(nu>1.0e15_dp) then
      q=normal_quantile(level); value=mu+sigma*normal_pdf(q)/(1.0_dp-level)
    else if(nu<=1.0_dp) then
      value=ieee_value(1.0_dp,ieee_positive_inf)
    else
      q=student_t_quantile(level,nu)
      value=mu+sigma*student_t_pdf(q,nu)*(nu+q*q)/((nu-1.0_dp)*(1.0_dp-level))
    end if
  end function es_t

  real(dp) function es_t01(level,df) result(value)
    real(dp), intent(in) :: level,df
    if(df<=2.0_dp) then; value=ieee_value(1.0_dp,ieee_quiet_nan)
    else if(df>1.0e15_dp) then; value=es_t(level)
    else; value=sqrt((df-2.0_dp)/df)*es_t(level,df=df); end if
  end function es_t01

  pure real(dp) function es_gpd(level,shape,scale) result(value)
    real(dp), intent(in) :: level,shape,scale
    real(dp) :: v
    if(shape>=1.0_dp) then; value=ieee_value(1.0_dp,ieee_positive_inf)
    else; v=qgpd(level,shape,scale); value=v+(scale+shape*v)/(1.0_dp-shape); end if
  end function es_gpd

  pure real(dp) function es_par(level,shape,scale) result(value)
    real(dp), intent(in) :: level,shape
    real(dp), intent(in), optional :: scale
    real(dp) :: k
    k=1.0_dp; if(present(scale))k=scale
    if(shape<=1.0_dp) then; value=ieee_value(1.0_dp,ieee_positive_inf)
    else; value=k*(shape/(shape-1.0_dp)*(1.0_dp-level)**(-1.0_dp/shape)-1.0_dp); end if
  end function es_par

  pure real(dp) function es_gpdtail(level,threshold,p_exceed,shape,scale) result(value)
    real(dp), intent(in) :: level,threshold,p_exceed,shape,scale
    real(dp) :: v
    if(shape>=1.0_dp) then; value=ieee_value(1.0_dp,ieee_positive_inf)
    else; v=qgpdtail(level,threshold,p_exceed,shape,scale); value=(v+scale-shape*threshold)/(1.0_dp-shape); end if
  end function es_gpdtail

  real(dp) function rvar_np(x,lower_level,upper_level) result(value)
    real(dp), intent(in) :: x(:),lower_level
    real(dp), intent(in), optional :: upper_level
    real(dp) :: upper,v1,v2
    integer :: n
    upper=1.0_dp; if(present(upper_level))upper=upper_level
    v1=var_np(x,lower_level); v2=var_np(x,upper)
    n=count(x>v1 .and. x<=v2)
    if(n>0) then; value=sum(x,mask=x>v1 .and. x<=v2)/real(n,dp)
    else; value=ieee_value(1.0_dp,ieee_quiet_nan); end if
  end function rvar_np

  pure real(dp) function geometric_objective(center,x,level,expectile) result(value)
    real(dp), intent(in) :: center(:),x(:,:),level(:)
    logical, intent(in), optional :: expectile
    real(dp) :: diff(size(x,2)),normv,s,term
    integer :: i
    logical :: exp_measure
    exp_measure=.false.; if(present(expectile))exp_measure=expectile
    value=0.0_dp
    do i=1,size(x,1)
      diff=x(i,:)-center; normv=sqrt(sum(diff*diff)); s=sum(diff*level)
      if(exp_measure) then; term=0.5_dp*normv*(normv+s)
      else; term=0.5_dp*(normv+s); end if
      value=value+term
    end do
    value=value/real(size(x,1),dp)
  end function geometric_objective

  function geometric_var(x,level,start,max_iterations) result(fit)
    real(dp), intent(in) :: x(:,:),level(:)
    real(dp), intent(in), optional :: start(:)
    integer, intent(in), optional :: max_iterations
    type(fit_result) :: fit
    real(dp), allocatable :: initial(:),optimum(:)
    real(dp) :: best
    integer :: iterations,evaluations,maxit
    logical :: converged
    if(size(x,2)/=size(level)) then; fit%message='Level dimension mismatch.'; return; end if
    allocate(initial(size(level))); initial=sum(x,dim=1)/real(size(x,1),dp); if(present(start))initial=start
    geom_data=x; geom_level=level; geom_expectile=.false.; maxit=2000; if(present(max_iterations))maxit=max_iterations
    call nelder_mead(geom_context_objective,initial,optimum,best,iterations,evaluations,converged,maxit,1.0e-9_dp)
    fit%parameters=optimum; fit%log_likelihood=-best; fit%iterations=iterations; fit%evaluations=evaluations
    fit%converged=converged; fit%ok=.true.
  end function geometric_var

  function geometric_expectile(x,level,start,max_iterations) result(fit)
    real(dp), intent(in) :: x(:,:),level(:)
    real(dp), intent(in), optional :: start(:)
    integer, intent(in), optional :: max_iterations
    type(fit_result) :: fit
    real(dp), allocatable :: initial(:),optimum(:)
    real(dp) :: best
    integer :: iterations,evaluations,maxit
    logical :: converged
    if(size(x,2)/=size(level)) then; fit%message='Level dimension mismatch.'; return; end if
    allocate(initial(size(level))); initial=sum(x,dim=1)/real(size(x,1),dp); if(present(start))initial=start
    geom_data=x; geom_level=level; geom_expectile=.true.; maxit=2000; if(present(max_iterations))maxit=max_iterations
    call nelder_mead(geom_context_objective,initial,optimum,best,iterations,evaluations,converged,maxit,1.0e-9_dp)
    fit%parameters=optimum; fit%log_likelihood=-best; fit%iterations=iterations; fit%evaluations=evaluations
    fit%converged=converged; fit%ok=.true.
  end function geometric_expectile

  real(dp) function geom_context_objective(center) result(value)
    real(dp), intent(in) :: center(:)
    value=geometric_objective(center,geom_data,geom_level,geom_expectile)
  end function geom_context_objective
end module qrmtools_risk
