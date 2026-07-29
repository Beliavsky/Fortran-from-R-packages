! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_evt
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use qrmtools_kinds, only : dp
  use qrmtools_types, only : fit_result, hill_result
  use qrmtools_stats, only : mean_value, variance_value, sort_increasing, quantile_type1, &
    normal_quantile, numerical_hessian, invert_matrix
  use qrmtools_distributions, only : dgev, dgpd, pgpd
  use qrmtools_optimization, only : nelder_mead
  implicit none
  private
  public :: fit_gev_quantile, fit_gev_pwm, loglik_gev, fit_gev_mle
  public :: fit_gpd_mom, fit_gpd_pwm, loglik_gpd, fit_gpd_mle
  public :: hill_estimator, mean_excess_np, mean_excess_gpd, tail_estimator_gpd

  real(dp), allocatable, save :: context_data(:)
  integer, save :: context_model=0
contains
  function fit_gev_quantile(x,p,cutoff) result(parameters)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: p(3),cutoff
    real(dp) :: parameters(3),prob(3),q(3),qd(2),y,l(3),a(2),ar,c,shape,scale,loc,m1,m2
    integer :: i
    prob=[0.25_dp,0.5_dp,0.75_dp]; if(present(p))prob=p
    c=3.0_dp; if(present(cutoff))c=cutoff
    do i=1,3; q(i)=quantile_type1(x,prob(i)); end do
    qd=[q(2)-q(1),q(3)-q(2)]
    if(minval(abs(qd))<=tiny(1.0_dp)) then
      parameters=ieee_value(1.0_dp,ieee_quiet_nan); return
    end if
    y=qd(1)/qd(2); l=log(-log(prob)); a=[l(1)-l(2),l(2)-l(3)]; ar=a(1)/a(2)
    if(y<1.0_dp/(exp(c/ar)-1.0_dp)) then
      shape=log(1.0_dp+1.0_dp/y)/a(2)
    else if(y<=ar) then
      m1=a(2)/c*log(ar/(exp(ar*c)-1.0_dp)); shape=log(y/ar)/m1
    else if(y<=(exp(ar*c)-1.0_dp)) then
      m2=-a(1)/c*log(ar*(exp(c/ar)-1.0_dp)); shape=log(y/ar)/m2
    else
      shape=-log(1.0_dp+y)/a(1)
    end if
    if(abs(shape)<=1.0e-10_dp) then
      scale=qd(1)/(-l(2)+l(1)); loc=q(2)+scale*l(2)
    else
      scale=shape*qd(1)/((-log(prob(2)))**(-shape)-(-log(prob(1)))**(-shape))
      loc=q(2)-scale/shape*((-log(prob(2)))**(-shape)-1.0_dp)
    end if
    parameters=[shape,loc,scale]
  end function fit_gev_quantile

  function fit_gev_pwm(x) result(parameters)
    real(dp), intent(in) :: x(:)
    real(dp) :: parameters(3),b0,b1,b2,y,a,b,c,shape,scale,loc
    real(dp), allocatable :: z(:)
    integer :: n,k
    n=size(x); z=x; call sort_increasing(z)
    b0=mean_value(x); b1=0.0_dp; b2=0.0_dp
    do k=1,n
      b1=b1+z(k)*real(k-1,dp)/real(n-1,dp)
      b2=b2+z(k)*real((k-1)*(k-2),dp)/real((n-1)*(n-2),dp)
    end do
    b1=b1/real(n,dp); b2=b2/real(n,dp)
    y=(3.0_dp*b2-b0)/(2.0_dp*b1-b0)
    if(y<=1.0_dp) then; parameters=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    a=-2.9554_dp; b=-4.1297_dp; c=3.782014_dp
    shape=a/y**2+b/y+c
    if(abs(shape)<=1.0e-10_dp) then
      scale=(2.0_dp*b1-b0)/log(2.0_dp); loc=b0-scale*0.5772156649015329_dp
    else
      scale=(2.0_dp*b1-b0)*shape/(gamma(1.0_dp-shape)*(2.0_dp**shape-1.0_dp))
      loc=b0-scale*(gamma(1.0_dp-shape)-1.0_dp)/shape
    end if
    parameters=[shape,loc,max(scale,tiny(1.0_dp))]
  end function fit_gev_pwm

  real(dp) function loglik_gev(parameters,x) result(value)
    real(dp), intent(in) :: parameters(:),x(:)
    integer :: i
    value=0.0_dp
    do i=1,size(x)
      value=value+dgev(x(i),parameters(1),parameters(2),parameters(3),.true.)
      if(.not.ieee_is_finite(value)) return
    end do
  end function loglik_gev

  function fit_gev_mle(x,initial,estimate_covariance,max_iterations) result(fit)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: initial(3)
    logical, intent(in), optional :: estimate_covariance
    integer, intent(in), optional :: max_iterations
    type(fit_result) :: fit
    real(dp) :: start(3),best
    real(dp), allocatable :: optimum(:),hessian(:,:),information(:,:),inverse(:,:)
    logical :: converged,estimate_cov,ok
    integer :: iterations,evaluations,maxit
    if(size(x)<3) then; fit%message='At least three observations are required.'; return; end if
    if(present(initial)) then; start=initial; else; start=fit_gev_pwm(x); end if
    if(.not.all(ieee_is_finite(start)) .or. start(3)<=0.0_dp) start=[0.0_dp,mean_value(x),sqrt(variance_value(x))]
    do while(.not.ieee_is_finite(loglik_gev(start,x))); start(3)=2.0_dp*start(3); end do
    context_data=x; context_model=1; maxit=2000; if(present(max_iterations))maxit=max_iterations
    call nelder_mead(negative_loglik,start,optimum,best,iterations,evaluations,converged,maxit,1.0e-9_dp)
    fit%parameters=optimum; fit%log_likelihood=-best; fit%iterations=iterations
    fit%evaluations=evaluations; fit%converged=converged; fit%ok=ieee_is_finite(best)
    if(.not.fit%ok)fit%message='Optimization failed.'
    estimate_cov=.true.; if(present(estimate_covariance))estimate_cov=estimate_covariance
    if(estimate_cov .and. fit%ok) then
      hessian=numerical_hessian(negative_loglik,optimum); information=hessian
      call invert_matrix(information,inverse,ok)
      if(ok) then
        fit%covariance=0.5_dp*(inverse+transpose(inverse)); allocate(fit%standard_errors(3))
        fit%standard_errors=sqrt(max([fit%covariance(1,1),fit%covariance(2,2),fit%covariance(3,3)],0.0_dp))
      end if
    end if
  end function fit_gev_mle

  function fit_gpd_mom(x) result(parameters)
    real(dp), intent(in) :: x(:)
    real(dp) :: parameters(2),m,v,shape
    m=mean_value(x); v=variance_value(x); shape=0.5_dp*(1.0_dp-m*m/v)
    parameters=[shape,max(m*(1.0_dp-shape),epsilon(1.0_dp))]
  end function fit_gpd_mom

  function fit_gpd_pwm(x) result(parameters)
    real(dp), intent(in) :: x(:)
    real(dp) :: parameters(2),a0,a1
    real(dp), allocatable :: z(:)
    integer :: n,k
    n=size(x); z=x; call sort_increasing(z); a0=mean_value(x); a1=0.0_dp
    do k=1,n; a1=a1+z(k)*real(n-k,dp)/real(n-1,dp); end do
    a1=a1/real(n,dp)
    parameters=[2.0_dp-a0/(a0-2.0_dp*a1),max(2.0_dp*a0*a1/(a0-2.0_dp*a1),epsilon(1.0_dp))]
  end function fit_gpd_pwm

  real(dp) function loglik_gpd(parameters,x) result(value)
    real(dp), intent(in) :: parameters(:),x(:)
    integer :: i
    value=0.0_dp
    do i=1,size(x)
      value=value+dgpd(x(i),parameters(1),parameters(2),.true.)
      if(.not.ieee_is_finite(value)) return
    end do
  end function loglik_gpd

  function fit_gpd_mle(x,initial,estimate_covariance,max_iterations) result(fit)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: initial(2)
    logical, intent(in), optional :: estimate_covariance
    integer, intent(in), optional :: max_iterations
    type(fit_result) :: fit
    real(dp) :: start(2),best,mx
    real(dp), allocatable :: optimum(:),hessian(:,:),inverse(:,:)
    logical :: converged,estimate_cov,ok
    integer :: iterations,evaluations,maxit
    if(any(x<0.0_dp) .or. size(x)<2) then; fit%message='GPD data must be nonnegative.'; return; end if
    if(present(initial)) then; start=initial; else; start=fit_gpd_pwm(x); end if
    if(.not.all(ieee_is_finite(start)) .or. start(2)<=0.0_dp) start=[0.0_dp,max(mean_value(x),epsilon(1.0_dp))]
    mx=maxval(x)
    if(start(1)<0.0_dp .and. mx>=-start(2)/start(1)) start(2)=-start(1)*mx*1.01_dp
    context_data=x; context_model=2; maxit=2000; if(present(max_iterations))maxit=max_iterations
    call nelder_mead(negative_loglik,start,optimum,best,iterations,evaluations,converged,maxit,1.0e-9_dp)
    fit%parameters=optimum; fit%log_likelihood=-best; fit%iterations=iterations
    fit%evaluations=evaluations; fit%converged=converged; fit%ok=ieee_is_finite(best)
    if(.not.fit%ok)fit%message='Optimization failed.'
    estimate_cov=.true.; if(present(estimate_covariance))estimate_cov=estimate_covariance
    if(estimate_cov .and. fit%ok) then
      hessian=numerical_hessian(negative_loglik,optimum); call invert_matrix(hessian,inverse,ok)
      if(ok) then
        fit%covariance=0.5_dp*(inverse+transpose(inverse)); allocate(fit%standard_errors(2))
        fit%standard_errors=sqrt(max([fit%covariance(1,1),fit%covariance(2,2)],0.0_dp))
      end if
    end if
  end function fit_gpd_mle

  real(dp) function negative_loglik(parameters) result(value)
    real(dp), intent(in) :: parameters(:)
    real(dp) :: ll
    if(context_model==1) then; ll=loglik_gev(parameters,context_data)
    else; ll=loglik_gpd(parameters,context_data); end if
    if(ieee_is_finite(ll)) then; value=-ll; else; value=1.0e150_dp; end if
  end function negative_loglik

  function hill_estimator(x,k_min,k_max,confidence_level) result(output)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: k_min,k_max
    real(dp), intent(in), optional :: confidence_level
    type(hill_result) :: output
    real(dp), allocatable :: z(:),logz(:),cumulative(:)
    real(dp) :: level,q,se
    integer :: n,klo,khi,k,i
    n=size(x); klo=10; if(present(k_min))klo=k_min; khi=n; if(present(k_max))khi=k_max
    level=0.95_dp; if(present(confidence_level))level=confidence_level
    if(klo<2 .or. khi>n .or. klo>khi) then; output%message='Invalid Hill order-statistic range.'; return; end if
    z=x; call sort_increasing(z); z=z(n:1:-1)
    if(z(khi)<=0.0_dp) then; output%message='Hill observations must be positive.'; return; end if
    logz=log(z(1:khi)); allocate(cumulative(khi)); cumulative(1)=logz(1)
    do i=2,khi; cumulative(i)=cumulative(i-1)+logz(i); end do
    allocate(output%k(khi-klo+1),output%probability(khi-klo+1),output%tail_index(khi-klo+1),&
      output%ci_low(khi-klo+1),output%ci_high(khi-klo+1))
    q=normal_quantile(1.0_dp-(1.0_dp-level)/2.0_dp)
    do k=klo,khi
      i=k-klo+1; output%k(i)=k; output%probability(i)=1.0_dp-real(k-1,dp)/real(n,dp)
      output%tail_index(i)=1.0_dp/(cumulative(k)/real(k,dp)-logz(k))
      se=output%tail_index(i)/sqrt(real(k,dp)); output%ci_low(i)=output%tail_index(i)-q*se
      output%ci_high(i)=output%tail_index(i)+q*se
    end do
    output%ok=.true.
  end function hill_estimator

  function mean_excess_np(x,omit) result(values)
    real(dp), intent(in) :: x(:)
    integer, intent(in), optional :: omit
    real(dp), allocatable :: values(:,:)
    real(dp), allocatable :: z(:),unique_values(:),mean_excess(:)
    integer, allocatable :: counts(:)
    integer :: n,o,nu,i,j,k,total
    n=size(x); o=3; if(present(omit))o=omit; z=x; call sort_increasing(z)
    allocate(unique_values(n),counts(n)); nu=0
    nu=1
    unique_values(1)=z(1)
    counts(1)=1
    do i=2,n
      if(z(i)>z(i-1)) then
        nu=nu+1
        unique_values(nu)=z(i)
        counts(nu)=1
      else
        counts(nu)=counts(nu)+1
      end if
    end do
    if(nu-o<1) then; allocate(values(0,2)); return; end if
    allocate(mean_excess(nu-o));
    do k=1,nu-o; mean_excess(k)=sum(unique_values(k+1:nu))/real(nu-k,dp)-unique_values(k); end do
    total=sum(counts(1:nu-o)); allocate(values(total,2)); j=0
    do k=1,nu-o; do i=1,counts(k); j=j+1; values(j,:)=[unique_values(k),mean_excess(k)]; end do; end do
  end function mean_excess_np

  pure real(dp) function mean_excess_gpd(x,shape,scale) result(value)
    real(dp), intent(in) :: x,shape,scale
    if(shape>=1.0_dp .or. scale+shape*x<=0.0_dp) then
      value=ieee_value(1.0_dp,ieee_quiet_nan)
    else
      value=(scale+shape*x)/(1.0_dp-shape)
    end if
  end function mean_excess_gpd

  pure real(dp) function tail_estimator_gpd(q,threshold,p_exceed,shape,scale) result(value)
    real(dp), intent(in) :: q,threshold,p_exceed,shape,scale
    if(q<threshold) then; value=ieee_value(1.0_dp,ieee_quiet_nan)
    else; value=p_exceed*pgpd(q-threshold,shape,scale,.false.); end if
  end function tail_estimator_gpd
end module qrmtools_evt
