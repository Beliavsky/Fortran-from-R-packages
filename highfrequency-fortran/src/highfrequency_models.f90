! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_models
  use highfrequency_kinds, only: dp, pi
  use highfrequency_stats, only: mean_value
  use highfrequency_linalg, only: least_squares
  use highfrequency_optimize, only: nelder_mead_bounded
  use highfrequency_types, only: har_model, heavy_model
  implicit none
  private
  public :: fit_har, har_forecast, fit_heavy, heavy_forecast
  public :: heavy_recursion

contains

  function fit_har(series, periods, log_transform, external) result(model)
    real(dp), intent(in) :: series(:)
    integer, intent(in) :: periods(:)
    logical, intent(in), optional :: log_transform
    real(dp), intent(in), optional :: external(:,:)
    type(har_model) :: model
    real(dp), allocatable :: x(:,:), y(:), transformed(:), beta(:)
    logical :: use_log, ok
    integer :: n, p, max_period, nobs, i, j, t, nexternal

    n=size(series)
    p=size(periods)
    use_log=.false.
    if(present(log_transform)) use_log=log_transform
    nexternal=0
    if(present(external)) then
      if(size(external,1)==n) nexternal=size(external,2)
    end if
    max_period=maxval(periods)
    nobs=n-max_period
    if(nobs<=p+nexternal .or. minval(periods)<1) return
    allocate(transformed(n))
    if(use_log)then
      transformed=log(max(series,tiny(1.0_dp)))
    else
      transformed=series
    end if
    allocate(x(nobs,1+p+nexternal),y(nobs),beta(1+p+nexternal))
    x(:,1)=1.0_dp
    do i=1,nobs
      t=max_period+i-1
      y(i)=transformed(t+1)
      do j=1,p
        x(i,1+j)=sum(transformed(t-periods(j)+1:t))/real(periods(j),dp)
      end do
      if(nexternal>0) x(i,2+p:1+p+nexternal)=external(t,1:nexternal)
    end do
    call least_squares(x,y,beta,ok)
    if(.not.ok) return
    model%nobs=nobs
    model%nreg=size(beta)
    allocate(model%coefficients(size(beta)),model%fitted(nobs),model%residuals(nobs))
    model%coefficients=beta
    model%fitted=matmul(x,beta)
    if(use_log) model%fitted=exp(model%fitted)
    if(use_log)then
      model%residuals=series(max_period+1:n)-model%fitted
    else
      model%residuals=y-model%fitted
    end if
    model%sigma2=sum(model%residuals**2)/real(max(1,nobs-size(beta)),dp)
    model%fitted_ok=.true.
  end function fit_har

  real(dp) function har_forecast(model, history, periods, log_transform, external_values) result(value)
    type(har_model), intent(in) :: model
    real(dp), intent(in) :: history(:)
    integer, intent(in) :: periods(:)
    logical, intent(in), optional :: log_transform
    real(dp), intent(in), optional :: external_values(:)
    real(dp), allocatable :: regressors(:), transformed(:)
    logical :: use_log
    integer :: j,p,n,ne
    value=0.0_dp
    if(.not.model%fitted_ok) return
    n=size(history)
    p=size(periods)
    if(n<maxval(periods)) return
    ne=0
    if(present(external_values)) ne=size(external_values)
    if(size(model%coefficients)/=1+p+ne) return
    use_log=.false.;if(present(log_transform))use_log=log_transform
    allocate(transformed(n),regressors(p+ne))
    if(use_log)then
      transformed=log(max(history,tiny(1.0_dp)))
    else
      transformed=history
    end if
    do j=1,p
      regressors(j)=sum(transformed(n-periods(j)+1:n))/real(periods(j),dp)
    end do
    if(ne>0)regressors(p+1:p+ne)=external_values
    value=model%predict_one(regressors)
    if(use_log)value=exp(value)
  end function har_forecast

  subroutine heavy_recursion(parameters, input, conditional)
    real(dp),intent(in)::parameters(3),input(:)
    real(dp),intent(out)::conditional(size(input))
    integer::t
    if(size(input)==0)return
    conditional(1)=max(tiny(1.0_dp),mean_value(input))
    do t=2,size(input)
      conditional(t)=parameters(1)+parameters(2)*input(t-1)+parameters(3)*conditional(t-1)
      conditional(t)=max(tiny(1.0_dp),conditional(t))
    end do
  end subroutine heavy_recursion

  function fit_heavy(returns, realized_measure, starting_values) result(model)
    real(dp),intent(in)::returns(:),realized_measure(:)
    real(dp),intent(in),optional::starting_values(6)
    type(heavy_model)::model
    real(dp),allocatable::ret(:),ret2(:),h(:),m(:)
    real(dp)::start1(3),start2(3),lower(3),upper(3),sol1(3),sol2(3),v1,v2
    logical::ok1,ok2
    integer::n
    n=min(size(returns),size(realized_measure))
    if(n<10)return
    allocate(ret(n),ret2(n),h(n),m(n))
    ret=returns(:n)-mean_value(returns(:n))
    ret2=ret*ret
    if(present(starting_values))then
      start1=starting_values(1:3)
      start2=starting_values(4:6)
    else
      start1=[max(tiny(1.0_dp),mean_value(ret2)*0.2_dp),0.3_dp,0.5_dp]
      start2=[max(tiny(1.0_dp),mean_value(realized_measure(:n))*0.1_dp),0.6_dp,0.3_dp]
    end if
    lower=[tiny(1.0_dp),0.0_dp,0.0_dp]
    upper=[max(1.0_dp,10.0_dp*mean_value(ret2)),0.999_dp,0.999_dp]
    call nelder_mead_bounded(objective_returns,start1,lower,upper,sol1,v1,ok1,1200,1.0e-8_dp)
    upper(1)=max(1.0_dp,10.0_dp*mean_value(realized_measure(:n)))
    call nelder_mead_bounded(objective_measure,start2,lower,upper,sol2,v2,ok2,1200,1.0e-8_dp)
    call heavy_recursion(sol1,realized_measure(:n),h)
    call heavy_recursion(sol2,realized_measure(:n),m)
    model%omega=sol1(1);model%alpha=sol1(2);model%beta=sol1(3)
    model%omega_rm=sol2(1);model%alpha_rm=sol2(2);model%beta_rm=sol2(3)
    allocate(model%variance(n),model%measure_variance(n))
    model%variance=h;model%measure_variance=m
    model%loglik_returns=-v1
    model%loglik_measure=-v2
    model%fitted_ok=ok1.and.ok2
  contains
    function objective_returns(par) result(value)
      real(dp),intent(in)::par(:)
      real(dp)::value
      real(dp)::cond(n)
      if(par(2)+par(3)>=0.999999_dp)then
        value=huge(1.0_dp)
        return
      end if
      call heavy_recursion(par(1:3),realized_measure(:n),cond)
      value=0.5_dp*sum(log(2.0_dp*pi)+log(cond)+ret2/cond)
    end function objective_returns
    function objective_measure(par) result(value)
      real(dp),intent(in)::par(:)
      real(dp)::value
      real(dp)::cond(n)
      if(par(2)+par(3)>=0.999999_dp)then
        value=huge(1.0_dp)
        return
      end if
      call heavy_recursion(par(1:3),realized_measure(:n),cond)
      value=0.5_dp*sum(log(2.0_dp*pi)+log(cond)+realized_measure(:n)/cond)
    end function objective_measure
  end function fit_heavy

  function heavy_forecast(model,last_measure,steps) result(forecast)
    type(heavy_model),intent(in)::model
    real(dp),intent(in)::last_measure
    integer,intent(in)::steps
    real(dp)::forecast(max(0,steps),2)
    real(dp)::h,m,nu
    integer::i
    if(steps<=0)return
    if(.not.model%fitted_ok)then
      forecast=0.0_dp
      return
    end if
    h=model%omega+model%alpha*last_measure+model%beta*model%variance(size(model%variance))
    m=model%omega_rm+model%alpha_rm*last_measure+model%beta_rm* &
      model%measure_variance(size(model%measure_variance))
    forecast(1,:)=[h,m]
    nu=model%alpha_rm+model%beta_rm
    do i=2,steps
      h=model%omega+model%alpha*m+model%beta*h
      m=model%omega_rm+nu*m
      forecast(i,:)=[h,m]
    end do
  end function heavy_forecast

end module highfrequency_models
