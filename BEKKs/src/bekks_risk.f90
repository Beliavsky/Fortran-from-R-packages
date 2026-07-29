! SPDX-License-Identifier: MIT
module bekks_risk
  use bekks_kinds, only: dp
  use bekks_types
  use bekks_math, only: normal_quantile, student_t_quantile, empirical_quantile, &
    sample_kurtosis, chi_square_cdf
  use bekks_forecast, only: forecast_bekk
  use bekks_estimation, only: bekk_fit
  implicit none
  private
  public :: var_bekk_fit, var_bekk_forecast
  public :: coverage_tests, backtest_forecasts, rolling_backtest

contains

  subroutine residual_quantiles(residuals,p,distribution,weights,q,status)
    real(dp), intent(in) :: residuals(:,:),p
    character(len=*), intent(in) :: distribution
    real(dp), intent(in), optional :: weights(:)
    real(dp), allocatable, intent(out) :: q(:)
    integer, intent(out) :: status
    real(dp), allocatable :: x(:)
    real(dp) :: kurt,df
    integer :: n,j
    n=size(residuals,2)
    if(p<=0.0_dp .or. p>=1.0_dp)then;status=bekk_invalid_input;return;end if
    if(present(weights))then
      if(size(weights)/=n)then;status=bekk_invalid_input;return;end if
      allocate(q(1),x(size(residuals,1)));x=matmul(residuals,weights)
      select case(trim(distribution))
      case('normal');q(1)=normal_quantile(1.0_dp-p)
      case('empirical');q(1)=empirical_quantile(x,1.0_dp-p)
      case('t')
        kurt=sample_kurtosis(x)-3.0_dp
        if(kurt<=0.0_dp)then;df=100.0_dp;else;df=max(4.001_dp,6.0_dp/kurt+4.0_dp);end if
        q(1)=student_t_quantile(1.0_dp-p,df)/sqrt(df/(df-2.0_dp))
      case default;status=bekk_invalid_input;return
      end select
    else
      allocate(q(n))
      do j=1,n
        select case(trim(distribution))
        case('normal');q(j)=normal_quantile(1.0_dp-p)
        case('empirical');q(j)=empirical_quantile(residuals(:,j),1.0_dp-p)
        case('t')
          kurt=sample_kurtosis(residuals(:,j))-3.0_dp
          if(kurt<=0.0_dp)then;df=100.0_dp;else;df=max(4.001_dp,6.0_dp/kurt+4.0_dp);end if
          q(j)=student_t_quantile(1.0_dp-p,df)/sqrt(df/(df-2.0_dp))
        case default;status=bekk_invalid_input;return
        end select
      end do
    end if
    status=bekk_ok
  end subroutine residual_quantiles

  subroutine var_bekk_fit(fit,p,result,portfolio_weights,distribution)
    type(bekk_fit_result), intent(in) :: fit
    real(dp), intent(in) :: p
    type(bekk_var_result), intent(out) :: result
    real(dp), intent(in), optional :: portfolio_weights(:)
    character(len=*), intent(in), optional :: distribution
    character(len=16) :: dist
    real(dp), allocatable :: q(:)
    integer :: t,n,i,j,st
    dist='empirical';if(present(distribution))dist=distribution
    t=size(fit%h,3);n=size(fit%h,1)
    call residual_quantiles(fit%residuals,p,trim(dist),portfolio_weights,q,st)
    if(st/=bekk_ok)then;result%status=st;return;end if
    result%quantiles=q
    if(present(portfolio_weights))then
      allocate(result%value(t,1))
      do i=1,t
        result%value(i,1)=q(1)*sqrt(max(dot_product(portfolio_weights, &
          matmul(fit%h(:,:,i),portfolio_weights)),0.0_dp))
      end do
    else
      allocate(result%value(t,n))
      do i=1,t
        do j=1,n
          result%value(i,j)=q(j)*sqrt(max(fit%h(j,j,i),0.0_dp))
        end do
      end do
    end if
    result%status=bekk_ok
  end subroutine var_bekk_fit

  subroutine var_bekk_forecast(fit,forecast,p,result,portfolio_weights,distribution)
    type(bekk_fit_result), intent(in) :: fit
    type(bekk_forecast_result), intent(in) :: forecast
    real(dp), intent(in) :: p
    type(bekk_var_result), intent(out) :: result
    real(dp), intent(in), optional :: portfolio_weights(:)
    character(len=*), intent(in), optional :: distribution
    character(len=16) :: dist
    real(dp), allocatable :: q(:)
    integer :: t,n,i,j,st
    dist='empirical';if(present(distribution))dist=distribution
    t=size(forecast%h,3);n=size(forecast%h,1)
    call residual_quantiles(fit%residuals,p,trim(dist),portfolio_weights,q,st)
    if(st/=bekk_ok)then;result%status=st;return;end if
    result%quantiles=q
    if(present(portfolio_weights))then
      allocate(result%value(t,1),result%lower(t,1),result%upper(t,1))
      do i=1,t
        result%value(i,1)=q(1)*sqrt(max(dot_product(portfolio_weights, &
          matmul(forecast%h(:,:,i),portfolio_weights)),0.0_dp))
        result%lower(i,1)=q(1)*sqrt(max(dot_product(portfolio_weights, &
          matmul(forecast%covariance_lower(:,:,i),portfolio_weights)),0.0_dp))
        result%upper(i,1)=q(1)*sqrt(max(dot_product(portfolio_weights, &
          matmul(forecast%covariance_upper(:,:,i),portfolio_weights)),0.0_dp))
      end do
    else
      allocate(result%value(t,n),result%lower(t,n),result%upper(t,n))
      do i=1,t
        do j=1,n
          result%value(i,j)=q(j)*sqrt(max(forecast%h(j,j,i),0.0_dp))
          result%lower(i,j)=q(j)*sqrt(max(forecast%covariance_lower(j,j,i),0.0_dp))
          result%upper(i,j)=q(j)*sqrt(max(forecast%covariance_upper(j,j,i),0.0_dp))
        end do
      end do
    end if
    result%status=bekk_ok
  end subroutine var_bekk_forecast

  subroutine coverage_tests(returns,var,p,kupiec_stat,kupiec_pvalue, &
      christoffersen_stat,christoffersen_pvalue)
    real(dp), intent(in) :: returns(:),var(:),p
    real(dp), intent(out) :: kupiec_stat,kupiec_pvalue,christoffersen_stat,christoffersen_pvalue
    integer :: n,x,i,n00,n01,n10,n11
    real(dp) :: alpha,phat,p01,p11,ll0,ll1,lrind
    logical, allocatable :: hit(:)
    n=size(returns);allocate(hit(n));hit=returns<var;x=count(hit);alpha=1.0_dp-p
    phat=real(x,dp)/real(max(1,n),dp)
    ll0=real(n-x,dp)*safe_log(1.0_dp-alpha)+real(x,dp)*safe_log(alpha)
    ll1=real(n-x,dp)*safe_log(1.0_dp-phat)+real(x,dp)*safe_log(phat)
    kupiec_stat=max(0.0_dp,-2.0_dp*(ll0-ll1));kupiec_pvalue=1.0_dp-chi_square_cdf(kupiec_stat,1)
    n00=0;n01=0;n10=0;n11=0
    do i=2,n
      if(.not.hit(i-1) .and. .not.hit(i))n00=n00+1
      if(.not.hit(i-1) .and. hit(i))n01=n01+1
      if(hit(i-1) .and. .not.hit(i))n10=n10+1
      if(hit(i-1) .and. hit(i))n11=n11+1
    end do
    p01=real(n01,dp)/real(max(1,n00+n01),dp)
    p11=real(n11,dp)/real(max(1,n10+n11),dp)
    ll0=real(n00+n10,dp)*safe_log(1.0_dp-phat)+real(n01+n11,dp)*safe_log(phat)
    ll1=real(n00,dp)*safe_log(1.0_dp-p01)+real(n01,dp)*safe_log(p01) &
      +real(n10,dp)*safe_log(1.0_dp-p11)+real(n11,dp)*safe_log(p11)
    lrind=max(0.0_dp,-2.0_dp*(ll0-ll1))
    christoffersen_stat=kupiec_stat+lrind
    christoffersen_pvalue=1.0_dp-chi_square_cdf(christoffersen_stat,2)
  contains
    pure real(dp) function safe_log(x)
      real(dp), intent(in) :: x
      safe_log=log(max(x,tiny(1.0_dp)))
    end function safe_log
  end subroutine coverage_tests

  subroutine backtest_forecasts(returns,var,p,result)
    real(dp), intent(in) :: returns(:,:),var(:,:),p
    type(bekk_backtest_result), intent(out) :: result
    integer :: n,m,j
    if(any(shape(returns)/=shape(var)))then;result%status=bekk_invalid_input;return;end if
    n=size(returns,1);m=size(returns,2);result%returns=returns;result%var=var
    allocate(result%hit_rate(m),result%kupiec_stat(m),result%kupiec_pvalue(m), &
      result%christoffersen_stat(m),result%christoffersen_pvalue(m))
    do j=1,m
      result%hit_rate(j)=real(count(returns(:,j)<var(:,j)),dp)/real(n,dp)
      call coverage_tests(returns(:,j),var(:,j),p,result%kupiec_stat(j), &
        result%kupiec_pvalue(j),result%christoffersen_stat(j), &
        result%christoffersen_pvalue(j))
    end do
    result%status=bekk_ok
  end subroutine backtest_forecasts

  subroutine rolling_backtest(data,spec,window_length,p,n_ahead,result, &
      portfolio_weights,distribution,max_iter)
    real(dp), intent(in) :: data(:,:)
    type(bekk_spec_type), intent(in) :: spec
    integer, intent(in) :: window_length,n_ahead
    real(dp), intent(in) :: p
    type(bekk_backtest_result), intent(out) :: result
    real(dp), intent(in), optional :: portfolio_weights(:)
    character(len=*), intent(in), optional :: distribution
    integer, intent(in), optional :: max_iter
    type(bekk_fit_result) :: fit
    type(bekk_forecast_result) :: fc
    type(bekk_var_result) :: vr
    real(dp), allocatable :: var(:,:),ret(:,:)
    character(len=16) :: dist
    integer :: t,n,nout,i,h,j,m,mi
    t=size(data,1);n=size(data,2);dist='empirical';if(present(distribution))dist=distribution
    mi=50;if(present(max_iter))mi=max_iter
    if(window_length<2 .or. window_length>=t .or. n_ahead<1)then
      result%status=bekk_invalid_input;return
    end if
    nout=t-window_length;m=merge(1,n,present(portfolio_weights))
    allocate(var(nout,m),ret(nout,m));i=1
    do while(i<=nout)
      h=min(n_ahead,nout-i+1)
      call bekk_fit(spec,data(i:i+window_length-1,:),fit,max_iter=mi)
      if(fit%status/=bekk_ok .and. fit%status/=bekk_no_convergence)then;result%status=fit%status;return;end if
      call forecast_bekk(fit,h,fc)
      call var_bekk_forecast(fit,fc,p,vr,portfolio_weights,trim(dist))
      var(i:i+h-1,:)=vr%value
      if(present(portfolio_weights))then
        do j=1,h;ret(i+j-1,1)=dot_product(data(window_length+i+j-1,:),portfolio_weights);end do
      else
        ret(i:i+h-1,:)=data(window_length+i:window_length+i+h-1,:)
      end if
      i=i+h
    end do
    call backtest_forecasts(ret,var,p,result)
  end subroutine rolling_backtest

end module bekks_risk
