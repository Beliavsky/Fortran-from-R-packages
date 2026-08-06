! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_portfolio
  use tvmvp_kinds, only : dp
  use tvmvp_status, only : tvmvp_error, clear_error, set_error, tvmvp_invalid_input
  use tvmvp_types, only : local_pca_result, factor_selection_result, poet_result, portfolio_result, &
                          portfolio_prediction_result, expanding_window_result, performance_metrics
  use tvmvp_linalg, only : solve_linear, solve_matrix, sample_mean, sample_variance
  use tvmvp_kernels, only : kernel_function
  use tvmvp_pca, only : silverman, determine_factors, local_pca_all
  use tvmvp_poet, only : estimate_residual_cov_poet_local
  use tvmvp_forecast, only : comp_expected_returns
  implicit none
  private
  public :: time_varying_cov, predict_portfolio, expanding_tvmvp
  public :: minimum_variance_portfolio, maximum_sharpe_portfolio, constrained_minimum_variance
contains
  subroutine default_rho_grid(grid)
    real(dp), allocatable, intent(out) :: grid(:)
    integer :: i
    allocate(grid(30))
    do i=1,30
      grid(i)=0.005_dp+real(i-1,dp)*(2.0_dp-0.005_dp)/29.0_dp
    end do
  end subroutine default_rho_grid

  subroutine time_varying_cov(returns,covariance,err,max_factors,m,bandwidth,kernel,m0,rho_grid, &
                              floor_value,epsilon2,full_result,source_compatible_boundary)
    real(dp), intent(in) :: returns(:,:)
    real(dp), allocatable, intent(out) :: covariance(:,:)
    type(tvmvp_error), intent(out) :: err
    integer, intent(in), optional :: max_factors,m,m0
    real(dp), intent(in), optional :: bandwidth,rho_grid(:),floor_value,epsilon2
    procedure(kernel_function), optional :: kernel
    type(poet_result), intent(out), optional :: full_result
    logical, intent(in), optional :: source_compatible_boundary
    type(factor_selection_result) :: selection
    type(local_pca_result) :: local_result
    type(poet_result) :: poet
    real(dp), allocatable :: grid(:)
    real(dp) :: h,floorv,eps
    integer :: nf,maxf,mm0
    logical :: compat
    call clear_error(err)
    if (size(returns,1)<4 .or. size(returns,2)<1) then
      allocate(covariance(0,0)); call set_error(err,tvmvp_invalid_input,'insufficient return data'); return
    end if
    compat=.true.; if (present(source_compatible_boundary)) compat=source_compatible_boundary
    h=silverman(returns); if (present(bandwidth)) h=bandwidth
    if (present(m)) then
      nf=m
    else
      maxf=min(3,min(size(returns,1),size(returns,2)))
      if (present(max_factors)) maxf=max_factors
      if (present(kernel)) then
        call determine_factors(returns,maxf,selection,h,kernel,compat)
      else
        call determine_factors(returns,maxf,selection,h,source_compatible_boundary=compat)
      end if
      if (selection%error%failed()) then
        allocate(covariance(0,0)); err=selection%error; return
      end if
      nf=selection%optimal_m
    end if
    if (present(kernel)) then
      call local_pca_all(returns,h,nf,local_result,kernel,compat)
    else
      call local_pca_all(returns,h,nf,local_result,source_compatible_boundary=compat)
    end if
    if (local_result%error%failed()) then
      allocate(covariance(0,0)); err=local_result%error; return
    end if
    mm0=10; if (present(m0)) mm0=m0
    floorv=1.0e-12_dp; if (present(floor_value)) floorv=floor_value
    eps=1.0e-6_dp; if (present(epsilon2)) eps=epsilon2
    if (present(rho_grid)) then
      allocate(grid(size(rho_grid))); grid=rho_grid
    else
      call default_rho_grid(grid)
    end if
    call estimate_residual_cov_poet_local(local_result,returns,poet,mm0,grid,floorv,eps)
    if (poet%error%failed()) then
      allocate(covariance(0,0)); err=poet%error; return
    end if
    allocate(covariance(size(poet%total_cov,1),size(poet%total_cov,2))); covariance=poet%total_cov
    if (present(full_result)) full_result=poet
  end subroutine time_varying_cov

  subroutine fill_portfolio_metrics(weights,mu,sigma,horizon,result)
    real(dp), intent(in) :: weights(:),mu(:),sigma(:,:)
    integer, intent(in) :: horizon
    type(portfolio_result), intent(inout) :: result
    real(dp) :: daily_risk
    result%expected_return=dot_product(weights,mu)*real(horizon,dp)
    daily_risk=sqrt(max(0.0_dp,dot_product(weights,matmul(sigma,weights))))
    result%risk=daily_risk*sqrt(real(horizon,dp))
    if (daily_risk>tiny(1.0_dp)) then
      result%sharpe=dot_product(weights,mu)/daily_risk
    else
      result%sharpe=0.0_dp
    end if
    result%available=.true.
  end subroutine fill_portfolio_metrics

  subroutine minimum_variance_portfolio(sigma,mu,horizon,result)
    real(dp), intent(in) :: sigma(:,:),mu(:)
    integer, intent(in) :: horizon
    type(portfolio_result), intent(out) :: result
    real(dp), allocatable :: ones(:),raw(:)
    call clear_error(result%error)
    if (size(sigma,1)/=size(sigma,2) .or. size(mu)/=size(sigma,1) .or. horizon<1) then
      call set_error(result%error,tvmvp_invalid_input,'invalid minimum-variance input'); return
    end if
    allocate(ones(size(mu))); ones=1.0_dp
    call solve_linear(sigma,ones,raw,result%error); if (result%error%failed()) return
    if (abs(sum(raw))<=tiny(1.0_dp)) then
      call set_error(result%error,tvmvp_invalid_input,'minimum-variance weights cannot be normalized'); return
    end if
    allocate(result%weights(size(mu))); result%weights=raw/sum(raw)
    call fill_portfolio_metrics(result%weights,mu,sigma,horizon,result)
  end subroutine minimum_variance_portfolio

  subroutine maximum_sharpe_portfolio(sigma,mu,horizon,result)
    real(dp), intent(in) :: sigma(:,:),mu(:)
    integer, intent(in) :: horizon
    type(portfolio_result), intent(out) :: result
    real(dp), allocatable :: raw(:)
    call clear_error(result%error)
    call solve_linear(sigma,mu,raw,result%error); if (result%error%failed()) return
    if (abs(sum(raw))<=tiny(1.0_dp)) then
      call set_error(result%error,tvmvp_invalid_input,'maximum-Sharpe weights cannot be normalized'); return
    end if
    allocate(result%weights(size(mu))); result%weights=raw/sum(raw)
    call fill_portfolio_metrics(result%weights,mu,sigma,horizon,result)
  end subroutine maximum_sharpe_portfolio

  subroutine constrained_minimum_variance(sigma,mu,horizon,target_return,result)
    real(dp), intent(in) :: sigma(:,:),mu(:),target_return
    integer, intent(in) :: horizon
    type(portfolio_result), intent(out) :: result
    real(dp), allocatable :: a(:,:),sinva(:,:),middle(:,:),b(:),coef(:)
    integer :: p
    call clear_error(result%error)
    p=size(mu)
    if (size(sigma,1)/=p .or. size(sigma,2)/=p .or. horizon<1) then
      call set_error(result%error,tvmvp_invalid_input,'invalid constrained portfolio input'); return
    end if
    allocate(a(p,2)); a(:,1)=1.0_dp; a(:,2)=mu
    call solve_matrix(sigma,a,sinva,result%error); if (result%error%failed()) return
    allocate(middle(2,2)); middle=matmul(transpose(a),sinva)
    allocate(b(2)); b=[1.0_dp,target_return/real(horizon,dp)]
    call solve_linear(middle,b,coef,result%error); if (result%error%failed()) return
    allocate(result%weights(p)); result%weights=matmul(sinva,coef)
    call fill_portfolio_metrics(result%weights,mu,sigma,horizon,result)
  end subroutine constrained_minimum_variance

  subroutine predict_portfolio(returns,horizon,result,max_factors,m,kernel,min_return,compute_max_sharpe, &
                               risk_free,m0,rho_grid,floor_value,epsilon2,source_compatible_boundary)
    real(dp), intent(in) :: returns(:,:)
    integer, intent(in) :: horizon
    type(portfolio_prediction_result), intent(out) :: result
    integer, intent(in), optional :: max_factors,m,m0
    procedure(kernel_function), optional :: kernel
    real(dp), intent(in), optional :: min_return,risk_free,rho_grid(:),floor_value,epsilon2
    logical, intent(in), optional :: compute_max_sharpe,source_compatible_boundary
    real(dp), allocatable :: mu(:),cov(:,:)
    type(tvmvp_error) :: err
    logical :: do_max,compat
    integer :: maxf,mm0,nf
    call clear_error(result%error)
    maxf=3; if (present(max_factors)) maxf=max_factors
    mm0=10; if (present(m0)) mm0=m0
    compat=.true.; if (present(source_compatible_boundary)) compat=source_compatible_boundary
    if (present(kernel)) then
      call time_varying_cov(returns,cov,err,maxf,m,kernel=kernel,m0=mm0,rho_grid=rho_grid, &
        floor_value=floor_value,epsilon2=epsilon2,source_compatible_boundary=compat)
    else
      call time_varying_cov(returns,cov,err,maxf,m,m0=mm0,rho_grid=rho_grid, &
        floor_value=floor_value,epsilon2=epsilon2,source_compatible_boundary=compat)
    end if
    if (err%failed()) then; result%error=err; return; end if
    call comp_expected_returns(returns,horizon,mu,err)
    if (err%failed()) then; result%error=err; return; end if
    if (present(risk_free)) mu=mu-risk_free
    call minimum_variance_portfolio(cov,mu,horizon,result%minimum_variance)
    if (result%minimum_variance%error%failed()) then; result%error=result%minimum_variance%error; return; end if
    do_max=.false.; if (present(compute_max_sharpe)) do_max=compute_max_sharpe
    if (do_max) then
      call maximum_sharpe_portfolio(cov,mu,horizon,result%maximum_sharpe)
      if (result%maximum_sharpe%error%failed()) then; result%error=result%maximum_sharpe%error; return; end if
    end if
    if (present(min_return)) then
      call constrained_minimum_variance(cov,mu,horizon,min_return,result%return_constrained)
      if (result%return_constrained%error%failed()) then; result%error=result%return_constrained%error; return; end if
    end if
    allocate(result%covariance(size(cov,1),size(cov,2)),result%expected_asset_returns(size(mu)))
    result%covariance=cov; result%expected_asset_returns=mu
    result%bandwidth=silverman(returns)
    if (present(m)) then
      result%factors=m
    else
      nf=0
      result%factors=nf
    end if
  end subroutine predict_portfolio

  subroutine calculate_metrics(x,annual_factor,metrics)
    real(dp), intent(in) :: x(:),annual_factor
    type(performance_metrics), intent(out) :: metrics
    metrics%cumulative_excess_return=sum(x)
    metrics%mean_excess_return=sample_mean(x)
    metrics%standard_deviation=sqrt(sample_variance(x))
    if (metrics%standard_deviation>tiny(1.0_dp)) &
      metrics%sharpe=metrics%mean_excess_return/metrics%standard_deviation
    metrics%annualized_mean=metrics%mean_excess_return*annual_factor*annual_factor
    metrics%annualized_standard_deviation=metrics%standard_deviation*annual_factor
  end subroutine calculate_metrics

  subroutine expanding_tvmvp(returns,initial_window,rebal_period,max_factors,result,return_type,kernel,rf,m0, &
                             rho_grid,floor_value,epsilon2,source_compatible_boundary,source_compatible_expanding)
    real(dp), intent(in) :: returns(:,:)
    integer, intent(in) :: initial_window,rebal_period,max_factors
    type(expanding_window_result), intent(out) :: result
    character(len=*), intent(in), optional :: return_type
    procedure(kernel_function), optional :: kernel
    real(dp), intent(in), optional :: rf(:),rho_grid(:),floor_value,epsilon2
    integer, intent(in), optional :: m0
    logical, intent(in), optional :: source_compatible_boundary,source_compatible_expanding
    integer, allocatable :: dates(:)
    real(dp), allocatable :: est_data(:,:),cov(:,:),ones(:),raw(:),w(:),rfvec(:)
    type(factor_selection_result) :: selection
    type(tvmvp_error) :: err
    real(dp) :: h,annual_factor
    integer :: n,p,rt,l,reb_t,hold_end,hold_len,pos,days_since,m,mm0
    logical :: update_m,compat_boundary,compat_expanding
    character(len=16) :: freq
    call clear_error(result%error)
    n=size(returns,1); p=size(returns,2)
    if (initial_window<4 .or. initial_window>=n .or. rebal_period<1 .or. max_factors<1) then
      call set_error(result%error,tvmvp_invalid_input,'invalid expanding-window input'); return
    end if
    compat_boundary=.true.; if (present(source_compatible_boundary)) compat_boundary=source_compatible_boundary
    compat_expanding=.true.; if (present(source_compatible_expanding)) compat_expanding=source_compatible_expanding
    mm0=10; if (present(m0)) mm0=m0
    rt=(n-initial_window+rebal_period-1)/rebal_period
    allocate(dates(rt))
    do l=1,rt
      dates(l)=initial_window+1+(l-1)*rebal_period
    end do
    allocate(result%rebalance_dates(rt),result%holding_lengths(rt),result%weights(p,rt))
    allocate(result%tvmvp_returns(n-initial_window),result%equal_returns(n-initial_window),rfvec(n-initial_window))
    result%rebalance_dates=dates; result%weights=0.0_dp
    if (present(rf)) then
      if (size(rf)==1) then
        rfvec=rf(1)
      else if (size(rf)==n-initial_window) then
        rfvec=rf
      else
        call set_error(result%error,tvmvp_invalid_input,'risk-free vector has the wrong length'); return
      end if
    else
      rfvec=0.0_dp
    end if
    pos=1; days_since=0; m=0
    do l=1,rt
      reb_t=dates(l)
      update_m=(l==1 .or. days_since>=252)
      if (update_m) then
        est_data=returns(1:reb_t-1,:)
        h=silverman(est_data)
        if (present(kernel)) then
          call determine_factors(est_data,min(max_factors,min(size(est_data,1),p)),selection,h,kernel,compat_boundary)
        else
          call determine_factors(est_data,min(max_factors,min(size(est_data,1),p)),selection,h, &
                                 source_compatible_boundary=compat_boundary)
        end if
        if (selection%error%failed()) then; result%error=selection%error; return; end if
        m=selection%optimal_m; days_since=0
      else
        days_since=days_since+rebal_period
        if (.not.compat_expanding) then
          est_data=returns(1:reb_t-1,:)
          h=silverman(est_data)
        end if
      end if
      if (present(kernel)) then
        call time_varying_cov(est_data,cov,err,m=m,bandwidth=h,kernel=kernel,m0=mm0,rho_grid=rho_grid, &
          floor_value=floor_value,epsilon2=epsilon2,source_compatible_boundary=compat_boundary)
      else
        call time_varying_cov(est_data,cov,err,m=m,bandwidth=h,m0=mm0,rho_grid=rho_grid, &
          floor_value=floor_value,epsilon2=epsilon2,source_compatible_boundary=compat_boundary)
      end if
      if (err%failed()) then; result%error=err; return; end if
      allocate(ones(p)); ones=1.0_dp
      call solve_linear(cov,ones,raw,err); if (err%failed()) then; result%error=err; return; end if
      allocate(w(p)); w=raw/sum(raw); result%weights(:,l)=w
      hold_end=min(reb_t+rebal_period-1,n); hold_len=hold_end-reb_t+1
      result%holding_lengths(l)=hold_len
      result%tvmvp_returns(pos:pos+hold_len-1)=matmul(returns(reb_t:hold_end,:),w)
      result%equal_returns(pos:pos+hold_len-1)=sum(returns(reb_t:hold_end,:),dim=2)/real(p,dp)
      pos=pos+hold_len
      deallocate(ones,raw,w,cov)
    end do
    result%tvmvp_returns=result%tvmvp_returns-rfvec
    result%equal_returns=result%equal_returns-rfvec
    freq='daily'; if (present(return_type)) freq=adjustl(return_type)
    select case(trim(freq))
    case('daily'); annual_factor=sqrt(252.0_dp)
    case('weekly'); annual_factor=sqrt(52.0_dp)
    case('monthly'); annual_factor=sqrt(12.0_dp)
    case default
      call set_error(result%error,tvmvp_invalid_input,'return_type must be daily, weekly, or monthly'); return
    end select
    call calculate_metrics(result%tvmvp_returns,annual_factor,result%tvmvp_metrics)
    call calculate_metrics(result%equal_returns,annual_factor,result%equal_metrics)
  end subroutine expanding_tvmvp
end module tvmvp_portfolio
