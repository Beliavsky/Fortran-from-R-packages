program compare_performanceanalytics
 use performanceanalytics_mod, only: dp, annualized_return, sd_value, &
  annualized_sharpe_ratio, sortino_ratio, max_drawdown, historical_var, &
  historical_es, capm_alpha, capm_beta, capm_beta_bull, correlation_value, tracking_error
 use comparison_asset_data, only: asset_price_data, asset_return_data, &
  read_asset_prices_binary, simple_returns, asset_data_ok
 implicit none
 type(asset_price_data) :: prices
 type(asset_return_data) :: returns
 real(dp), allocatable :: spy(:), efa(:), rf(:)
 real(dp) :: value, t0, t1
 integer :: status, spy_col, efa_col, reps, i, unit
 character(512) :: output
 character(256) :: message

 call get_command_argument(1, output)
 if (len_trim(output)==0) output='fortran_results.csv'
 call read_asset_prices_binary('../../asset_class_etf_prices.bin', prices, status, message)
 if (status/=asset_data_ok) error stop trim(message)
 returns=simple_returns(prices)
 spy_col=find_symbol(prices%symbols,'SPY')
 efa_col=find_symbol(prices%symbols,'EFA')
 if (spy_col==0 .or. efa_col==0) error stop 'required ETF column not found'
 allocate(spy(size(returns%returns,1)),efa(size(returns%returns,1)),rf(size(returns%returns,1)))
 spy=returns%returns(:,spy_col); efa=returns%returns(:,efa_col); rf=0.0_dp
 open(newunit=unit,file=trim(output),status='replace')
 write(unit,'(a)') 'case,value,seconds,abs_tol,rel_tol'
 reps=100
 call timed_scalar('etf_spy_annualized_return', reps, 1)
 call timed_scalar('etf_spy_annualized_sd', reps, 2)
 call timed_scalar('etf_spy_sharpe', reps, 3)
 call timed_scalar('etf_spy_sortino', reps, 4)
 call timed_scalar('etf_spy_max_drawdown', reps, 5)
 call timed_scalar('etf_spy_historical_var95', reps, 6)
 call timed_scalar('etf_spy_historical_es95', reps, 7)
 call timed_scalar('etf_efa_capm_alpha', reps, 8)
 call timed_scalar('etf_efa_capm_beta', reps, 9)
 call timed_scalar('etf_efa_capm_beta_bull', reps, 10)
 call timed_scalar('etf_efa_spy_correlation', reps, 11)
 call timed_scalar('etf_efa_tracking_error', reps, 12)
 close(unit)
contains
 subroutine timed_scalar(name, count, operation)
  character(*), intent(in) :: name
  integer, intent(in) :: count, operation
  call cpu_time(t0)
  do i=1,count
   select case(operation)
   case(1); value=annualized_return(spy,252.0_dp,.true.)
   case(2); value=sd_value(spy)*sqrt(252.0_dp)
   case(3); value=annualized_sharpe_ratio(spy,rf,252.0_dp,.false.)
   case(4); value=sortino_ratio(spy,0.0_dp)
   case(5); value=max_drawdown(spy)
   case(6); value=historical_var(spy,.95_dp)
   case(7); value=historical_es(spy,.95_dp)
   case(8); value=capm_alpha(efa,spy)
   case(9); value=capm_beta(efa,spy)
   case(10); value=capm_beta_bull(efa,spy)
   case(11); value=correlation_value(efa,spy)
   case(12); value=tracking_error(efa,spy,252.0_dp)
   end select
  end do
  call cpu_time(t1)
  call emit(name,value,t1-t0)
 end subroutine
 subroutine emit(name, result, seconds)
  character(*), intent(in) :: name
  real(dp), intent(in) :: result, seconds
  write(unit,'(a,",",es25.16e3,",",es16.8,",",es12.4,",",es12.4)') &
   trim(name),result,seconds,1e-10_dp,1e-8_dp
 end subroutine
 integer function find_symbol(symbols,target)
  character(*), intent(in) :: symbols(:),target
  integer :: j
  find_symbol=0
  do j=1,size(symbols)
   if (trim(symbols(j))==trim(target)) then
    find_symbol=j; return
   end if
  end do
 end function
end program
