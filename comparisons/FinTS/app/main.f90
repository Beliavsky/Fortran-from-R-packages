program compare_fints_market
 use fints, only: dp, acf, AutocorTest, ArchTest, acf_result, test_result
 use comparison_asset_data, only: asset_price_data, asset_return_data, &
  read_asset_prices_binary, log_returns, asset_data_ok
 implicit none
 type(asset_price_data) :: prices
 type(asset_return_data) :: returns
 type(acf_result) :: acf_values
 type(test_result) :: test
 real(dp), allocatable :: spy(:), spy2(:)
 real(dp) :: value, t0, t1
 integer :: status, spy_col, i, j, unit, reps
 character(512) :: output
 character(256) :: message
 call get_command_argument(1,output)
 if(len_trim(output)==0)output='fortran_results.csv'
 call read_asset_prices_binary('../../asset_class_etf_prices.bin',prices,status,message)
 if(status/=asset_data_ok)error stop trim(message)
 returns=log_returns(prices);spy_col=find_symbol(prices%symbols,'SPY')
 if(spy_col==0)error stop 'SPY column not found'
 allocate(spy(size(returns%returns,1)),spy2(size(returns%returns,1)))
 spy=returns%returns(:,spy_col);spy2=(spy-sum(spy)/real(size(spy),dp))**2
 open(newunit=unit,file=trim(output),status='replace')
 write(unit,'(a)')'case,value,seconds,abs_tol,rel_tol'
 call timed_acf('etf_spy_acf20',spy,'correlation')
 call timed_acf('etf_spy_squared_acf20',spy2,'correlation')
 call timed_acf('etf_spy_pacf20',spy,'partial')
 reps=100;call cpu_time(t0)
 do i=1,reps;call AutocorTest(spy,test,lag=20);value=test%statistic;end do
 call cpu_time(t1);call emit('etf_spy_ljung_box_stat',value,t1-t0)
 call cpu_time(t0)
 do i=1,reps;call AutocorTest(spy,test,lag=20);value=test%p_value;end do
 call cpu_time(t1);call emit('etf_spy_ljung_box_pvalue',value,t1-t0)
 reps=20;call cpu_time(t0)
 do i=1,reps;call ArchTest(spy,test,lags=12);value=test%statistic;end do
 call cpu_time(t1);call emit('etf_spy_arch12_stat',value,t1-t0)
 call cpu_time(t0)
 do i=1,reps;call ArchTest(spy,test,lags=12);value=test%p_value;end do
 call cpu_time(t1);call emit('etf_spy_arch12_pvalue',value,t1-t0)
 close(unit)
contains
 subroutine timed_acf(name,x,kind)
  character(*),intent(in)::name,kind
  real(dp),intent(in)::x(:)
  call cpu_time(t0)
  do i=1,100
   call acf(x,acf_values,lag_max=20,acf_type=kind)
   value=sum(acf_values%value*[(real(j,dp),j=1,size(acf_values%value))])
  end do
  call cpu_time(t1);call emit(name,value,t1-t0)
 end subroutine
 subroutine emit(name,result,seconds)
  character(*),intent(in)::name
  real(dp),intent(in)::result,seconds
  write(unit,'(a,",",es25.16e3,",",es16.8,",",es12.4,",",es12.4)') &
   trim(name),result,seconds,1e-9_dp,1e-8_dp
 end subroutine
 integer function find_symbol(symbols,target)
  character(*),intent(in)::symbols(:),target
  integer::k
  find_symbol=0
  do k=1,size(symbols)
   if(trim(symbols(k))==trim(target))then;find_symbol=k;return;end if
  end do
 end function
end program
