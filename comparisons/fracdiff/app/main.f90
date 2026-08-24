program compare_fracdiff
 use fracdiff, only: dp,diffseries,fd_gph,fd_sperio,fractional_d_estimate
 use comparison_asset_data, only: asset_price_data,asset_return_data,read_asset_prices_binary,log_returns,asset_data_ok
 implicit none
 integer,parameter::n=4096; real(dp)::x(n),dx(n),v,t0,t1
 real(dp),allocatable::spy(:),spy_squared(:),market_dx(:)
 integer::i,j,u,st,market_status,spy_column
 type(fractional_d_estimate)::est; character(512)::out
 type(asset_price_data)::prices
 type(asset_return_data)::market_returns
 character(256)::market_message
 call get_command_argument(1,out); if(len_trim(out)==0)out='fortran_results.csv'
 call read_asset_prices_binary('../../asset_class_etf_prices.bin',prices,market_status,market_message)
 if(market_status/=asset_data_ok)error stop trim(market_message)
 market_returns=log_returns(prices)
 spy_column=find_symbol(prices%symbols,'SPY')
 if(spy_column==0)error stop 'SPY column not found'
 allocate(spy(size(market_returns%returns,1)),spy_squared(size(market_returns%returns,1)), &
  market_dx(size(market_returns%returns,1)))
 spy=market_returns%returns(:,spy_column)
 spy_squared=(spy-sum(spy)/real(size(spy),dp))**2
 open(newunit=u,file=trim(out),status='replace'); write(u,'(a)')'case,value,seconds,abs_tol,rel_tol'
 do i=1,n; x(i)=sin(.017_dp*i)+.35_dp*cos(.0031_dp*i)+.0002_dp*i; end do
 call cpu_time(t0)
 do j=1,100
  call diffseries(x,.25_dp,dx,st); v=wcheck(dx)
 end do
 call cpu_time(t1); call emit('diffseries_d025',v,t1-t0,1e-5_dp)
 call cpu_time(t0)
 do j=1,100
  call diffseries(x,.70_dp,dx,st); v=wcheck(dx)
 end do
 call cpu_time(t1); call emit('diffseries_d070',v,t1-t0,1e-5_dp)
 call cpu_time(t0); do j=1,30; est=fd_gph(x); v=est%d; end do; call cpu_time(t1); call emit('gph_d',v,t1-t0,1e-7_dp)
 call cpu_time(t0); do j=1,30; est=fd_sperio(x); v=est%d; end do; call cpu_time(t1); call emit('sperio_d',v,t1-t0,1e-7_dp)
 call cpu_time(t0); do j=1,30; call diffseries(spy,.25_dp,market_dx,st); v=wcheck(market_dx); end do
 call cpu_time(t1); call emit('etf_spy_diffseries_d025',v,t1-t0,1e-7_dp)
 call cpu_time(t0); do j=1,10; est=fd_gph(spy); v=est%d; end do
 call cpu_time(t1); call emit('etf_spy_gph_d',v,t1-t0,1e-7_dp)
 call cpu_time(t0); do j=1,10; est=fd_sperio(spy); v=est%d; end do
 call cpu_time(t1); call emit('etf_spy_sperio_d',v,t1-t0,1e-7_dp)
 call cpu_time(t0); do j=1,10; est=fd_gph(spy_squared); v=est%d; end do
 call cpu_time(t1); call emit('etf_spy_squared_gph_d',v,t1-t0,1e-7_dp)
 call cpu_time(t0); do j=1,10; est=fd_sperio(spy_squared); v=est%d; end do
 call cpu_time(t1); call emit('etf_spy_squared_sperio_d',v,t1-t0,1e-7_dp)
 close(u)
contains
 real(dp) function wcheck(q)
  real(dp),intent(in)::q(:); integer::k
  wcheck=sum(q*[(real(k,dp),k=1,size(q))])
 end function
 subroutine emit(name,value,secs,atol)
  character(*),intent(in)::name; real(dp),intent(in)::value,secs,atol
  write(u,'(a,",",es24.16,",",es16.8,",",es12.4,",",es12.4)')trim(name),value,secs,atol,1e-8_dp
 end subroutine
 integer function find_symbol(symbols,target)
  character(*),intent(in)::symbols(:),target; integer::k
  find_symbol=0
  do k=1,size(symbols)
   if(trim(symbols(k))==trim(target))then; find_symbol=k; return; end if
  end do
 end function
end program
