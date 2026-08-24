program compare_rugarch
 use rugarch, only: dp, distribution_pdf, distribution_cdf, distribution_quantile, &
  dist_norm,dist_std,dist_ged,dist_snorm,dist_sstd,dist_sged,dist_jsu, &
  garch_spec,make_garch_spec,garch_filter,garch_log_likelihood,model_sgarch,model_gjrgarch
 use comparison_asset_data, only: asset_price_data,asset_return_data, &
  read_asset_prices_binary,log_returns,asset_data_ok
 implicit none
 integer,parameter::n=1001,nd=7
 integer,parameter::kinds(nd)=[dist_norm,dist_std,dist_ged,dist_snorm,dist_sstd,dist_sged,dist_jsu]
 character(5),parameter::names(nd)=[character(5)::'norm','std','ged','snorm','sstd','sged','jsu']
 real(dp),parameter::shapes(nd)=[5._dp,7._dp,1.6_dp,5._dp,7._dp,1.6_dp,1.8_dp]
 real(dp),parameter::skews(nd)=[1._dp,1._dp,1._dp,1.3_dp,1.3_dp,1.3_dp,.5_dp]
 real(dp)::x(n),p(n),y(n),t0,t1,value
 real(dp),allocatable::spy(:),market_sigma(:),market_residuals(:)
 integer::i,j,k,u,reps,market_status,spy_column
 character(512)::out
 character(256)::market_message
 type(asset_price_data)::prices
 type(asset_return_data)::market_returns
 type(garch_spec)::market_spec
 call get_command_argument(1,out); if(len_trim(out)==0)out='fortran_results.csv'
 call read_asset_prices_binary('../../asset_class_etf_prices.bin',prices,market_status,market_message)
 if(market_status/=asset_data_ok)error stop trim(market_message)
 market_returns=log_returns(prices);spy_column=find_symbol(prices%symbols,'SPY')
 if(spy_column==0)error stop 'SPY column not found'
 allocate(spy(size(market_returns%returns,1)),market_sigma(size(market_returns%returns,1)), &
  market_residuals(size(market_returns%returns,1)))
 spy=100._dp*market_returns%returns(:,spy_column)
 open(newunit=u,file=trim(out),status='replace'); write(u,'(a)')'case,value,seconds,abs_tol,rel_tol'
 do i=1,n; x(i)=-4._dp+8._dp*real(i-1,dp)/real(n-1,dp); p(i)=.001_dp+.998_dp*real(i-1,dp)/real(n-1,dp); end do
 do k=1,nd
  reps=1000
  call cpu_time(t0)
  do j=1,reps
   do i=1,n
    y(i)=distribution_pdf(x(i),kinds(k),shapes(k),skews(k))
   end do
   value=wcheck(y)
  end do
  call cpu_time(t1); call emit(trim(names(k))//'_density',value,t1-t0)
  call cpu_time(t0)
  do j=1,reps
   do i=1,n
    y(i)=distribution_cdf(x(i),kinds(k),shapes(k),skews(k))
   end do
   value=wcheck(y)
  end do
  call cpu_time(t1); call emit(trim(names(k))//'_cdf',value,t1-t0)
  reps=500; if(k==1 .or. k==4)reps=2000
  call cpu_time(t0)
  do j=1,reps
   do i=1,n
    y(i)=distribution_quantile(p(i),kinds(k),shapes(k),skews(k))
   end do
   value=wcheck(y)
  end do
  call cpu_time(t1); call emit(trim(names(k))//'_quantile',value,t1-t0)
 end do
 market_spec=make_garch_spec(1,1,model_sgarch,dist_norm)
 market_spec%mean=.04_dp;market_spec%omega=.02_dp
 market_spec%alpha=[.08_dp];market_spec%beta=[.90_dp]
 reps=20;call cpu_time(t0);do j=1,reps
  call garch_filter(spy,market_spec,market_residuals,market_sigma)
  value=wcheck(market_sigma)
 end do;call cpu_time(t1);call emit('etf_spy_sgarch_sigma',value,t1-t0)
 call cpu_time(t0);do j=1,reps
  value=garch_log_likelihood(spy,market_spec)
 end do;call cpu_time(t1);call emit('etf_spy_sgarch_loglik',value,t1-t0)
 market_spec%model=model_gjrgarch;market_spec%gamma=[.04_dp]
 call cpu_time(t0);do j=1,reps
  call garch_filter(spy,market_spec,market_residuals,market_sigma)
  value=wcheck(market_sigma)
 end do;call cpu_time(t1);call emit('etf_spy_gjrgarch_sigma',value,t1-t0)
 call cpu_time(t0);do j=1,reps
  value=garch_log_likelihood(spy,market_spec)
 end do;call cpu_time(t1);call emit('etf_spy_gjrgarch_loglik',value,t1-t0)
 close(u)
contains
 real(dp) function wcheck(q)
  real(dp),intent(in)::q(:); integer::ii
  wcheck=sum(q*[(real(ii,dp),ii=1,size(q))])
 end function
 subroutine emit(name,v,s)
  character(*),intent(in)::name; real(dp),intent(in)::v,s
  write(u,'(a,",",es24.16,",",es16.8,",",es12.4,",",es12.4)')trim(name),v,s,5e-5_dp,2e-7_dp
 end subroutine
 integer function find_symbol(symbols,target)
  character(*),intent(in)::symbols(:),target;integer::ii
  find_symbol=0
  do ii=1,size(symbols)
   if(trim(symbols(ii))==trim(target))then;find_symbol=ii;return;end if
  end do
 end function
end program
