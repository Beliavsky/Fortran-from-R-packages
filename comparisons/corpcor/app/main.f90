program compare_corpcor_market
 use corpcor, only: dp, matrix_shrinkage_result, correlation_shrinkage, &
  covariance_shrinkage, inverse_covariance_shrinkage, partial_correlation_shrinkage
 use comparison_asset_data, only: asset_price_data, asset_return_data, &
  read_asset_prices_binary, log_returns, asset_data_ok
 implicit none
 type(asset_price_data)::prices
 type(asset_return_data)::returns
 type(matrix_shrinkage_result)::result
 real(dp)::value,t0,t1
 integer::status,i,unit,reps
 character(512)::output
 character(256)::message
 call get_command_argument(1,output);if(len_trim(output)==0)output='fortran_results.csv'
 call read_asset_prices_binary('../../asset_class_etf_prices.bin',prices,status,message)
 if(status/=asset_data_ok)error stop trim(message)
 returns=log_returns(prices)
 open(newunit=unit,file=trim(output),status='replace')
 write(unit,'(a)')'case,value,seconds,abs_tol,rel_tol'
 reps=30
 call timed_case('etf_cor_shrink_fixed',1,reps)
 call timed_case('etf_cov_shrink_fixed',2,reps)
 call timed_case('etf_invcov_shrink_fixed',3,reps)
 call timed_case('etf_pcor_shrink_fixed',4,reps)
 call timed_case('etf_cov_shrink_estimated',5,10)
 call timed_case('etf_cov_lambda_estimated',6,10)
 call timed_case('etf_var_lambda_estimated',7,10)
 close(unit)
contains
 subroutine timed_case(name,operation,count)
  character(*),intent(in)::name
  integer,intent(in)::operation,count
  call cpu_time(t0)
  do i=1,count
   select case(operation)
   case(1);result=correlation_shrinkage(returns%returns,lambda=.2_dp)
   case(2);result=covariance_shrinkage(returns%returns,lambda=.2_dp,lambda_var=.1_dp)
   case(3);result=inverse_covariance_shrinkage(returns%returns,lambda=.2_dp,lambda_var=.1_dp)
   case(4);result=partial_correlation_shrinkage(returns%returns,lambda=.2_dp)
   case default;result=covariance_shrinkage(returns%returns)
   end select
   select case(operation)
   case(6);value=result%lambda
   case(7);value=result%lambda_var
   case default;value=wcheck(result%value)
   end select
  end do
  call cpu_time(t1);call emit(name,value,t1-t0)
 end subroutine
 real(dp) function wcheck(x)
  real(dp),intent(in)::x(:,:)
  integer::row,column,k
  wcheck=0.0_dp;k=0
  do column=1,size(x,2)
   do row=1,size(x,1)
    k=k+1;wcheck=wcheck+real(k,dp)*x(row,column)
   end do
  end do
 end function
 subroutine emit(name,x,seconds)
  character(*),intent(in)::name
  real(dp),intent(in)::x,seconds
  write(unit,'(a,",",es25.16e3,",",es16.8,",",es12.4,",",es12.4)') &
   trim(name),x,seconds,1e-9_dp,1e-8_dp
 end subroutine
end program
