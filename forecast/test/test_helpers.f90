program test_helpers
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use forecast, only : dp, dshw_model, decomposition_result, dshw_fit, dshw_forecast, mstl_decompose, &
      ts_clean, ts_outliers, bld_mbb_bootstrap, forecast_result
   implicit none
   real(dp) :: y(96),x(20)
   real(dp),allocatable::samples(:,:),repl(:),cleaned(:)
   integer,allocatable::idx(:)
   integer::i
   type(dshw_model)::dm
   type(decomposition_result)::dec
   type(forecast_result)::fc
   do i=1,96
      y(i)=(15.0_dp+0.02_dp*i)*(1.0_dp+0.10_dp*sin(2.0_dp*acos(-1.0_dp)*i/4.0_dp)) * &
         (1.0_dp+0.06_dp*cos(2.0_dp*acos(-1.0_dp)*i/12.0_dp))
   end do
   dm=dshw_fit(y,4,12,optimize=.false.)
   fc=dshw_forecast(dm,4)
   call check(all(fc%mean>0.0_dp),'DSHW forecasts')
   dec=mstl_decompose(y,[4,12],2)
   call check(maxval(abs(y-dec%trend-sum(dec%seasonal,dim=2)-dec%remainder))<1.0e-10_dp,'MSTL reconstruction')
   x=[(real(i,dp),i=1,20)]
   x(8)=100.0_dp
   x(12)=ieee_value(0.0_dp,ieee_quiet_nan)
   call ts_outliers(x,idx,repl)
   call check(any(idx==8),'outlier detection')
   cleaned=ts_clean(x)
   call check(all(cleaned==cleaned),'cleaning fills missing')
   samples=bld_mbb_bootstrap(y,3,12)
   call check(size(samples,1)==96 .and. size(samples,2)==3,'bootstrap dimensions')
   call check(maxval(abs(samples(:,1)-y))<1e-12_dp,'bootstrap retains original first')
   print '(a)','test_helpers: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      write(*,'(a)')'FAIL: '//trim(msg)
      error stop 1
      end if
   end subroutine
end program
