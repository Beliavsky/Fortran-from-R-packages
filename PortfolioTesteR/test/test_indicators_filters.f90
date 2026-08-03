! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program test_indicators_filters
  use portfolio_tester
  implicit none
  real(dp),allocatable::prices(:,:),mom(:,:),ma(:,:),rsi(:,:),sel(:,:),w(:,:)
  real(dp),allocatable::vstd(:,:),vrange(:,:),vmad(:,:),vabs(:,:),known(:,:)
  integer,allocatable::counts(:)
  integer::t
  call generate_sample_prices(80,8,prices,123_i8)
  call calc_momentum(prices,12,mom)
  call calc_moving_average(prices,10,ma)
  call calc_rsi(prices,6,rsi)
  call calc_rolling_volatility(prices,12,vstd,annualization=52.0_dp,method=vol_std)
  call calc_rolling_volatility(prices,12,vrange,method=vol_range)
  call calc_rolling_volatility(prices,12,vmad,annualization=52.0_dp,method=vol_mad)
  call calc_rolling_volatility(prices,12,vabs,method=vol_abs_return)
  call filter_top_n(mom,3,sel)
  call selection_counts(sel,counts)
  call weight_equally(sel,w)
  call assert_true(all(.not.is_finite(mom(:12,:))), 'momentum warmup')
  call assert_true(all(is_finite(mom(13:,:))), 'momentum finite')
  call assert_true(all(is_finite(ma(10:,:))), 'moving average finite')
  call assert_true(minval(rsi,mask=is_finite(rsi))>=0.0_dp .and. &
    maxval(rsi,mask=is_finite(rsi))<=100.0_dp,'rsi bounds')
  call assert_true(all(is_finite(vstd(13:,:))), 'standard volatility finite')
  call assert_true(all(is_finite(vrange(13:,:))), 'range volatility finite')
  call assert_true(all(is_finite(vmad(13:,:))), 'MAD volatility finite')
  call assert_true(all(is_finite(vabs(13:,:))), 'absolute-return volatility finite')
  allocate(known(8,1)); known(:,1) = [100.0_dp,101.0_dp,102.01_dp,103.0301_dp, &
    104.060401_dp,105.10100501_dp,106.1520150601_dp,107.213535210701_dp]
  call calc_rolling_volatility(known,4,vstd,method=vol_std)
  call calc_rolling_volatility(known,4,vmad,method=vol_mad)
  call calc_rolling_volatility(known,4,vabs,method=vol_abs_return)
  call assert_close(vstd(8,1),0.0_dp,1.0e-12_dp,'constant-return standard volatility')
  call assert_close(vmad(8,1),0.0_dp,1.0e-12_dp,'constant-return MAD volatility')
  call assert_close(vabs(8,1),0.01_dp,1.0e-12_dp,'constant-return absolute volatility')
  call assert_true(all(counts(:12)==0),'filter warmup empty')
  call assert_true(all(counts(13:)==3),'top n counts')
  do t=13,size(w,1)
    call assert_close(sum(w(t,:)),1.0_dp,1.0e-12_dp,'equal weights sum')
  end do
  print '(a)','test_indicators_filters: PASS'
contains
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol
    character(len=*),intent(in)::msg
    call assert_true(abs(x-y)<=tol,msg)
  end subroutine
end program test_indicators_filters
