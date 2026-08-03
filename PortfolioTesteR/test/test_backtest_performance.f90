! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
program test_backtest_performance
  use portfolio_tester
  implicit none
  real(dp)::prices(3,1),weights(3,1)
  real(dp),allocatable::dd(:)
  type(backtest_result)::bt
  type(performance_result)::perf
  prices(:,1)=[100.0_dp,110.0_dp,121.0_dp];weights=1.0_dp
  call run_backtest(prices,weights,1000.0_dp,bt,integer_shares=.false.,frequency=1.0_dp)
  call assert_close(bt%portfolio_value(1),1000.0_dp,1.0e-10_dp,'initial value')
  call assert_close(bt%portfolio_value(2),1100.0_dp,1.0e-10_dp,'second value')
  call assert_close(bt%portfolio_value(3),1210.0_dp,1.0e-10_dp,'third value')
  call assert_close(bt%total_return,0.21_dp,1.0e-12_dp,'total return')
  call analyze_performance(bt,1.0_dp,perf)
  call assert_close(perf%total_return,0.21_dp,1.0e-12_dp,'performance total')
  call calculate_drawdown_series([100.0_dp,120.0_dp,90.0_dp,121.0_dp],dd)
  call assert_close(dd(3),-0.25_dp,1.0e-12_dp,'drawdown')
  call assert_true(calculate_recovery_time(dd)==1,'recovery time')
  print '(a)','test_backtest_performance: PASS'
contains
  subroutine assert_true(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;write(*,'(a)')'FAIL: '//msg;error stop 1;end if
  end subroutine
  subroutine assert_close(x,y,tol,msg)
    real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::msg
    call assert_true(abs(x-y)<=tol,msg)
  end subroutine
end program test_backtest_performance
