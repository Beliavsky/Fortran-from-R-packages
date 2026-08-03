program test_backtest
   use pmwr, only : dp, backtest_result, run_backtest
   implicit none
   type(backtest_result) :: bt
   real(dp) :: prices(4,1), targets(4,1)

   prices(:,1) = [100.0_dp,110.0_dp,120.0_dp,115.0_dp]
   targets(:,1) = [1.0_dp,1.0_dp,0.0_dp,0.0_dp]
   call run_backtest(prices, targets, bt, initial_cash=1000.0_dp, lag=1)
   call assert_close(bt%position(2,1), 1.0_dp, 1.0e-12_dp, "position buy")
   call assert_close(bt%wealth(3), 1010.0_dp, 1.0e-12_dp, "wealth mark")
   call assert_close(bt%wealth(4), 1005.0_dp, 1.0e-12_dp, "final wealth")
   if (bt%journal%n /= 2) error stop "journal trade count"

   call run_backtest(prices, reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp],[4,1]), bt, &
                     initial_cash=1000.0_dp, lag=1, convert_weights=.true.)
   call assert_close(bt%position(2,1), 1000.0_dp/100.0_dp, 1.0e-12_dp, "weight conversion")

   print *, "test_backtest: PASS"
contains
   subroutine assert_close(a,b,tol,label)
      real(dp), intent(in) :: a,b,tol
      character(len=*), intent(in) :: label
      if (abs(a-b) > tol) then
         print *, trim(label), a, b
         error stop 1
      end if
   end subroutine assert_close
end program test_backtest
