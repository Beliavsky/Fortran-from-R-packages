program test_portfolio
   use pmwr, only : dp, rebalance_result, unit_price_result, rebalance_portfolio, unit_prices, &
                    dividend_adjust, split_adjust
   implicit none
   type(rebalance_result) :: rb
   type(unit_price_result) :: up
   real(dp), allocatable :: x(:)
   real(dp) :: current(2), target(2), price(2), nav(3), cf(2)
   integer :: idx(2)

   current = [10.0_dp, 20.0_dp]
   target = [0.75_dp, 0.25_dp]
   price = [10.0_dp, 5.0_dp]
   call rebalance_portfolio(current, target, price, rb)
   call assert_close(rb%notional, 200.0_dp, 1.0e-12_dp, "notional")
   call assert_close(rb%target(1), 15.0_dp, 1.0e-12_dp, "target 1")
   call assert_close(rb%target(2), 10.0_dp, 1.0e-12_dp, "target 2")
   call assert_close(rb%turnover, 100.0_dp, 1.0e-12_dp, "turnover")

   call dividend_adjust([100.0_dp,95.0_dp,100.0_dp], [2], [5.0_dp], x, backward=.false.)
   call assert_close(x(2), 100.0_dp, 1.0e-12_dp, "dividend adjusted 2")
   call assert_close(x(3), 100.0_dp*100.0_dp/95.0_dp, 1.0e-12_dp, "dividend adjusted 3")

   call split_adjust([100.0_dp,50.0_dp,55.0_dp], [2], [2.0_dp], x)
   call assert_close(x(1), 50.0_dp, 1.0e-12_dp, "split adjusted")

   nav = [1000.0_dp, 1100.0_dp, 1650.0_dp]
   idx = [1,3]; cf = [1000.0_dp,500.0_dp]
   call unit_prices(nav, idx, cf, up)
   call assert_close(up%price(1), 100.0_dp, 1.0e-12_dp, "unit price 1")
   call assert_close(up%price(2), 110.0_dp, 1.0e-12_dp, "unit price 2")
   call assert_close(up%price(3), 115.0_dp, 1.0e-12_dp, "unit price 3")

   print *, "test_portfolio: PASS"
contains
   subroutine assert_close(a,b,tol,label)
      real(dp), intent(in) :: a,b,tol
      character(len=*), intent(in) :: label
      if (abs(a-b) > tol) then
         print *, trim(label), a, b
         error stop 1
      end if
   end subroutine assert_close
end program test_portfolio
