program test_trades
   use pmwr, only : dp, journal_type, pl_summary_result, pl_path_result, make_journal, positions_at, &
                    pl_summary, pl_path, split_trade_runs, time_weighted_exposure
   implicit none
   type(journal_type) :: j
   type(pl_summary_result) :: s
   type(pl_path_result) :: path
   real(dp), allocatable :: pos(:,:), a(:), p(:), t(:)
   integer, allocatable :: first(:), last(:)
   real(dp) :: exposure

   call make_journal([1.0_dp,2.0_dp,3.0_dp], [10.0_dp,-4.0_dp,-6.0_dp], &
                     [100.0_dp,110.0_dp,120.0_dp], [1,1,1], j)
   call positions_at(j, [1.0_dp,2.0_dp,3.0_dp], 1, pos)
   call assert_close(pos(1,1), 10.0_dp, 1.0e-12_dp, "position 1")
   call assert_close(pos(2,1), 6.0_dp, 1.0e-12_dp, "position 2")
   call assert_close(pos(3,1), 0.0_dp, 1.0e-12_dp, "position 3")

   call pl_summary(j%amount, j%price, s)
   if (.not. s%closed) error stop "trade should be closed"
   call assert_close(s%total_pl, 160.0_dp, 1.0e-12_dp, "total pl")
   call assert_close(s%average_buy, 100.0_dp, 1.0e-12_dp, "average buy")
   call assert_close(s%average_sell, 116.0_dp, 1.0e-12_dp, "average sell")

   call pl_path(j%amount, j%price, path)
   call assert_close(path%realized(2), 40.0_dp, 1.0e-12_dp, "realized 2")
   call assert_close(path%realized(3), 160.0_dp, 1.0e-12_dp, "realized 3")

   call split_trade_runs([10.0_dp,-15.0_dp], [100.0_dp,110.0_dp], [1.0_dp,2.0_dp], a,p,t,first,last)
   if (size(a) /= 3 .or. size(first) /= 2) error stop "split trade dimensions"
   call assert_close(a(2), -10.0_dp, 1.0e-12_dp, "split close")
   call assert_close(a(3), -5.0_dp, 1.0e-12_dp, "split open")

   exposure = time_weighted_exposure([10.0_dp,-10.0_dp], [1.0_dp,3.0_dp], start_time=0.0_dp, end_time=4.0_dp)
   call assert_close(exposure, 5.0_dp, 1.0e-12_dp, "exposure")

   print *, "test_trades: PASS"
contains
   subroutine assert_close(a,b,tol,label)
      real(dp), intent(in) :: a,b,tol
      character(len=*), intent(in) :: label
      if (abs(a-b) > tol) then
         print *, trim(label), a, b
         error stop 1
      end if
   end subroutine assert_close
end program test_trades
