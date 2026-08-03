program test_analysis
   use pmwr, only : dp, drawdown_table, nav_summary_result, attribution_result, compute_drawdowns, &
                    link_return_contributions, summarize_nav, return_attribution, link_geometric1
   implicit none
   type(drawdown_table) :: d
   type(nav_summary_result) :: s
   type(attribution_result) :: ar
   real(dp), allocatable :: linked(:)
   real(dp) :: c(2,2), r(2), ret(2,2), w(2,2), br(2,2), bw(2,2), active(2)

   call compute_drawdowns([100.0_dp,120.0_dp,90.0_dp,100.0_dp,120.0_dp,110.0_dp], d)
   if (size(d%depth) /= 2) error stop "drawdown count"
   if (d%peak(1) /= 2 .or. d%trough(1) /= 3 .or. d%recover(1) /= 5) error stop "drawdown indexes"
   call assert_close(d%depth(1), 0.25_dp, 1.0e-12_dp, "drawdown depth")
   if (d%recover(2) /= 0) error stop "unrecovered drawdown"

   c(1,:) = [0.06_dp,0.04_dp]; c(2,:) = [-0.03_dp,-0.02_dp]
   r = [0.10_dp,-0.05_dp]
   call link_return_contributions(c, r, linked, link_geometric1)
   call assert_close(linked(1), 0.027_dp, 1.0e-12_dp, "linked 1")
   call assert_close(linked(2), 0.018_dp, 1.0e-12_dp, "linked 2")

   ret(1,:) = [0.02_dp,0.01_dp]; ret(2,:) = [-0.01_dp,0.03_dp]
   br(1,:) = [0.01_dp,0.00_dp]; br(2,:) = [0.00_dp,0.02_dp]
   w = 0.5_dp; bw = 0.5_dp
   call return_attribution(ret,w,br,bw,ar)
   active = sum(w*ret,dim=2)-sum(bw*br,dim=2)
   call assert_close(sum(ar%allocation(1,:)+ar%selection(1,:)+ar%interaction(1,:)), active(1), 1.0e-12_dp, "attrib 1")

   call summarize_nav([100.0_dp,110.0_dp,99.0_dp,121.0_dp], s, periods_per_year=3.0_dp)
   call assert_close(s%total_return, 0.21_dp, 1.0e-12_dp, "nav return")
   if (s%max_drawdown <= 0.0_dp) error stop "nav drawdown"

   print *, "test_analysis: PASS"
contains
   subroutine assert_close(a,b,tol,label)
      real(dp), intent(in) :: a,b,tol
      character(len=*), intent(in) :: label
      if (abs(a-b) > tol) then
         print *, trim(label), a, b
         error stop 1
      end if
   end subroutine assert_close
end program test_analysis
