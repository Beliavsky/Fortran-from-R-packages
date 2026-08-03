program test_returns
   use pmwr, only : dp, simple_returns, portfolio_returns_weights, portfolio_returns_positions
   implicit none
   real(dp), allocatable :: r(:), h(:,:), c(:,:)
   real(dp) :: p1(3), p2(3,2), w(3,2), pos(3,2)
   logical :: reb(3)

   p1 = [100.0_dp, 110.0_dp, 99.0_dp]
   call simple_returns(p1, r)
   call assert_close(r(1), 0.1_dp, 1.0e-12_dp, "simple return 1")
   call assert_close(r(2), -0.1_dp, 1.0e-12_dp, "simple return 2")

   p2(:,1) = [100.0_dp, 110.0_dp, 121.0_dp]
   p2(:,2) = [100.0_dp, 100.0_dp, 100.0_dp]
   w = 0.5_dp
   reb = .true.
   call portfolio_returns_weights(p2, w, reb, r, h, c)
   call assert_close(r(1), 0.05_dp, 1.0e-12_dp, "weighted return 1")
   call assert_close(r(2), 0.05_dp, 1.0e-12_dp, "weighted return 2")
   call assert_close(sum(c(1,:)), r(1), 1.0e-12_dp, "contributions 1")

   pos = 0.0_dp
   pos(:,1) = 1.0_dp
   call portfolio_returns_positions(p2, pos, reb, r)
   call assert_close(r(1), 0.1_dp, 1.0e-12_dp, "position return 1")
   call assert_close(r(2), 0.1_dp, 1.0e-12_dp, "position return 2")

   print *, "test_returns: PASS"
contains
   subroutine assert_close(x, y, tol, label)
      real(dp), intent(in) :: x, y, tol
      character(len=*), intent(in) :: label
      if (abs(x-y) > tol) then
         print *, trim(label), x, y
         error stop 1
      end if
   end subroutine assert_close
end program test_returns
