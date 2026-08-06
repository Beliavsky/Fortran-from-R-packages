program test_projection
   use sparse_index_tracking, only : dp, project_capped_simplex, bisection
   implicit none

   real(dp), parameter :: tolerance = 2.0e-12_dp
   real(dp) :: z(4), expected(4)
   real(dp), allocatable :: weights(:), weights_from_c(:)
   integer :: info

   z = [0.8_dp, 0.4_dp, 0.1_dp, -0.2_dp]
   expected = [0.45_dp, 0.425_dp, 0.125_dp, 0.0_dp]

   call project_capped_simplex(z, 0.45_dp, weights, info)
   call assert_true(info == 0, 'projection returned a nonzero status')
   call assert_close(maxval(abs(weights - expected)), 0.0_dp, tolerance, 'projection fixture')
   call assert_close(sum(weights), 1.0_dp, tolerance, 'projection sum')
   call assert_true(all(weights >= 0.0_dp) .and. all(weights <= 0.45_dp), 'projection bounds')

   call bisection(-2.0_dp * z, 0.45_dp, weights_from_c, info)
   call assert_true(info == 0, 'bisection returned a nonzero status')
   call assert_close(maxval(abs(weights_from_c - expected)), 0.0_dp, tolerance, 'bisection fixture')

   call project_capped_simplex(z, 0.20_dp, weights, info)
   call assert_true(info /= 0, 'infeasible upper bound was accepted')

   print '(a)', 'test_projection: PASS'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print '(a)', 'FAIL: ' // trim(message)
         error stop 1
      end if
   end subroutine assert_true


   subroutine assert_close(actual, expected_value, tol, message)
      real(dp), intent(in) :: actual, expected_value, tol
      character(len=*), intent(in) :: message

      call assert_true(abs(actual - expected_value) <= tol, message)
   end subroutine assert_close

end program test_projection
