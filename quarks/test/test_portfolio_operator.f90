program test_portfolio_operator
   use quarks
   implicit none
   real(dp), parameter :: tol = 2.0e-14_dp
   real(dp) :: x(2,2), weights(2), time_weights(2,2)
   type(pl_result) :: linear, exact_vector, exact_matrix, corrected

   x(1,:) = [0.01_dp, 0.02_dp]
   x(2,:) = [-0.01_dp, 0.03_dp]
   weights = [0.4_dp, 0.6_dp]
   time_weights(1,:) = weights
   time_weights(2,:) = [0.7_dp, 0.2_dp]

   linear = plop(x, weights, 1)
   call assert_close(linear%pl(1), 0.016_dp, tol, 'linear row 1')
   call assert_close(linear%pl(2), 0.014_dp, tol, 'linear row 2')

   exact_vector = plop(x, weights, 0)
   call assert_close(exact_vector%pl(1), &
      0.4_dp * exp(0.01_dp) + 0.6_dp * exp(0.02_dp) - 1.0_dp, &
      tol, 'exact vector')

   exact_matrix = plop_time_varying(x, time_weights, 0, .true.)
   call assert_close(exact_matrix%pl(2), &
      sum(exp(x(2,:)) * time_weights(2,:) - 1.0_dp), tol, 'upstream matrix')

   corrected = plop_time_varying(x, time_weights, 0, .false.)
   call assert_close(corrected%pl(2), &
      sum((exp(x(2,:)) - 1.0_dp) * time_weights(2,:)), tol, 'corrected matrix')
   print *, 'test_portfolio_operator: PASS'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print *, trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_portfolio_operator
