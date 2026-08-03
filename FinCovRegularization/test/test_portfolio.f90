! SPDX-License-Identifier: GPL-2.0-only
program test_portfolio
   use fincovregularization
   implicit none
   real(dp) :: covariance3(3,3), covariance2(2,2), weights3(3), weights2(2), expected3(3), objective
   integer :: status

   covariance3 = 0.0_dp
   covariance3(1,1) = 1.0_dp
   covariance3(2,2) = 4.0_dp
   covariance3(3,3) = 9.0_dp
   expected3 = [36.0_dp/49.0_dp, 9.0_dp/49.0_dp, 4.0_dp/49.0_dp]

   weights3 = gmvp(covariance3, .true., status)
   call assert_true(status == fincov_ok, 'GMVP short status')
   call assert_vector_close(weights3, expected3, 1.0e-10_dp, 'GMVP short')

   weights3 = gmvp(covariance3, .false., status)
   call assert_true(status == fincov_ok, 'GMVP long-only status')
   call assert_vector_close(weights3, expected3, 1.0e-10_dp, 'GMVP long-only')

   covariance2 = reshape([1.0_dp,2.0_dp,2.0_dp,5.0_dp],[2,2])
   weights2 = gmvp(covariance2, .false., status)
   call assert_true(status == fincov_ok, 'GMVP boundary status')
   call assert_vector_close(weights2, [1.0_dp,0.0_dp], 1.0e-10_dp, 'GMVP boundary')

   covariance2 = 0.0_dp
   covariance2(1,1) = 1.0_dp
   covariance2(2,2) = 4.0_dp
   weights2 = risk_parity(covariance2, status, objective_value=objective)
   call assert_true(status == fincov_ok, 'risk parity status')
   call assert_close(sum(weights2), 1.0_dp, 1.0e-12_dp, 'risk parity sum')
   call assert_true(objective < 1.0e-12_dp, 'risk parity objective')
   call assert_vector_close(weights2, [2.0_dp/3.0_dp,1.0_dp/3.0_dp], 2.0e-5_dp, 'risk parity weights')

   print '(a)', 'test_portfolio: PASS'
contains
   subroutine assert_close(actual_value, expected_value, tolerance, label)
      real(dp), intent(in) :: actual_value, expected_value, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual_value-expected_value) > tolerance) then
         print '(a,2es24.14)', trim(label)//' failed: ', actual_value, expected_value
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_vector_close(actual_vector, expected_vector, tolerance, label)
      real(dp), intent(in) :: actual_vector(:), expected_vector(:), tolerance
      character(len=*), intent(in) :: label
      if (maxval(abs(actual_vector-expected_vector)) > tolerance) then
         print '(a,es24.14)', trim(label)//' failed, max error: ', maxval(abs(actual_vector-expected_vector))
         print '(100f12.6)', actual_vector
         error stop 1
      end if
   end subroutine assert_vector_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print '(a)', trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true
end program test_portfolio
