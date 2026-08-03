! SPDX-License-Identifier: GPL-2.0-only
program test_regularization
   use fincovregularization
   implicit none
   real(dp) :: diagonal_matrix(2,2), sigma(3,3), expected(3,3), actual(3,3), value
   integer :: status

   diagonal_matrix = 0.0_dp
   diagonal_matrix(1,1) = 3.0_dp
   diagonal_matrix(2,2) = 4.0_dp
   call assert_close(f_norm2(diagonal_matrix), 25.0_dp, 1.0e-12_dp, 'F.norm2')
   value = o_norm2(diagonal_matrix, status)
   call assert_true(status == fincov_ok, 'O.norm2 status')
   call assert_close(value, 16.0_dp, 1.0e-10_dp, 'O.norm2')

   sigma = reshape([&
      2.0_dp, 0.4_dp, 0.1_dp, &
      0.4_dp, 3.0_dp, -0.3_dp, &
      0.1_dp, -0.3_dp, 4.0_dp], [3,3])

   expected = sigma
   expected(1,3) = 0.0_dp
   expected(3,1) = 0.0_dp
   actual = banding(sigma, 1, status)
   call assert_true(status == fincov_ok, 'banding status')
   call assert_matrix_close(actual, expected, 1.0e-12_dp, 'banding')

   actual = tapering(sigma, 2.0_dp, 0.5_dp, status)
   call assert_true(status == fincov_ok, 'tapering status')
   call assert_matrix_close(actual, expected, 1.0e-12_dp, 'tapering')

   expected = sigma
   expected(1,3) = 0.0_dp
   expected(3,1) = 0.0_dp
   actual = hard_thresholding(sigma, 0.25_dp, status)
   call assert_matrix_close(actual, expected, 1.0e-12_dp, 'hard thresholding')

   expected = 0.0_dp
   expected(1,1) = 2.0_dp
   expected(2,2) = 3.0_dp
   expected(3,3) = 4.0_dp
   expected(1,2) = 0.15_dp
   expected(2,1) = 0.15_dp
   expected(2,3) = -0.05_dp
   expected(3,2) = -0.05_dp
   actual = soft_thresholding(sigma, 0.25_dp, status)
   call assert_matrix_close(actual, expected, 1.0e-12_dp, 'soft thresholding')

   diagonal_matrix = reshape([1.0_dp,1.2_dp,1.2_dp,1.0_dp],[2,2])
   value = threshold_min(diagonal_matrix, 'soft', status=status)
   call assert_true(status == fincov_ok, 'threshold_min status')
   call assert_close(value, 0.2_dp, 2.0e-7_dp, 'threshold_min soft')

   actual = ind_cov(sigma, status)
   call assert_true(status == fincov_ok, 'Ind.Cov status')
   call assert_true(maxval(abs(actual - transpose(actual))) < 1.0e-12_dp, 'Ind.Cov symmetry')
   call assert_true(maxval(abs(actual - reshape([actual(1,1),0.0_dp,0.0_dp, &
      0.0_dp,actual(2,2),0.0_dp,0.0_dp,0.0_dp,actual(3,3)],[3,3]))) <= tiny(1.0_dp), 'Ind.Cov diagonal')

   print '(a)', 'test_regularization: PASS'
contains
   subroutine assert_close(actual_value, expected_value, tolerance, label)
      real(dp), intent(in) :: actual_value, expected_value, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual_value-expected_value) > tolerance) then
         print '(a,2es24.14)', trim(label)//' failed: ', actual_value, expected_value
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_matrix_close(actual_matrix, expected_matrix, tolerance, label)
      real(dp), intent(in) :: actual_matrix(:,:), expected_matrix(:,:), tolerance
      character(len=*), intent(in) :: label
      if (maxval(abs(actual_matrix-expected_matrix)) > tolerance) then
         print '(a,es24.14)', trim(label)//' failed, max error: ', maxval(abs(actual_matrix-expected_matrix))
         error stop 1
      end if
   end subroutine assert_matrix_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print '(a)', trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true
end program test_regularization
