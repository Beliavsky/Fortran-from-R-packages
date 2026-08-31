! SPDX-License-Identifier: GPL-2.0-only
program test_mi_combine
   use mitools, only : dp, mi_combine, mi_result, mitools_success
   implicit none
   real(dp) :: estimates(2, 5)
   real(dp) :: variances(2, 2, 5)
   type(mi_result) :: result
   integer :: status

   estimates(:, 1) = [1.00_dp, 2.00_dp]
   estimates(:, 2) = [1.10_dp, 1.90_dp]
   estimates(:, 3) = [0.95_dp, 2.05_dp]
   estimates(:, 4) = [1.05_dp, 2.10_dp]
   estimates(:, 5) = [0.90_dp, 1.95_dp]

   variances(:, :, 1) = reshape([0.040_dp, 0.006_dp, 0.006_dp, 0.090_dp], [2, 2])
   variances(:, :, 2) = reshape([0.042_dp, 0.005_dp, 0.005_dp, 0.088_dp], [2, 2])
   variances(:, :, 3) = reshape([0.038_dp, 0.004_dp, 0.004_dp, 0.092_dp], [2, 2])
   variances(:, :, 4) = reshape([0.041_dp, 0.006_dp, 0.006_dp, 0.091_dp], [2, 2])
   variances(:, :, 5) = reshape([0.039_dp, 0.005_dp, 0.005_dp, 0.089_dp], [2, 2])

   call mi_combine(estimates, variances, result, status, 120.0_dp)
   if (status /= mitools_success) error stop "mi_combine returned failure"
   call assert_close(result%coefficients(1), 1.0_dp, 1.0e-13_dp, "coefficient 1")
   call assert_close(result%coefficients(2), 2.0_dp, 1.0e-13_dp, "coefficient 2")
   call assert_close(result%variance(1, 1), 0.0475_dp, 1.0e-13_dp, "variance 11")
   call assert_close(result%variance(1, 2), 0.00445_dp, 1.0e-13_dp, "variance 12")
   call assert_close(result%variance(2, 2), 0.0975_dp, 1.0e-13_dp, "variance 22")
   call assert_close(result%df(1), 62.39328123_dp, 2.0e-8_dp, "df 1")
   call assert_close(result%df(2), 94.88896372_dp, 2.0e-8_dp, "df 2")
   call assert_close(result%missinfo(1), 0.18364984_dp, 2.0e-8_dp, "missinfo 1")
   call assert_close(result%missinfo(2), 0.09578275_dp, 2.0e-8_dp, "missinfo 2")
   print *, "test_mi_combine: PASS"

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed scalar value to compare with the deterministic reference.
      real(dp), intent(in) :: expected !! Reference scalar value derived independently from the upstream formulas.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Short diagnostic label printed if the assertion fails.

      if (abs(actual - expected) > tolerance) then
         print *, trim(label), actual, expected
         error stop "assert_close failed"
      end if
   end subroutine assert_close

end program test_mi_combine
