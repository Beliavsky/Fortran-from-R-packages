! SPDX-License-Identifier: GPL-2.0-only
program test_mi_summary
   use mitools, only : dp, mi_combine, mi_result, mi_summary, mitools_success
   implicit none
   real(dp) :: estimates_in(2, 5)
   real(dp) :: variances(2, 2, 5)
   real(dp), allocatable :: estimates(:)
   real(dp), allocatable :: lower(:)
   real(dp), allocatable :: missinfo(:)
   real(dp), allocatable :: se(:)
   real(dp), allocatable :: upper(:)
   type(mi_result) :: result
   integer :: status

   estimates_in(:, 1) = [1.00_dp, 2.00_dp]
   estimates_in(:, 2) = [1.10_dp, 1.90_dp]
   estimates_in(:, 3) = [0.95_dp, 2.05_dp]
   estimates_in(:, 4) = [1.05_dp, 2.10_dp]
   estimates_in(:, 5) = [0.90_dp, 1.95_dp]
   variances(:, :, 1) = reshape([0.040_dp, 0.006_dp, 0.006_dp, 0.090_dp], [2, 2])
   variances(:, :, 2) = reshape([0.042_dp, 0.005_dp, 0.005_dp, 0.088_dp], [2, 2])
   variances(:, :, 3) = reshape([0.038_dp, 0.004_dp, 0.004_dp, 0.092_dp], [2, 2])
   variances(:, :, 4) = reshape([0.041_dp, 0.006_dp, 0.006_dp, 0.091_dp], [2, 2])
   variances(:, :, 5) = reshape([0.039_dp, 0.005_dp, 0.005_dp, 0.089_dp], [2, 2])

   call mi_combine(estimates_in, variances, result, status, 120.0_dp)
   if (status /= mitools_success) error stop "combine failed"
   call mi_summary(result, 0.05_dp, estimates, se, lower, upper, missinfo, status)
   if (status /= mitools_success) error stop "summary failed"
   call assert_close(se(1), 0.21794495_dp, 2.0e-8_dp, "se 1")
   call assert_close(se(2), 0.31224990_dp, 2.0e-8_dp, "se 2")
   call assert_close(lower(1), 0.56438888_dp, 3.0e-8_dp, "lower 1")
   call assert_close(upper(1), 1.43561112_dp, 3.0e-8_dp, "upper 1")
   call assert_close(lower(2), 1.38009622_dp, 3.0e-8_dp, "lower 2")
   call assert_close(upper(2), 2.61990378_dp, 3.0e-8_dp, "upper 2")
   print *, "test_mi_summary: PASS"

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed scalar value to compare with the deterministic reference.
      real(dp), intent(in) :: expected !! Independently calculated reference value.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference.
      character(len=*), intent(in) :: label !! Short diagnostic label printed if the assertion fails.

      if (abs(actual - expected) > tolerance) then
         print *, trim(label), actual, expected
         error stop "assert_close failed"
      end if
   end subroutine assert_close

end program test_mi_summary
