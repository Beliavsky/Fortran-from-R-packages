! SPDX-License-Identifier: GPL-2.0-only
program test_mi_scalar
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mitools, only : dp, mi_combine, mi_result, mitools_success
   implicit none
   real(dp) :: estimates(4)
   real(dp) :: variances(4)
   type(mi_result) :: result
   integer :: status

   estimates = [2.5_dp, 2.5_dp, 2.5_dp, 2.5_dp]
   variances = [0.04_dp, 0.04_dp, 0.04_dp, 0.04_dp]
   call mi_combine(estimates, variances, result, status)
   if (status /= mitools_success) error stop "scalar combine failed"
   if (result%nimp /= 4) error stop "wrong imputation count"
   if (abs(result%coefficients(1) - 2.5_dp) > 1.0e-14_dp) error stop "wrong scalar estimate"
   if (abs(result%variance(1, 1) - 0.04_dp) > 1.0e-14_dp) error stop "wrong scalar variance"
   if (ieee_is_finite(result%df(1))) error stop "zero between variance should give infinite df"
   if (abs(result%missinfo(1)) > 1.0e-14_dp) error stop "zero between variance should give zero missinfo"
   print *, "test_mi_scalar: PASS"
end program test_mi_scalar
