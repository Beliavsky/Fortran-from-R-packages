! SPDX-License-Identifier: GPL-2.0-only
program test_imputation
   use mitools, only : dp, imputation_cbind, imputation_dimensions, imputation_get
   use mitools, only : imputation_list, imputation_list_from_array, imputation_rbind, mitools_success
   implicit none
   real(dp) :: a(2, 2, 2)
   real(dp) :: b(1, 2, 2)
   real(dp) :: c(2, 1, 2)
   real(dp), allocatable :: dataset(:, :)
   type(imputation_list) :: col_bound
   type(imputation_list) :: first
   type(imputation_list) :: rows
   type(imputation_list) :: row_bound
   type(imputation_list) :: second
   integer :: ncol
   integer :: nimp
   integer :: nrow
   integer :: status

   a(:, :, 1) = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2])
   a(:, :, 2) = reshape([11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp], [2, 2])
   b(:, :, 1) = reshape([5.0_dp, 6.0_dp], [1, 2])
   b(:, :, 2) = reshape([15.0_dp, 16.0_dp], [1, 2])
   c(:, :, 1) = reshape([7.0_dp, 8.0_dp], [2, 1])
   c(:, :, 2) = reshape([17.0_dp, 18.0_dp], [2, 1])

   call imputation_list_from_array(a, first, status)
   if (status /= mitools_success) error stop "first constructor failed"
   call imputation_list_from_array(b, rows, status)
   if (status /= mitools_success) error stop "row constructor failed"
   call imputation_rbind(first, rows, row_bound, status)
   if (status /= mitools_success) error stop "rbind failed"
   call imputation_dimensions(row_bound, nrow, ncol, nimp, status)
   if (status /= mitools_success .or. nrow /= 3 .or. ncol /= 2 .or. nimp /= 2) error stop "rbind dimensions"
   call imputation_get(row_bound, 2, dataset, status)
   if (status /= mitools_success) error stop "imputation_get failed"
   if (abs(dataset(3, 2) - 16.0_dp) > 1.0e-14_dp) error stop "rbind value mismatch"

   call imputation_list_from_array(c, second, status)
   if (status /= mitools_success) error stop "second constructor failed"
   call imputation_cbind(first, second, col_bound, status)
   if (status /= mitools_success) error stop "cbind failed"
   if (abs(col_bound%values(2, 3, 2) - 18.0_dp) > 1.0e-14_dp) error stop "cbind value mismatch"
   print *, "test_imputation: PASS"
end program test_imputation
