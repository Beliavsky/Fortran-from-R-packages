! SPDX-License-Identifier: GPL-2.0-only
program test_errors
   use mitools, only : dp, imputation_get, imputation_list, imputation_list_from_array
   use mitools, only : mi_combine, mi_confidence_intervals, mi_result
   use mitools, only : mitools_invalid_index, mitools_invalid_probability, mitools_success
   implicit none
   real(dp) :: cube(1, 1, 2)
   real(dp) :: estimates(2)
   real(dp) :: variances(2)
   real(dp), allocatable :: dataset(:, :)
   real(dp), allocatable :: lower(:)
   real(dp), allocatable :: upper(:)
   type(imputation_list) :: list
   type(mi_result) :: result
   integer :: status

   estimates = [1.0_dp, 1.1_dp]
   variances = [0.1_dp, 0.1_dp]
   call mi_combine(estimates, variances, result, status)
   if (status /= mitools_success) error stop "setup combine failed"
   call mi_confidence_intervals(result, 0.0_dp, lower, upper, status)
   if (status /= mitools_invalid_probability) error stop "invalid alpha status mismatch"

   cube(1, 1, 1) = 1.0_dp
   cube(1, 1, 2) = 2.0_dp
   call imputation_list_from_array(cube, list, status)
   if (status /= mitools_success) error stop "list setup failed"
   call imputation_get(list, 3, dataset, status)
   if (status /= mitools_invalid_index) error stop "invalid index status mismatch"
   print *, "test_errors: PASS"
end program test_errors
