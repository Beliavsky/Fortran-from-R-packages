! SPDX-License-Identifier: GPL-2.0-only
program test_pv
   use mitools, only : dp, mitools_success, pv_materialize, pv_select
   implicit none
   real(dp) :: base(3, 2)
   real(dp) :: pv(3, 2, 3)
   real(dp), allocatable :: data(:, :)
   real(dp), allocatable :: selected(:, :)
   integer :: status

   base = reshape([1.0_dp, 2.0_dp, 3.0_dp, 10.0_dp, 20.0_dp, 30.0_dp], [3, 2])
   pv(:, 1, 1) = [101.0_dp, 102.0_dp, 103.0_dp]
   pv(:, 1, 2) = [111.0_dp, 112.0_dp, 113.0_dp]
   pv(:, 1, 3) = [121.0_dp, 122.0_dp, 123.0_dp]
   pv(:, 2, 1) = [201.0_dp, 202.0_dp, 203.0_dp]
   pv(:, 2, 2) = [211.0_dp, 212.0_dp, 213.0_dp]
   pv(:, 2, 3) = [221.0_dp, 222.0_dp, 223.0_dp]

   call pv_select(pv, 2, selected, status)
   if (status /= mitools_success) error stop "pv_select failed"
   if (abs(selected(3, 2) - 213.0_dp) > 1.0e-14_dp) error stop "pv_select value mismatch"
   call pv_materialize(base, pv, 3, data, status)
   if (status /= mitools_success) error stop "pv_materialize failed"
   if (size(data, 1) /= 3 .or. size(data, 2) /= 4) error stop "pv_materialize shape mismatch"
   if (abs(data(2, 1) - 2.0_dp) > 1.0e-14_dp) error stop "base data changed"
   if (abs(data(2, 3) - 122.0_dp) > 1.0e-14_dp) error stop "first PV mismatch"
   if (abs(data(2, 4) - 222.0_dp) > 1.0e-14_dp) error stop "second PV mismatch"
   print *, "test_pv: PASS"
end program test_pv
