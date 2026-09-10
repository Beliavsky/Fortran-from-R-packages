! SPDX-License-Identifier: GPL-2.0-or-later
! Test/example code for the modern Fortran translation of the upstream Rssa package.
program test_gap_group
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use rssa, only : cadzow, dp, grouping_auto_wcor_ssa, grouping_result, hmatr, igapfill
   use rssa, only : gap_summary, ssa_1d, ssa_result, summarize_gaps
   implicit none
   real(dp) :: x(8), incomplete(8)
   real(dp), allocatable :: filled(:), filtered(:), heterogeneity(:, :)
   type(ssa_result) :: fit
   type(grouping_result) :: groups
   type(gap_summary) :: summary
   integer :: i, info
   x = [(2.0_dp**real(i - 1, dp), i=1, 8)]
   incomplete = x
   incomplete(4) = ieee_value(0.0_dp, ieee_quiet_nan)
   filled = igapfill(incomplete, 4, 1, maxiter=50, info=info)
   call assert_true(info == 0 .and. abs(filled(4) - x(4)) < 1.0e-4_dp, 'iterative gap fill')
   summary = summarize_gaps(incomplete, 4)
   call assert_true(summary%n_missing == 1, 'gap summary')
   filtered = cadzow(x, 4, 1, iterations=2, info=info)
   call assert_true(info == 0 .and. size(filtered) == size(x), 'Cadzow')
   fit = ssa_1d(x, 4, 2)
   groups = grouping_auto_wcor_ssa(fit)
   call assert_true(groups%info == 0, 'automatic grouping')
   heterogeneity = hmatr(fit)
   call assert_true(all(shape(heterogeneity) == [2, 2]), 'heterogeneity matrix')
   print *, 'test_gap_group: PASS'
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition required to be true.
      character(len=*), intent(in) :: message !! Failure diagnostic.
      if (.not. condition) error stop message
   end subroutine assert_true
end program test_gap_group
