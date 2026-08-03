! SPDX-License-Identifier: GPL-2.0-or-later
program tests_example
   use icsnp, only : dp, icsnp_ok, HotellingsT2, rank_ctest, ind_ctest, test_result
   implicit none
   integer, parameter :: n = 60
   real(dp) :: x(n, 4), y(n, 4)
   type(test_result) :: result
   integer :: i

   do i = 1, n
      x(i, 1) = sin(0.17_dp * real(i, dp))
      x(i, 2) = cos(0.29_dp * real(i, dp))
      x(i, 3) = sin(0.41_dp * real(i, dp))
      x(i, 4) = cos(0.13_dp * real(i, dp))
   end do
   y = x
   y(:, 1) = y(:, 1) + 0.30_dp

   call HotellingsT2(x, result, y=y)
   call print_test(result)
   call rank_ctest(x, result, y=y, scores='rank')
   call print_test(result)
   call ind_ctest(x, [1, 2], result, scores='normal')
   call print_test(result)
contains
   subroutine print_test(value)
      type(test_result), intent(in) :: value
      if (value%status /= icsnp_ok) error stop 'test failed'
      write(*, '(a)') trim(value%method)
      write(*, '(a,f12.6,a,f12.6)') '  statistic = ', value%statistic, &
         ', p-value = ', value%p_value
   end subroutine print_test
end program tests_example
