! SPDX-License-Identifier: GPL-2.0-or-later
program test_independence_hp
   use icsnp, only : dp, icsnp_ok, ind_ctest, ind_ictest, HP_loc_test, test_result
   implicit none
   integer, parameter :: n = 72, p = 4
   real(dp) :: x(n, p), mu(p)
   integer :: i
   type(test_result) :: result

   do i = 1, n
      x(i, 1) = sin(0.17_dp * real(i, dp))
      x(i, 2) = cos(0.29_dp * real(i, dp)) + 0.15_dp * x(i, 1)
      x(i, 3) = sin(0.41_dp * real(i, dp) + 0.3_dp)
      x(i, 4) = cos(0.13_dp * real(i, dp)) - 0.10_dp * x(i, 3)
   end do
   mu = 0.0_dp

   call ind_ctest(x, [1, 2], result, scores='rank')
   call check_test(result, 'ind_ctest')
   call ind_ictest(x, [1, 2], result, scores='sign')
   call check_test(result, 'ind_ictest approximation')
   call ind_ictest(x, [1, 2], result, scores='rank', method='permutation', n_simu=20, seed=7)
   call check(result%status == icsnp_ok .and. result%replications == 20, 'ind_ictest permutation')
   call HP_loc_test(x(:, 1:3), result, mu=mu(1:3), scores='rank')
   call check_test(result, 'HP_loc_test approximation')
   call HP_loc_test(x(:, 1:3), result, mu=mu(1:3), scores='sign', method='permutation', &
      n_perm=20, seed=9)
   call check(result%status == icsnp_ok .and. result%replications == 20, 'HP_loc_test permutation')
   print '(a)', 'test_independence_hp: PASS'
contains
   subroutine check_test(value, label)
      type(test_result), intent(in) :: value
      character(len=*), intent(in) :: label
      call check(value%status == icsnp_ok .and. value%statistic >= 0.0_dp .and. &
         value%p_value >= 0.0_dp .and. value%p_value <= 1.0_dp, label)
   end subroutine check_test

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // trim(label)
         error stop 1
      end if
   end subroutine check
end program test_independence_hp
