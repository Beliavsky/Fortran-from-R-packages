! SPDX-License-Identifier: GPL-3.0-only
program test_optimizer
   use mcrp_module
   implicit none

   real(dp) :: r(18, 4), start(4), lower(4), lambda(3)
   real(dp) :: before, after, scaled
   type(mcrp_result) :: result, result_all
   logical :: variance_only(3)
   integer :: i, failures

   failures = 0
   do i = 1, 18
      r(i, 1) = 0.006_dp * sin(0.61_dp * real(i, dp)) + &
         0.001_dp * real(mod(i, 3), dp)
      r(i, 2) = 0.004_dp * cos(0.47_dp * real(i, dp)) + &
         0.25_dp * r(i, 1) + 0.0003_dp * real(i, dp)
      r(i, 3) = 0.003_dp * sin(0.83_dp * real(i, dp)) - &
         0.20_dp * r(i, 1) + 0.00015_dp * real(i * i, dp)
      r(i, 4) = 0.005_dp * cos(0.29_dp * real(i, dp)) + &
         0.15_dp * r(i, 2) - 0.0002_dp * real(i, dp)
   end do

   start = [0.40_dp, 0.30_dp, 0.20_dp, 0.10_dp]
   lower = 0.0_dp
   lambda = 1.0_dp
   variance_only = [.true., .false., .false.]

   before = mcrp_objective_value(r, start, lambda=lambda, active=variance_only)
   call close(before, 1.4577912183570263e-2_dp, 1.0e-13_dp, &
      'fixed variance objective reference')
   call close(mcrp_objective_value(r, start, lambda=lambda), &
      2.4575340794956735_dp, 1.0e-11_dp, &
      'fixed combined objective reference')
   scaled = mcrp_objective_value(r, 7.0_dp * start, lambda=lambda, &
      active=variance_only)
   call close(before, scaled, 1.0e-14_dp, 'objective scale invariance')

   call mcrp(start, r, result, lambda=lambda, active=variance_only, &
      lower=lower, max_iterations=5000, tolerance=1.0e-11_dp)
   call check(result%status == mcrp_success, 'variance optimizer convergence')
   call close(sum(abs(result%weights)), 1.0_dp, 1.0e-12_dp, &
      'normalized weights')
   call check(all(result%weights >= -1.0e-14_dp), 'long-only weights')
   after = mcrp_objective_value(r, result%raw_parameters, lambda=lambda, &
      active=variance_only)
   call check(after < before * 1.0e-4_dp, 'variance objective reduction')
   call check(maxval(result%variance_contributions) - &
      minval(result%variance_contributions) < 2.0e-5_dp, &
      'equal variance contributions')

   before = mcrp_objective_value(r, start, lambda=lambda)
   call mcrp(start, r, result_all, lambda=lambda, lower=lower, &
      max_iterations=7000, tolerance=1.0e-10_dp)
   after = mcrp_objective_value(r, result_all%raw_parameters, lambda=lambda)
   call check(ieee_finite(after), 'all-moment finite objective')
   call check(after < before, 'all-moment objective reduction')
   call close(sum(abs(result_all%weights)), 1.0_dp, 1.0e-12_dp, &
      'all-moment normalized weights')
   call close(sum(result_all%variance_contributions), 1.0_dp, 1.0e-8_dp, &
      'variance contribution identity')
   call close(sum(result_all%skewness_contributions), 1.0_dp, 1.0e-7_dp, &
      'skew contribution identity')
   call close(sum(result_all%kurtosis_contributions), 1.0_dp, 1.0e-7_dp, &
      'kurt contribution identity')

   if (failures /= 0) then
      write(*, '(a, i0)') 'test_optimizer failures: ', failures
      error stop 1
   end if
   write(*, '(a)') 'test_optimizer: all tests passed'

contains

   logical function ieee_finite(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x
      ieee_finite = ieee_is_finite(x)
   end function ieee_finite

   subroutine close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         failures = failures + 1
         write(*, '(a, 2es24.15)') trim(label)//': ', actual, expected
      end if
   end subroutine close

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') trim(label)//': failed'
      end if
   end subroutine check

end program test_optimizer
