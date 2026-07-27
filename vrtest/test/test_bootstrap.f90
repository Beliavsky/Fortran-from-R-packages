! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
program test_bootstrap
   use vrtest, only : dp, seed_random, bootstrap_result, auto_bootstrap_result, &
      panel_vr_result, subsample_result, wright_critical_result, &
      variance_ratio_bootstrap, automatic_vr_bootstrap, panel_variance_ratio, &
      subsample_variance_ratio, wright_critical_values, joint_wright_critical_values
   implicit none
   integer, parameter :: n = 60
   real(dp) :: y(n), panel(n,3)
   integer, parameter :: kvec(3) = [2,5,10]
   integer :: i
   type(bootstrap_result) :: boot
   type(auto_bootstrap_result) :: autoboot
   type(panel_vr_result) :: pvr
   type(subsample_result) :: sub
   type(wright_critical_result) :: wc, jwc

   do i = 1, n
      y(i) = 0.01_dp*sin(0.37_dp*real(i,dp))+0.006_dp*cos(0.11_dp*real(i,dp)) + &
         0.002_dp*sin(1.31_dp*real(i,dp))
      panel(i,1) = y(i)
      panel(i,2) = 0.8_dp*y(i)+0.004_dp*cos(0.43_dp*real(i,dp))
      panel(i,3) = -0.3_dp*y(i)+0.007_dp*sin(0.71_dp*real(i,dp))
   end do

   call seed_random(20260725)
   boot = variance_ratio_bootstrap(y,kvec,40,'Mammen')
   call assert_probabilities(boot%lm_p_values,'LM bootstrap p-values')
   call assert_probability(boot%cd_p_value,'CD bootstrap p-value')
   call assert_true(all(shape(boot%confidence_intervals) == [3,2]),'bootstrap CI shape')

   call seed_random(20260725)
   autoboot = automatic_vr_bootstrap(y,40,'Rademacher')
   call assert_probability(autoboot%p_value,'automatic bootstrap p-value')
   call assert_true(size(autoboot%statistic_interval) == 2,'automatic CI shape')

   call seed_random(20260725)
   pvr = panel_variance_ratio(panel,30)
   call assert_probability(pvr%max_absolute_p_value,'panel max p-value')
   call assert_probability(pvr%sum_square_p_value,'panel sum-square p-value')
   call assert_probability(pvr%mean_p_value,'panel mean p-value')

   sub = subsample_variance_ratio(y,kvec)
   call assert_true(all(sub%block_lengths >= maxval(kvec)),'valid subsample lengths')
   call assert_probabilities(sub%p_values,'subsample p-values')

   call seed_random(99)
   wc = wright_critical_values(40,5,50)
   call assert_true(all(shape(wc%critical_values) == [6,3]),'Wright critical shape')
   call assert_true(all(wc%critical_values(1,:) <= wc%critical_values(6,:)),'ordered Wright tails')

   call seed_random(99)
   jwc = joint_wright_critical_values(40,kvec,50)
   call assert_true(all(shape(jwc%critical_values) == [3,3]),'joint critical shape')
   call assert_true(all(jwc%critical_values(1,:) <= jwc%critical_values(3,:)),'ordered joint quantiles')

   print '(a)', 'test_bootstrap: PASS'
contains
   subroutine assert_probability(value,label)
      real(dp), intent(in) :: value
      character(len=*), intent(in) :: label
      call assert_true(value >= 0.0_dp .and. value <= 1.0_dp,label)
   end subroutine assert_probability

   subroutine assert_probabilities(value,label)
      real(dp), intent(in) :: value(:)
      character(len=*), intent(in) :: label
      call assert_true(all(value >= 0.0_dp .and. value <= 1.0_dp),label)
   end subroutine assert_probabilities

   subroutine assert_true(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print '(a)', trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true
end program test_bootstrap
