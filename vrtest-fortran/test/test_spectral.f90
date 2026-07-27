! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
program test_spectral
   use vrtest, only : dp, seed_random, auto_q_result, average_exponential_result, &
      spectral_shape_result, generalized_spectral_result, dl_result, chen_deo_result, &
      automatic_portmanteau, average_exponential_test, spectral_shape_test, &
      generalized_spectral_test, dominguez_lobato_statistic, dominguez_lobato_test, &
      chen_deo_test
   implicit none
   integer, parameter :: n = 80
   integer, parameter :: kvec(3) = [2,5,10]
   real(dp) :: y(n)
   integer :: i
   type(auto_q_result) :: aq
   type(average_exponential_result) :: ae
   type(spectral_shape_result) :: ss
   type(generalized_spectral_result) :: gs
   type(dl_result) :: dl0, dl
   type(chen_deo_result) :: cd

   do i = 1, n
      y(i) = 0.013_dp*sin(0.41_dp*real(i,dp))+0.009_dp*cos(0.17_dp*real(i,dp)) + &
         0.004_dp*sin(1.07_dp*real(i,dp))+0.00003_dp*real(i-40,dp)
   end do

   aq = automatic_portmanteau(y,10)
   call assert_true(aq%selected_lag >= 1 .and. aq%selected_lag <= 10,'Auto.Q selected lag')
   call assert_probability(aq%p_value,'Auto.Q p-value')

   ae = average_exponential_test(y)
   call assert_true(ae%exponential_lm >= 0.0_dp,'average exponential LM')
   call assert_true(ae%exponential_lr >= 0.0_dp,'average exponential LR')

   ss = spectral_shape_test(y)
   call assert_true(ss%anderson_darling >= 0.0_dp,'spectral AD')
   call assert_true(ss%cramer_von_mises >= 0.0_dp,'spectral CvM')
   call assert_true(ss%maximum >= 0.0_dp,'spectral M')

   call seed_random(777)
   gs = generalized_spectral_test(y,20)
   call assert_true(gs%statistic >= 0.0_dp,'generalized spectral statistic')
   call assert_probability(gs%p_value,'generalized spectral p-value')
   call assert_true(all(gs%bootstrap_critical_values(1:2) <= gs%bootstrap_critical_values(2:3)), &
      'generalized spectral critical values')

   dl0 = dominguez_lobato_statistic(y,2)
   call assert_true(dl0%cp_statistic >= 0.0_dp,'DL Cp statistic')
   call assert_true(dl0%kp_statistic >= 0.0_dp,'DL Kp statistic')
   call seed_random(777)
   dl = dominguez_lobato_test(y,20,2)
   call assert_probability(dl%cp_p_value,'DL Cp p-value')
   call assert_probability(dl%kp_p_value,'DL Kp p-value')

   cd = chen_deo_test(y,kvec)
   call assert_true(cd%solve_info == 0,'Chen-Deo covariance solve')
   call assert_true(cd%qp_statistic >= 0.0_dp,'Chen-Deo QP')
   call assert_true(all(cd%chi_square_upper_quantiles > 0.0_dp),'Chen-Deo quantiles')

   print '(a)', 'test_spectral: PASS'
contains
   subroutine assert_probability(value,label)
      real(dp), intent(in) :: value
      character(len=*), intent(in) :: label
      call assert_true(value >= 0.0_dp .and. value <= 1.0_dp,label)
   end subroutine assert_probability

   subroutine assert_true(condition,label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         print '(a)', trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true
end program test_spectral
