! SPDX-License-Identifier: Artistic-2.0
program test_simulation_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
   use ldhmm
   implicit none
   type(ldhmm_model) :: model, simulated, advanced
   real(dp) :: param(2,3), gamma_matrix(2,2), delta(2), x(10)
   real(dp), allocatable :: moving_average(:), acf(:), simulated_acf(:), cleaned(:)
   integer :: status

   param(1,:) = [0.002_dp,0.015_dp,1.0_dp]
   param(2,:) = [-0.003_dp,0.040_dp,1.4_dp]
   gamma_matrix(1,:) = [0.97_dp,0.03_dp]
   gamma_matrix(2,:) = [0.08_dp,0.92_dp]
   delta = [0.70_dp,0.30_dp]
   model = ldhmm_create(2,param,gamma_matrix,delta,status=status)

   x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp,8.0_dp,9.0_dp,10.0_dp]
   moving_average = ldhmm_sma(x,4)
   call assert_close(moving_average(4),2.5_dp,1.0e-14_dp,'sma first full window')
   call assert_close(moving_average(10),8.5_dp,1.0e-14_dp,'sma last window')

   x(5) = ieee_value(0.0_dp,ieee_quiet_nan)
   moving_average = ldhmm_sma(x,3,na_backfill=.true.)
   call assert_true(all(ieee_is_finite(moving_average)),'sma backfill')
   cleaned = ldhmm_drop_outliers([1.0_dp,-9.0_dp,2.0_dp,3.0_dp],1)
   call assert_true(size(cleaned) == 3,'drop outlier count')
   call assert_close(maxval(abs(cleaned)), 3.0_dp, 1.0e-15_dp, 'drop outlier value')

   acf = ldhmm_abs_acf([1.0_dp,-2.0_dp,3.0_dp,-4.0_dp,5.0_dp],2)
   call assert_true(size(acf) == 2,'acf size')
   call assert_true(all(ieee_is_finite(acf)),'acf finite')

   call seed_random(12345)
   simulated = ldhmm_simulate_state_transition(model,init=2000,status=status)
   call assert_true(status == LDHMM_SUCCESS,'simulation status')
   call assert_true(size(simulated%observations) == 2000,'simulation size')
   call assert_true(all(simulated%states_local >= 1 .and. &
      simulated%states_local <= model%m),'simulation states')
   call assert_true(all(ieee_is_finite(simulated%observations)),'simulation observations')
   advanced = ldhmm_simulate_state_transition(simulated,status=status)
   call assert_true(size(advanced%states_local) == 2000,'advanced simulation size')
   simulated_acf = ldhmm_simulate_abs_acf(model,n=1000,lag_max=3,status=status)
   call assert_true(size(simulated_acf) == 3,'simulated acf size')
   call assert_true(all(ieee_is_finite(simulated_acf)),'simulated acf finite')
   print '(a)', 'test_simulation_stats: PASS'

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//message
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: message
      call assert_true(abs(actual-expected) <= tolerance, message)
   end subroutine assert_close

end program test_simulation_stats
