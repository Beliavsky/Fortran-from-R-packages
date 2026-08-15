program test_mle
   use countdm
   implicit none
   integer :: failures
   integer, allocatable :: x(:), xb(:)
   real(dp) :: meanx, alpha_exact
   type(bell_closed_result_t) :: bc
   type(mle_result_t) :: fit

   failures = 0
   call get_data_sbirth(x)
   meanx = real(sum(x), dp) / real(size(x), dp)

   bc = bell_mle_closed(x)
   call check_close('closed Bell theta', bc%theta, 0.3170451070377209_dp, 2.0e-12_dp)
   call check_true('closed Bell loglik finite', bc%loglik > -huge(1.0_dp) / 10.0_dp)

   fit = mle_poisson(x, 0.5_dp)
   call check_true('Poisson converged', fit%converged)
   call check_close('Poisson MLE equals mean', fit%estimate(1), meanx, 3.0e-6_dp)
   call check_true('Poisson SE positive', fit%se(1) > 0.0_dp)

   fit = mle_bell(x, 0.3_dp)
   call check_true('Bell converged', fit%converged)
   call check_close('numeric Bell matches closed', fit%estimate(1), bc%theta, 3.0e-6_dp)

   allocate(xb(10))
   xb = [1, 1, 1, 2, 2, 2, 3, 3, 4, 5]
   alpha_exact = real(sum(xb) - size(xb), dp) / real(sum(xb), dp)
   fit = mle_borel(xb, 0.4_dp)
   call check_true('Borel converged', fit%converged)
   call check_close('Borel MLE', fit%estimate(1), alpha_exact, 5.0e-6_dp)

   fit = mle_bt(x, 0.3_dp, 1.0_dp)
   call check_true('Bell-Touchard converged', fit%converged)
   call check_true('Bell-Touchard parameter domain', all(fit%estimate > 0.0_dp))
   call check_true('Bell-Touchard finite AIC', fit%aic < huge(1.0_dp) / 10.0_dp)
   call check_close('Bell-Touchard lambda oracle', fit%estimate(1), 1.1924964407_dp, 2.0e-5_dp)
   call check_close('Bell-Touchard theta oracle', fit%estimate(2), 0.1107797240_dp, 2.0e-5_dp)

   fit = mle_zip(x, 0.3_dp, 0.8_dp)
   call check_true('ZIP converged', fit%converged)
   call check_true('ZIP parameter domain', fit%estimate(1) > 0.0_dp .and. fit%estimate(1) < 1.0_dp &
      .and. fit%estimate(2) > 0.0_dp)
   call check_close('ZIP alpha oracle', fit%estimate(1), 0.7241909474_dp, 3.0e-5_dp)
   call check_close('ZIP theta oracle', fit%estimate(2), 1.5783506518_dp, 3.0e-5_dp)

   fit = mle_zibell(x, 0.2_dp, 0.8_dp)
   call check_true('ZIBell converged', fit%converged)
   call check_true('ZIBell parameter domain', fit%estimate(1) > 0.0_dp .and. fit%estimate(1) < 1.0_dp &
      .and. fit%estimate(2) > 0.0_dp)
   call check_close('ZIBell alpha oracle', fit%estimate(1), 0.6179014377_dp, 3.0e-5_dp)
   call check_close('ZIBell lambda oracle', fit%estimate(2), 0.6155879719_dp, 3.0e-5_dp)

   fit = mle_zibellt(x, 0.5_dp, 1.0_dp, 0.2_dp)
   call check_true('ZIBT converged', fit%converged)
   call check_true('ZIBT parameter order/domain', fit%estimate(1) > 0.0_dp .and. fit%estimate(2) > 0.0_dp &
      .and. fit%estimate(3) > 0.0_dp .and. fit%estimate(3) < 1.0_dp)
   call check_close('ZIBT lambda oracle', fit%estimate(1), 0.7627038500_dp, 5.0e-5_dp)
   call check_close('ZIBT theta oracle', fit%estimate(2), 0.6129463523_dp, 5.0e-5_dp)
   call check_close('ZIBT pi oracle', fit%estimate(3), 0.5656943003_dp, 5.0e-5_dp)

   fit = mle_zoip(x, 0.2_dp, 0.1_dp, 0.8_dp)
   call check_true('ZOIP converged', fit%converged)
   call check_true('ZOIP simplex', fit%estimate(1) > 0.0_dp .and. fit%estimate(2) > 0.0_dp &
      .and. sum(fit%estimate(1:2)) < 1.0_dp .and. fit%estimate(3) > 0.0_dp)
   call check_close('ZOIP alpha oracle', fit%estimate(1), 0.7685519980_dp, 6.0e-5_dp)
   call check_close('ZOIP beta oracle', fit%estimate(2), 0.0889193932_dp, 6.0e-5_dp)
   call check_close('ZOIP theta oracle', fit%estimate(3), 2.4304171827_dp, 6.0e-5_dp)

   fit = mle_zoibell(x, 0.2_dp, 0.1_dp, 0.8_dp)
   call check_true('ZOIBell converged', fit%converged)
   call check_true('ZOIBell simplex', fit%estimate(1) > 0.0_dp .and. fit%estimate(2) > 0.0_dp &
      .and. sum(fit%estimate(1:2)) < 1.0_dp .and. fit%estimate(3) > 0.0_dp)
   call check_close('ZOIBell alpha oracle', fit%estimate(1), 0.7810945259_dp, 8.0e-5_dp)
   call check_close('ZOIBell beta oracle', fit%estimate(2), 0.1194029812_dp, 8.0e-5_dp)
   call check_close('ZOIBell theta oracle', fit%estimate(3), 0.7627195388_dp, 8.0e-5_dp)

   if (failures == 0) then
      print '(a)', 'test_mle: PASS'
   else
      print '(a,i0)', 'test_mle: FAIL ', failures
      error stop 1
   end if

contains
   subroutine check_close(name, got, expected, tol)
      character(*), intent(in) :: name
      real(dp), intent(in) :: got, expected, tol
      if (abs(got - expected) > tol * max(1.0_dp, abs(expected))) then
         failures = failures + 1
         print '(a,2(1x,es24.16))', 'failure: '//trim(name), got, expected
      end if
   end subroutine check_close

   subroutine check_true(name, condition)
      character(*), intent(in) :: name
      logical, intent(in) :: condition
      if (.not. condition) then
         failures = failures + 1
         print '(a)', 'failure: '//trim(name)
      end if
   end subroutine check_true
end program test_mle
