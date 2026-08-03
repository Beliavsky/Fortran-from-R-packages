! SPDX-License-Identifier: GPL-2.0-or-later
program test_normal
   use jdmbs
   implicit none

   type(jdmbs_control) :: control
   type(jdmbs_result) :: fit, repeat_fit
   real(dp) :: spot(2), mu(2), sigma(2), strike(2), expected(2), discount

   spot = [100.0_dp, 80.0_dp]
   mu = [0.05_dp, -0.02_dp]
   sigma = 0.0_dp
   strike = [95.0_dp, 82.0_dp]
   control%day = 73
   control%monte_carlo = 5
   control%seed = 987654321_int64
   control%discount_rate = 0.03_dp
   control%legacy_mode = .false.
   control%store_paths = .true.

   call normal_bs(spot, mu, sigma, strike, fit, control)
   expected = spot * exp(mu * real(control%day, dp) / control%days_per_year)
   discount = exp(-control%discount_rate * real(control%day, dp) / &
      control%days_per_year)
   call check(fit%status == jdmbs_success, 'standard status')
   call check(maxval(abs(fit%terminal_price(:, 1) - expected)) < 1.0e-12_dp, &
      'standard deterministic terminal')
   call check(maxval(abs(fit%paths(:, 1, 0) - spot)) < 1.0e-14_dp, &
      'standard path starts at spot')
   call check(maxval(abs(fit%call_price - discount * max(expected - strike, 0.0_dp))) &
      < 1.0e-12_dp, 'discounted call')
   call check(maxval(abs(fit%put_price - discount * max(strike - expected, 0.0_dp))) &
      < 1.0e-12_dp, 'discounted put')
   call check(maxval(abs(fit%call_se)) < tiny(1.0_dp), 'zero deterministic standard error')

   control%legacy_mode = .true.
   call normal_bs(spot, mu, sigma, strike, fit, control)
   expected = spot * exp(mu * real(control%day - 1, dp) / control%days_per_year)
   call check(maxval(abs(fit%terminal_price(:, 1) - expected)) < 1.0e-12_dp, &
      'legacy off-by-one terminal')

   sigma = [0.2_dp, 0.3_dp]
   call normal_bs(spot, mu, sigma, strike, fit, control)
   call normal_bs(spot, mu, sigma, strike, repeat_fit, control)
   call check(maxval(abs(fit%terminal_price - repeat_fit%terminal_price)) < tiny(1.0_dp), &
      'fixed seed reproducibility')
   call check(all(fit%terminal_price >= 0.0_dp), 'nonnegative terminal prices')

   print '(a)', 'test_normal: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_normal
