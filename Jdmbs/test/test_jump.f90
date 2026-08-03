! SPDX-License-Identifier: GPL-2.0-or-later
program test_jump
   use jdmbs
   implicit none

   type(jdmbs_control) :: control
   type(jdmbs_result) :: diffusion, no_jump, jump_fit
   real(dp) :: spot(2), mu(2), sigma(2), strike(2), parity(2)

   spot = [100.0_dp, 120.0_dp]
   mu = [0.03_dp, 0.01_dp]
   sigma = [0.2_dp, 0.15_dp]
   strike = [100.0_dp, 110.0_dp]
   control%day = 60
   control%monte_carlo = 800
   control%seed = 246813579_int64
   control%legacy_mode = .false.

   call normal_bs(spot, mu, sigma, strike, diffusion, control)
   call jdm_bs(spot, mu, sigma, 0.0_dp, strike, no_jump, control)
   call check(no_jump%status == jdmbs_success, 'zero-jump status')
   call check(maxval(abs(diffusion%terminal_price - no_jump%terminal_price)) < tiny(1.0_dp), &
      'zero lambda equals diffusion')

   call jdm_bs(spot, mu, sigma, 3.0_dp, strike, jump_fit, control)
   call check(jump_fit%status == jdmbs_success, 'jump status')
   call check(jump_fit%jump_events > 0, 'jump events generated')
   call check(all(jump_fit%terminal_price >= 0.0_dp), 'jump prices nonnegative')
   parity = exp(-control%discount_rate * real(control%day, dp) / &
      control%days_per_year) * &
      (sum(jump_fit%terminal_price, dim=2) / real(control%monte_carlo, dp) - strike)
   call check(maxval(abs(jump_fit%call_price - jump_fit%put_price - parity)) &
      < 2.0e-12_dp, 'pathwise put-call identity')

   print '(a)', 'test_jump: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_jump
