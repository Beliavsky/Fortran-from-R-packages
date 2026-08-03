! SPDX-License-Identifier: GPL-2.0-or-later
program test_correlated_jump
   use jdmbs
   implicit none

   type(jdmbs_control) :: control
   type(jdmbs_result) :: fit, zero_transmission
   real(dp) :: transmission(3, 3), zero_matrix(3, 3)
   real(dp) :: spot(3), mu(3), sigma(3), strike(3), event_mean

   spot = 100.0_dp
   mu = 0.0_dp
   sigma = 0.0_dp
   strike = 100.0_dp
   transmission = 1.0_dp
   zero_matrix = 0.0_dp
   control%day = 90
   control%monte_carlo = 4000
   control%seed = 135792468_int64
   control%legacy_mode = .false.
   control%store_paths = .true.

   call jdm_new_bs(transmission, spot, mu, sigma, 2.5_dp, strike, fit, control)
   call check(fit%status == jdmbs_success, 'correlated status')
   call check(maxval(abs(fit%terminal_price(1, :) - fit%terminal_price(2, :))) &
      < 1.0e-12_dp, 'all-one transmission asset 1 and 2')
   call check(maxval(abs(fit%terminal_price(1, :) - fit%terminal_price(3, :))) &
      < 1.0e-12_dp, 'all-one transmission asset 1 and 3')
   event_mean = real(fit%jump_events, dp) / real(control%monte_carlo, dp)
   call check(abs(event_mean - 2.5_dp) < 0.12_dp, 'mean event count')
   call check(maxval(abs(fit%paths(:, :, 0) - spread(spot, 2, control%monte_carlo))) &
      < 1.0e-14_dp, 'stored initial values')

   call jdm_new_bs(zero_matrix, spot, mu, sigma, 4.0_dp, strike, &
      zero_transmission, control)
   call check(maxval(abs(zero_transmission%terminal_price - 100.0_dp)) &
      < 1.0e-12_dp, 'zero transmission removes jump effect')

   print '(a)', 'test_correlated_jump: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_correlated_jump
