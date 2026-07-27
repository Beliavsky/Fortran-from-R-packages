! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
program test_hedging
   use opthedging, only : dp, hedging_iid, hedging_result
   implicit none

   integer :: i
   integer :: j
   real(dp) :: c2
   real(dp) :: expected_a
   real(dp) :: expected_c
   real(dp) :: m1
   real(dp) :: m2
   real(dp) :: payoff
   real(dp) :: rho
   real(dp) :: xi(4)
   real(dp) :: log_returns(4)
   type(hedging_result) :: result

   xi = [-0.20_dp, -0.05_dp, 0.05_dp, 0.20_dp]
   log_returns = log(1.0_dp + xi)
   result = hedging_iid(log_returns, 1.0_dp, 100.0_dp, 0.0_dp, .true., &
      1, 3, 80.0_dp, 120.0_dp)
   if (.not. result%ok) then
      print '(a)', result%message
      error stop 1
   end if

   m1 = sum(xi) / real(size(xi), dp)
   m2 = dot_product(xi, xi) / real(size(xi), dp)
   rho = m1 / m2
   c2 = 1.0_dp - m1 * rho
   call check_close(result%rho, rho, 1.0e-15_dp, "rho")

   do i = 1, size(result%s)
      expected_a = 0.0_dp
      expected_c = 0.0_dp
      do j = 1, size(xi)
         payoff = max(100.0_dp - result%s(i) * (1.0_dp + xi(j)), 0.0_dp)
         expected_a = expected_a + xi(j) * payoff
         expected_c = expected_c + (1.0_dp - rho * xi(j)) * payoff / c2
      end do
      expected_a = expected_a / real(size(xi), dp) / m2
      expected_c = expected_c / real(size(xi), dp)
      call check_close(result%a(1, i), expected_a, 1.0e-12_dp, "a")
      call check_close(result%c(1, i), expected_c, 1.0e-12_dp, "c")
      call check_close(result%phi1(i), &
         (expected_a - expected_c * rho) / result%s(i), 1.0e-12_dp, "phi1")
   end do

   call check_close(result%option_value_at(1, 100.0_dp), result%c(1, 2), &
      0.0_dp, "grid option lookup")
   call check_close(result%initial_hedge_at(100.0_dp), result%phi1(2), &
      0.0_dp, "initial hedge lookup")

   print '(a)', 'test_hedging: PASS'

contains

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual
      real(dp), intent(in) :: expected
      real(dp), intent(in) :: tolerance
      character(len=*), intent(in) :: label

      if (abs(actual - expected) > tolerance) then
         print '(a)', 'failed: ' // label
         print '(a,es24.16)', 'actual:   ', actual
         print '(a,es24.16)', 'expected: ', expected
         error stop 1
      end if
   end subroutine check_close

end program test_hedging
