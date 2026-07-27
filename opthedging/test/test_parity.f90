! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
program test_parity
   use opthedging, only : dp, hedging_iid, hedging_result
   implicit none

   integer :: k
   real(dp) :: factors(9)
   real(dp) :: log_returns(9)
   real(dp) :: parity_error
   real(dp) :: share_error
   type(hedging_result) :: call_result
   type(hedging_result) :: put_result
   type(hedging_result) :: singular_result

   factors = [0.78_dp, 0.86_dp, 0.94_dp, 0.99_dp, 1.01_dp, 1.05_dp, &
      1.11_dp, 1.19_dp, 1.28_dp]
   log_returns = log(factors)
   call_result = hedging_iid(log_returns, 0.75_dp, 105.0_dp, 0.04_dp, &
      .false., 5, 501, 30.0_dp, 200.0_dp)
   put_result = hedging_iid(log_returns, 0.75_dp, 105.0_dp, 0.04_dp, &
      .true., 5, 501, 30.0_dp, 200.0_dp)
   if (.not. call_result%ok .or. .not. put_result%ok) error stop 1

   do k = 1, size(call_result%c, 1)
      parity_error = maxval(abs(call_result%c(k, :) - put_result%c(k, :) - &
         (call_result%s - call_result%discounted_strike)))
      if (parity_error > 2.0e-11_dp) then
         print '(a,i0,a,es12.4)', 'put-call parity failed at period ', k, &
            ', error = ', parity_error
         error stop 1
      end if
   end do

   share_error = maxval(abs(call_result%phi1 - put_result%phi1 - 1.0_dp))
   if (share_error > 2.0e-11_dp) then
      print '(a,es12.4)', 'hedge parity error = ', share_error
      error stop 1
   end if

   singular_result = hedging_iid(log([1.10_dp, 1.10_dp, 1.10_dp]), &
      1.0_dp, 100.0_dp, 0.0_dp, .true., 2, 11, 50.0_dp, 150.0_dp)
   if (singular_result%ok) then
      print '(a)', 'singular return sample was not rejected'
      error stop 1
   end if

   print '(a)', 'test_parity: PASS'
end program test_parity
