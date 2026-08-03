! SPDX-License-Identifier: MIT
program test_special
   use invgamstochvol
   implicit none

   real(dp) :: value, expected
   integer :: status

   value = ourgeo(1.5_dp, 1.2_dp, 1.2_dp, 0.3_dp, 200, status)
   expected = 0.7_dp**(-1.5_dp)
   call check(status == invgam_success, 'hypergeometric status')
   call check(abs(value - expected) < 2.0e-13_dp, 'hypergeometric identity')

   value = ourgeo(-2.0_dp, 3.0_dp, 4.0_dp, 0.5_dp, 20, status)
   call check(abs(value - 0.4_dp) < 2.0e-14_dp, 'terminating hypergeometric series')

   value = log_rising_factorial(2.5_dp, 3, status)
   expected = log(2.5_dp * 3.5_dp * 4.5_dp)
   call check(status == invgam_success, 'rising factorial status')
   call check(abs(value - expected) < 1.0e-14_dp, 'log rising factorial')

   print '(a)', 'test_special: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_special
