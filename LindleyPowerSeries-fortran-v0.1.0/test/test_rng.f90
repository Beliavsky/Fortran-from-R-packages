! SPDX-License-Identifier: GPL-2.0-or-later
program test_rng
   use lindley_power_series
   implicit none

   real(dp), allocatable :: x(:)
   integer :: failures

   failures = 0

   x = rlindleygeometric(200, 1.0_dp, 0.5_dp)
   call check_sample(x, 'geometric')
   x = rlindleylogarithmic(200, 1.0_dp, 0.5_dp)
   call check_sample(x, 'logarithmic')
   x = rlindleynb(200, 1.0_dp, 0.5_dp, 3)
   call check_sample(x, 'negative-binomial')
   x = rlindleybinomial(200, 1.0_dp, 1.5_dp, 3)
   call check_sample(x, 'binomial')
   x = rlindleypoisson(200, 1.0_dp, 2.0_dp)
   call check_sample(x, 'poisson')

   if (failures /= 0) error stop 1
   print '(a)', 'test_rng: PASS'

contains

   subroutine check_sample(y, name)
      real(dp), intent(in) :: y(:)
      character(len=*), intent(in) :: name

      if (size(y) /= 200 .or. any(y < 0.0_dp) .or. any(y >= huge(1.0_dp))) then
         print '(a)', trim(name)//': FAIL'
         failures = failures + 1
      end if
   end subroutine check_sample

end program test_rng
