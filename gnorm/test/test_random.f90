! SPDX-License-Identifier: GPL-2.0-or-later
program test_random
   use gnorm
   implicit none
   integer, parameter :: n = 150000
   real(dp), allocatable :: x(:), y(:)
   real(dp) :: sample_mean, sample_variance
   integer :: status

   x = rgnorm(n, mu=1.5_dp, alpha=sqrt(2.0_dp), beta=2.0_dp, &
      seed=123456_i8, status=status)
   if (status /= gnorm_success) error stop 'FAIL: random status'
   y = rgnorm(n, mu=1.5_dp, alpha=sqrt(2.0_dp), beta=2.0_dp, &
      seed=123456_i8)
   if (any(abs(x - y) > tiny(1.0_dp))) error stop 'FAIL: seeded reproducibility'

   sample_mean = sum(x) / real(n, dp)
   sample_variance = sum((x - sample_mean)**2) / real(n - 1, dp)
   if (abs(sample_mean - 1.5_dp) > 0.012_dp) error stop 'FAIL: random mean'
   if (abs(sample_variance - 1.0_dp) > 0.025_dp) error stop 'FAIL: random variance'

   call rgnorm_fill(y(1:1000), alpha=0.0_dp, status=status)
   if (status /= gnorm_invalid_argument) error stop 'FAIL: invalid RNG status'

   print '(a)', 'test_random: PASS'
end program test_random
