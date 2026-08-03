! SPDX-License-Identifier: GPL-2.0-or-later
program quantile_example
   use gnorm
   implicit none
   real(dp), parameter :: p(5) = [0.01_dp, 0.10_dp, 0.50_dp, 0.90_dp, 0.99_dp]
   real(dp) :: q(size(p))
   integer :: i

   q = qgnorm(p, mu=2.0_dp, alpha=1.5_dp, beta=1.2_dp)
   do i = 1, size(p)
      print '(f7.3,2es20.10)', p(i), q(i), &
         pgnorm(q(i), mu=2.0_dp, alpha=1.5_dp, beta=1.2_dp)
   end do
end program quantile_example
