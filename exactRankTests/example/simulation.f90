! SPDX-License-Identifier: GPL-2.0-or-later
program simulation_example
   use exact_rank_tests
   implicit none
   type(rank_distribution) :: law
   real(dp) :: uniforms(3, 8)
   integer :: i, j, status
   ! Enumerating these uniforms realizes each sign configuration exactly once.
   do j = 1, 8
      do i = 1, 3
         uniforms(i, j) = merge(0.25_dp, 0.75_dp, btest(j - 1, i - 1))
      end do
   end do
   call simulate_distribution([1.0_dp, 2.0_dp, 3.0_dp], 3, uniforms, law, status)
   if (status /= rank_ok) error stop 'Simulation failed'
   print '(a,f8.5)', 'P(T=3) = ', dperm(3.0_dp, law)
   print '(a,f8.5)', 'P(T>=5) = ', pperm(5.0_dp, law, 'greater')
   print '(a,f8.5)', 'Median = ', qperm(0.5_dp, law)
end program simulation_example
