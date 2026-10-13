! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program inequality_demo
   use ineq_api, only : dp, gini, atkinson, theil, foster, sst, lc, lorenz_curve, ineq_ok, lc_lognorm
   implicit none
   real(dp), parameter :: incomes(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 10.0_dp]
   type(lorenz_curve) :: curve
   integer :: i
   print '(a,f10.6)', 'Gini:               ', gini(incomes)
   print '(a,f10.6)', 'Atkinson (0.5):     ', atkinson(incomes)
   print '(a,f10.6)', 'Theil T:            ', theil(incomes)
   print '(a,f10.6)', 'Headcount below 5:  ', foster(incomes, 5.0_dp)
   print '(a,f10.6)', 'SST poverty below 5:', sst(incomes, 5.0_dp)
   curve = lc(incomes)
   if (curve%status /= ineq_ok) error stop 'Lorenz calculation failed'
   print '(a)', 'Population share, empirical Lorenz, lognormal Lorenz'
   do i = 1, size(curve%p)
      print '(3f12.6)', curve%p(i), curve%l(i), lc_lognorm(curve%p(i))
   end do
end program inequality_demo
