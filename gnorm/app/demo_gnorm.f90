! SPDX-License-Identifier: GPL-2.0-or-later
program demo_gnorm
   use gnorm
   implicit none
   real(dp), allocatable :: draws(:)
   real(dp) :: p, q

   p = 0.975_dp
   q = qgnorm(p, alpha=sqrt(2.0_dp), beta=2.0_dp)
   print '(a,f10.6)', 'standard-normal 97.5% quantile: ', q
   print '(a,f10.6)', 'CDF at that quantile:           ', &
      pgnorm(q, alpha=sqrt(2.0_dp), beta=2.0_dp)

   draws = rgnorm(8, mu=1.0_dp, alpha=0.8_dp, beta=1.3_dp, seed=2026_i8)
   print '(a,8f10.5)', 'seeded draws:', draws
end program demo_gnorm
