! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Risk 1.0 by Saralees Nadarajah and Stephen Chan.
! Copyright (c) 2017 Saralees Nadarajah and Stephen Chan.
program demo_risk
   use risk
   implicit none

   type(normal_distribution) :: standard_normal
   real(dp), parameter :: inf = huge(1.0_dp)
   real(dp), parameter :: alpha = 0.95_dp

   standard_normal = normal_distribution(mu=0.0_dp,sigma=1.0_dp)

   write(*,'(a,f12.6)') 'VaR:                    ',varg(standard_normal,alpha)
   write(*,'(a,f12.6)') 'Expected shortfall:     ',esg(standard_normal,alpha)
   write(*,'(a,f12.6)') 'Tail conditional median:',tcm(standard_normal,alpha)
   write(*,'(a,f12.6)') 'Expectation:            ',expect(standard_normal,-inf,inf)
   write(*,'(a,f12.6)') 'Omega at zero:          ',omegag(standard_normal,0.0_dp,-inf,inf)
end program demo_risk
