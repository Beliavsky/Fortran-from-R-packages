! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
program basic_pricing
   use bcc1997
   implicit none

   type(bcc_result) :: result

   result = bcc(kappav=0.5_dp, kappar=0.0_dp, thetav=0.025_dp, &
      thetar=0.0_dp, sigmav=0.09_dp, sigmar=1.0e-7_dp, muj=0.0_dp, &
      sigmaj=1.0e-7_dp, rho=0.1_dp, lambda=0.0_dp, s0=100.0_dp, &
      k=100.0_dp, v0=0.04_dp, r0=0.01_dp, t=1.0_dp)

   print '(a,f12.6)', 'call: ', result%call
   print '(a,f12.6)', 'put:  ', result%put
end program basic_pricing
