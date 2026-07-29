! SPDX-License-Identifier: GPL-3.0-or-later
program second_order_cone
   use cccp
   implicit none
   real(dp) :: q(3), fmat(3,3), g(3), d(3)
   type(cone_constraint) :: cone(1)
   type(cccp_solution) :: sol

   q = [1.0_dp,0.0_dp,0.0_dp]
   fmat = 0.0_dp
   fmat(1,2) = 1.0_dp
   fmat(2,3) = 1.0_dp
   g = 0.0_dp
   d = [1.0_dp,0.0_dp,0.0_dp]
   cone(1) = socc(fmat,g,d,0.0_dp)
   sol = cccp_solve(q,cones=cone)
   print '(a,a)','Status: ',trim(sol%status)
   print '(a,3f12.6)','x: ',sol%x
end program second_order_cone
