! SPDX-License-Identifier: GPL-3.0-or-later
program semidefinite_program
   use cccp
   implicit none
   real(dp) :: q(2), flist(2,2,2), f0(2,2)
   type(cone_constraint) :: cone(1)
   type(cccp_solution) :: sol

   q = [-1.0_dp,-1.0_dp]
   flist(:,:,1) = reshape([1.0_dp,0.0_dp,0.0_dp,0.0_dp],[2,2])
   flist(:,:,2) = reshape([0.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
   f0 = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
   cone(1) = psdc(flist,f0)
   sol = cccp_solve(q,cones=cone)
   print '(a,a)','Status: ',trim(sol%status)
   print '(a,2f12.6)','x: ',sol%x
end program semidefinite_program
