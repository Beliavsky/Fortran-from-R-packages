! SPDX-License-Identifier: GPL-3.0-or-later
program demo_cccp
   use cccp_api
   implicit none
   real(dp) :: p(2,2), q(2), a(1,2), b(1), cov(3,3), x0(3), budgets(3)
   type(cccp_solution) :: sol

   p = 2.0_dp * reshape([2.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2])
   q = [1.0_dp,1.0_dp]
   a = reshape([1.0_dp,1.0_dp],[1,2])
   b = 1.0_dp
   sol = cccp(p,q,a,b)
   print '(a,2f12.6)','Equality-constrained QP solution: ',sol%x
   print '(a,f12.6)','Objective: ',sol%state%pobj

   cov = reshape([0.04_dp,0.006_dp,0.004_dp,0.006_dp,0.09_dp,0.01_dp, &
      0.004_dp,0.01_dp,0.16_dp],[3,3])
   x0 = [0.4_dp,0.35_dp,0.25_dp]
   budgets = 1.0_dp/3.0_dp
   sol = rp(x0,cov,budgets)
   print '(a,3f12.6)','Risk-parity weights: ',sol%x
end program demo_cccp
