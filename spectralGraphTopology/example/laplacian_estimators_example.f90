! SPDX-License-Identifier: GPL-3.0-only
program laplacian_estimators_example
   use spectral_graph_topology, only : dp, graph_result, L, learn_laplacian_gle_mm, &
      learn_laplacian_gle_admm, learn_combinatorial_graph_laplacian, relative_error
   use sgt_linalg, only : symmetric_pseudoinverse
   implicit none
   real(dp) :: w(6),mask(4,4)
   real(dp), allocatable :: truth(:,:),covariance(:,:)
   type(graph_result) :: mm,admm,cgl
   integer :: status

   w=[1.0_dp,0.8_dp,0.0_dp,0.7_dp,0.0_dp,0.9_dp]
   truth=L(w)
   call symmetric_pseudoinverse(truth,covariance,status)
   mask=merge(1.0_dp,0.0_dp,truth<0.0_dp)

   call learn_laplacian_gle_mm(covariance,mm,a_mask=mask,maxiter=2000)
   call learn_laplacian_gle_admm(covariance,admm,a_mask=mask,maxiter=3000)
   call learn_combinatorial_graph_laplacian(covariance,cgl,a_mask=mask,max_cycle=500)

   print '(a,es12.4)', 'GLE-MM relative error: ',relative_error(mm%laplacian,truth)
   print '(a,es12.4)', 'GLE-ADMM relative error: ',relative_error(admm%laplacian,truth)
   print '(a,es12.4)', 'CGL relative error: ',relative_error(cgl%laplacian,truth)
end program laplacian_estimators_example
