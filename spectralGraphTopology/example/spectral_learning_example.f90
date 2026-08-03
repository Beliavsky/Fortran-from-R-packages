! SPDX-License-Identifier: GPL-3.0-only
program spectral_learning_example
   use spectral_graph_topology, only : dp, graph_result, L, learn_k_component_graph, relative_error
   use sgt_linalg, only : symmetric_pseudoinverse
   implicit none
   real(dp) :: w(6)
   real(dp), allocatable :: truth(:,:),covariance(:,:)
   type(graph_result) :: estimate
   integer :: status

   w=[1.0_dp,0.8_dp,0.0_dp,0.7_dp,0.0_dp,0.9_dp]
   truth=L(w)
   call symmetric_pseudoinverse(truth,covariance,status)
   call learn_k_component_graph(covariance,estimate,beta=100.0_dp,maxiter=3000, &
      record_objective=.true.)

   print '(a,l1)', 'Converged: ',estimate%convergence
   print '(a,i0)', 'Iterations: ',estimate%iterations
   print '(a,es12.4)', 'Relative Laplacian error: ',relative_error(estimate%laplacian,truth)
   print '(a,*(f8.4,1x))', 'Estimated weights: ',estimate%weights
end program spectral_learning_example
