! SPDX-License-Identifier: GPL-3.0-only
program demo_spectral_graph_topology
   use spectral_graph_topology, only : dp, graph_result, L, learn_k_component_graph, &
      learn_smooth_approx_graph
   use sgt_linalg, only : symmetric_pseudoinverse
   implicit none
   real(dp) :: w(6), data(4,3)
   real(dp), allocatable :: laplacian(:,:)
   type(graph_result) :: graph, smooth
   real(dp), allocatable :: covariance(:,:)
   integer :: status

   w = [1.0_dp, 0.8_dp, 0.0_dp, 0.7_dp, 0.6_dp, 0.9_dp]
   laplacian = L(w)
   call symmetric_pseudoinverse(laplacian, covariance, status)
   call learn_k_component_graph(covariance, graph, beta=100.0_dp, maxiter=3000, reltol=1e-7_dp)
   print '(a,l1)', 'spectral learner converged: ', graph%convergence
   print '(a,i0)', 'iterations: ', graph%iterations
   print '(a,es12.4)', 'estimated Laplacian norm: ', sqrt(sum(graph%laplacian**2))

   data = reshape([0.0_dp,0.1_dp,2.0_dp,2.1_dp, &
                   0.2_dp,0.0_dp,1.9_dp,2.2_dp, &
                   0.1_dp,0.2_dp,2.1_dp,2.0_dp], shape(data))
   call learn_smooth_approx_graph(data, 1, smooth)
   print '(a,es12.4)', 'smooth graph row-sum error: ', maxval(abs(sum(smooth%laplacian,dim=2)))
end program demo_spectral_graph_topology
