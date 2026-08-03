! SPDX-License-Identifier: GPL-3.0-only
program student_t_graph_example
   use fingraph, only : dp, fingraph_result, L, learn_regular_heavytail_graph
   use fingraph_linalg, only : symmetric_pseudoinverse
   use fingraph_rng, only : random_mvt
   implicit none
   real(dp) :: w(6), nu
   real(dp), allocatable :: laplacian(:,:), scatter(:,:), returns(:,:)
   type(fingraph_result) :: graph
   integer :: status

   w = 1.0_dp/3.0_dp
   nu = 4.0_dp
   laplacian = L(w)
   call symmetric_pseudoinverse(laplacian,scatter,status)
   scatter = ((nu-2.0_dp)/nu)*scatter
   call random_mvt(scatter,nu,1000,returns,seed=778)

   call learn_regular_heavytail_graph(returns,graph,heavy_type='student',nu=nu,rho=1.0_dp)
   write(*,'(a,l1,a,i0)') 'Converged: ',graph%convergence,', iterations: ',graph%iterations
   write(*,'(a,es12.4)') 'Final primal Laplacian residual: ', &
      graph%primal_lap_residual(graph%iterations)
end program student_t_graph_example
