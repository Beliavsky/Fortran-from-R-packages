! SPDX-License-Identifier: GPL-3.0-only
program k_component_graph_example
   use fingraph, only : dp, fingraph_result, L, block_diag, &
      learn_kcomp_heavytail_graph
   use fingraph_linalg, only : symmetric_pseudoinverse, symmetric_eigen_jacobi
   use fingraph_rng, only : random_mvn
   implicit none
   real(dp) :: w(6)
   real(dp), allocatable :: block(:,:), laplacian(:,:), covariance(:,:), returns(:,:)
   real(dp), allocatable :: eigenvalues(:), eigenvectors(:,:)
   type(fingraph_result) :: graph
   integer :: status

   w = 1.0_dp/3.0_dp
   block = L(w)
   laplacian = block_diag(block,block)
   call symmetric_pseudoinverse(laplacian,covariance,status)
   call random_mvn(covariance,1000,returns,seed=779)

   call learn_kcomp_heavytail_graph(returns,graph,k=2,rho=100.0_dp, &
      reltol=1.0e-4_dp,record_objective=.true.)
   call symmetric_eigen_jacobi(graph%laplacian,eigenvalues,eigenvectors,status)
   write(*,'(a,l1,a,i0)') 'Converged: ',graph%convergence,', iterations: ',graph%iterations
   write(*,'(a,*(es11.3,1x))') 'Smallest eigenvalues: ',eigenvalues(1:3)
end program k_component_graph_example
