! SPDX-License-Identifier: GPL-3.0-only
program demo_fingraph
   use fingraph, only : dp, fingraph_result, L, block_diag, &
      learn_connected_graph, learn_regular_heavytail_graph, &
      learn_kcomp_heavytail_graph, relative_error
   use fingraph_linalg, only : symmetric_pseudoinverse
   use fingraph_rng, only : random_mvt
   implicit none

   real(dp) :: w(6), nu
   real(dp), allocatable :: laplacian(:,:), covariance(:,:), x(:,:), two_blocks(:,:)
   type(fingraph_result) :: fit
   integer :: status

   w = 1.0_dp/3.0_dp
   laplacian = L(w)
   call symmetric_pseudoinverse(laplacian,covariance,status)

   call learn_connected_graph(covariance,fit,rho=100.0_dp)
   write(*,'(a,l1,a,i0,a,es11.3)') 'Connected Gaussian: converged=',fit%convergence, &
      ', iterations=',fit%iterations,', relative error=',relative_error(fit%laplacian,laplacian)

   nu = 4.0_dp
   call random_mvt(((nu-2.0_dp)/nu)*covariance,nu,1200,x,seed=2026)
   call learn_regular_heavytail_graph(x,fit,heavy_type='student',nu=nu,rho=1.0_dp)
   write(*,'(a,l1,a,i0,a,es11.3)') 'Connected Student-t: converged=',fit%convergence, &
      ', iterations=',fit%iterations,', relative error=',relative_error(fit%laplacian,laplacian)

   two_blocks = block_diag(laplacian,laplacian)
   call symmetric_pseudoinverse(two_blocks,covariance,status)
   call random_mvt(((nu-2.0_dp)/nu)*covariance,nu,1600,x,seed=2027)
   call learn_kcomp_heavytail_graph(x,fit,k=2,heavy_type='student',nu=nu, &
      rho=100.0_dp,reltol=1.0e-4_dp)
   write(*,'(a,l1,a,i0,a,es11.3)') 'Two-component Student-t: converged=',fit%convergence, &
      ', iterations=',fit%iterations,', relative error=',relative_error(fit%laplacian,two_blocks)
end program demo_fingraph
