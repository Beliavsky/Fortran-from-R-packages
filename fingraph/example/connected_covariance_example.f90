! SPDX-License-Identifier: GPL-3.0-only
program connected_covariance_example
   use fingraph, only : dp, fingraph_result, learn_connected_graph
   implicit none
   real(dp) :: covariance(4,4)
   type(fingraph_result) :: graph
   integer :: i

   covariance = reshape([ &
      0.75_dp,-0.25_dp,-0.25_dp,-0.25_dp, &
     -0.25_dp, 0.75_dp,-0.25_dp,-0.25_dp, &
     -0.25_dp,-0.25_dp, 0.75_dp,-0.25_dp, &
     -0.25_dp,-0.25_dp,-0.25_dp, 0.75_dp], [4,4])

   call learn_connected_graph(covariance,graph,rho=100.0_dp)
   write(*,'(a,l1,a,i0)') 'Converged: ',graph%convergence,', iterations: ',graph%iterations
   write(*,'(a)') 'Estimated Laplacian:'
   do i = 1,size(graph%laplacian,1)
      write(*,'(*(f10.5,1x))') graph%laplacian(i,:)
   end do
end program connected_covariance_example
