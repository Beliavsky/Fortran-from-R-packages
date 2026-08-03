! SPDX-License-Identifier: GPL-3.0-only
program operators_example
   use spectral_graph_topology, only : dp, L, A, Lstar, metrics, graph_metrics
   implicit none
   real(dp) :: w(6)
   real(dp), allocatable :: laplacian(:,:), adjacency(:,:), adjoint(:)
   type(graph_metrics) :: score

   w=[1.0_dp,0.0_dp,0.5_dp,0.8_dp,0.0_dp,1.2_dp]
   laplacian=L(w)
   adjacency=A(w)
   adjoint=Lstar(laplacian)
   score=metrics(adjacency,adjacency)

   print '(a)', 'Laplacian:'
   call print_matrix(laplacian)
   print '(a)', 'Adjacency:'
   call print_matrix(adjacency)
   print '(a,*(f8.3,1x))', 'Lstar(L(w)): ',adjoint
   print '(a,f6.3)', 'Self-comparison F-score: ',score%fscore
contains
   subroutine print_matrix(x)
      real(dp), intent(in) :: x(:,:)
      integer :: i
      do i=1,size(x,1)
         print '(*(f9.4,1x))',x(i,:)
      end do
   end subroutine print_matrix
end program operators_example
