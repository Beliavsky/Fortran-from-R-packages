! SPDX-License-Identifier: GPL-3.0-only
program test_connected_graph
   use fingraph, only : dp, fingraph_result, L, learn_connected_graph, &
      relative_error, fscore, fg_ok
   use fingraph_linalg, only : symmetric_pseudoinverse, frobenius_norm
   implicit none
   real(dp) :: w(6)
   real(dp), allocatable :: truth(:,:), covariance(:,:)
   type(fingraph_result) :: result
   integer :: status

   w = 1.0_dp/3.0_dp
   truth = L(w)
   call symmetric_pseudoinverse(truth,covariance,status)
   call check(status == fg_ok,'pseudoinverse')

   call learn_connected_graph(covariance,result,rho=100.0_dp,maxiter=2000)
   call check(result%status == fg_ok,'connected status')
   call check(result%convergence,'connected convergence')
   call check(result%iterations > 5,'connected iterations')
   call check(relative_error(result%laplacian,truth) < 1.0e-2_dp,'connected recovery')
   call check(fscore(truth,result%laplacian,0.1_dp) > 0.99_dp,'connected support')
   call check(maxval(abs(sum(result%laplacian,dim=2))) < 1.0e-10_dp,'zero row sums')
   call check(frobenius_norm(result%laplacian-transpose(result%laplacian)) < 1.0e-12_dp, &
      'symmetry')

   print '(a)', 'test_connected_graph: PASS'
contains
   subroutine check(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine check
end program test_connected_graph
