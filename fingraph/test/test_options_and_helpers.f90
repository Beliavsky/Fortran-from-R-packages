! SPDX-License-Identifier: GPL-3.0-only
program test_options_and_helpers
   use fingraph, only : dp, fingraph_result, L, compute_student_weight, &
      learn_connected_graph, learn_regular_heavytail_graph, fg_ok, fg_invalid_input, &
      fg_no_convergence
   use fingraph_linalg, only : symmetric_pseudoinverse
   use fingraph_rng, only : random_mvn
   implicit none
   real(dp) :: w(6), degrees(4), expected
   real(dp), allocatable :: truth(:,:), covariance(:,:), x(:,:), z(:,:)
   type(fingraph_result) :: connected, regular, invalid
   integer :: status, j

   w = 1.0_dp/3.0_dp
   truth = L(w)
   call symmetric_pseudoinverse(truth,covariance,status)
   call check(status == fg_ok,'pseudoinverse')
   degrees = 1.0_dp

   call learn_connected_graph(covariance,connected,initialization='qp',degrees=degrees, &
      rho=100.0_dp,maxiter=3000)
   call check(connected%status == fg_ok,'qp initialization')

   call random_mvn(covariance,500,x,seed=331)
   allocate(z(size(x,1),size(x,2)))
   do j = 1,size(x,2)
      z(:,j) = (x(:,j)-sum(x(:,j))/real(size(x,1),dp)) &
         / sqrt(sum((x(:,j)-sum(x(:,j))/real(size(x,1),dp))**2) &
         / real(size(x,1)-1,dp))
   end do
   call learn_regular_heavytail_graph(z,regular,heavy_type='gaussian',rho=100.0_dp, &
      update_rho=.false.,maxiter=3000)
   call check(regular%status == fg_ok,'fixed rho')

   call learn_regular_heavytail_graph(x,invalid,heavy_type='student',nu=2.0_dp)
   call check(invalid%status == fg_invalid_input,'invalid nu')
   call learn_connected_graph(covariance,invalid,rho=100.0_dp,maxiter=1)
   call check(invalid%status == fg_no_convergence,'maxiter status')
   call check(invalid%iterations == 1,'maxiter count')
   expected = 8.0_dp/5.0_dp
   call check(abs(compute_student_weight([1.0_dp,1.0_dp],[0.5_dp,0.5_dp],4,4.0_dp) &
      - expected) < 1.0e-14_dp,'student weight')

   print '(a)', 'test_options_and_helpers: PASS'
contains
   subroutine check(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine check
end program test_options_and_helpers
