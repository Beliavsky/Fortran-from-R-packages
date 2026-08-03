! SPDX-License-Identifier: GPL-3.0-only
program test_regular_heavytail
   use fingraph, only : dp, fingraph_result, L, learn_regular_heavytail_graph, &
      relative_error, fscore, fg_ok
   use fingraph_linalg, only : symmetric_pseudoinverse
   use fingraph_rng, only : random_mvn, random_mvt
   implicit none
   real(dp) :: w(6), nu
   real(dp), allocatable :: truth(:,:), covariance(:,:), x(:,:)
   type(fingraph_result) :: result
   integer :: status

   w = 1.0_dp/3.0_dp
   truth = L(w)
   call symmetric_pseudoinverse(truth,covariance,status)
   call check(status == fg_ok,'pseudoinverse')

   call random_mvn(covariance,600,x,seed=124)
   call learn_regular_heavytail_graph(x,result,heavy_type='gaussian',rho=100.0_dp, &
      maxiter=3000)
   call check(result%status == fg_ok,'gaussian status')
   call check(result%convergence,'gaussian convergence')
   call check(relative_error(result%laplacian,truth) < 0.18_dp,'gaussian recovery')
   call check(fscore(truth,result%laplacian,0.1_dp) > 0.99_dp,'gaussian support')
   call check(size(result%lagrangian) == result%iterations,'gaussian history')

   nu = 4.0_dp
   call random_mvt(((nu-2.0_dp)/nu)*covariance,nu,1200,x,seed=125)
   call learn_regular_heavytail_graph(x,result,heavy_type='student',nu=nu,rho=1.0_dp, &
      maxiter=4000)
   call check(result%status == fg_ok,'student status')
   call check(result%convergence,'student convergence')
   call check(relative_error(result%laplacian,truth) < 0.18_dp,'student recovery')
   call check(fscore(truth,result%laplacian,0.02_dp) > 0.99_dp,'student support')
   call check(all(result%primal_lap_residual >= 0.0_dp),'primal history')
   call check(all(result%elapsed_time >= 0.0_dp),'elapsed history')

   print '(a)', 'test_regular_heavytail: PASS'
contains
   subroutine check(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine check
end program test_regular_heavytail
