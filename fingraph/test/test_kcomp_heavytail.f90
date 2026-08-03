! SPDX-License-Identifier: GPL-3.0-only
program test_kcomp_heavytail
   use fingraph, only : dp, fingraph_result, L, block_diag, &
      learn_kcomp_heavytail_graph, relative_error, fg_ok
   use fingraph_linalg, only : symmetric_pseudoinverse, symmetric_eigen_jacobi
   use fingraph_rng, only : random_mvn, random_mvt
   implicit none
   real(dp) :: w(6), nu
   real(dp), allocatable :: block(:,:), truth(:,:), covariance(:,:), x(:,:)
   real(dp), allocatable :: values(:), vectors(:,:)
   type(fingraph_result) :: result
   integer :: status

   w = 1.0_dp/3.0_dp
   block = L(w)
   truth = block_diag(block,block)
   call symmetric_pseudoinverse(truth,covariance,status)
   call check(status == fg_ok,'pseudoinverse')

   call random_mvn(covariance,1000,x,seed=226)
   call learn_kcomp_heavytail_graph(x,result,k=2,rho=100.0_dp,reltol=1.0e-4_dp, &
      maxiter=4000)
   call check(result%status == fg_ok,'gaussian kcomp status')
   call check(result%convergence,'gaussian kcomp convergence')
   call check(relative_error(result%laplacian,truth) < 0.20_dp,'gaussian kcomp recovery')
   call validate_k_components(result%laplacian,2)

   nu = 4.0_dp
   call random_mvt(((nu-2.0_dp)/nu)*covariance,nu,1600,x,seed=227)
   call learn_kcomp_heavytail_graph(x,result,k=2,heavy_type='student',nu=nu, &
      rho=100.0_dp,reltol=1.0e-4_dp,maxiter=5000,record_objective=.true.)
   call check(result%status == fg_ok,'student kcomp status')
   call check(result%convergence,'student kcomp convergence')
   call check(relative_error(result%laplacian,truth) < 0.20_dp,'student kcomp recovery')
   call check(size(result%lagrangian) == result%iterations,'objective history')
   call check(size(result%beta_seq) == result%iterations,'beta history')
   call validate_k_components(result%laplacian,2)

   print '(a)', 'test_kcomp_heavytail: PASS'
contains
   subroutine validate_k_components(laplacian,k)
      real(dp), intent(in) :: laplacian(:,:)
      integer, intent(in) :: k
      real(dp) :: scale
      call symmetric_eigen_jacobi(laplacian,values,vectors,status)
      call check(status == fg_ok,'eigen status')
      scale = max(1.0_dp,maxval(abs(values)))
      call check(count(abs(values) < 2.0e-3_dp*scale) >= k,'component eigenvalues')
      call check(maxval(abs(sum(laplacian,dim=2))) < 1.0e-10_dp,'component row sums')
      call check(maxval(laplacian-diagonal_matrix(laplacian)) <= 1.0e-10_dp, &
         'nonpositive off diagonal')
   end subroutine validate_k_components

   pure function diagonal_matrix(a) result(d)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: d(size(a,1),size(a,2))
      integer :: i
      d = 0.0_dp
      do i = 1,min(size(a,1),size(a,2))
         d(i,i) = a(i,i)
      end do
   end function diagonal_matrix

   subroutine check(condition,message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message
      if (.not. condition) then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine check
end program test_kcomp_heavytail
