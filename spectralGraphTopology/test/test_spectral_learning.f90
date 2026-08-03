! SPDX-License-Identifier: GPL-3.0-only
program test_spectral_learning
   use spectral_graph_topology, only : dp, graph_result, L, relative_error, &
      learn_k_component_graph, learn_cospectral_graph, learn_bipartite_graph, &
      learn_bipartite_k_component_graph, sgt_ok, sgt_no_convergence
   use sgt_linalg, only : symmetric_pseudoinverse, symmetric_eigen_jacobi
   implicit none
   real(dp) :: w(6)
   real(dp), allocatable :: ltrue(:,:),covariance(:,:),eig(:),vec(:,:),fixed(:)
   type(graph_result) :: res
   integer :: status

   w=[1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp]
   ltrue=L(w)
   call symmetric_pseudoinverse(ltrue,covariance,status)

   call learn_k_component_graph(covariance,res,beta=100.0_dp,maxiter=3000,reltol=1e-7_dp,abstol=1e-9_dp)
   call check(res%status==sgt_ok .or. res%status==sgt_no_convergence,'k-component status')
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-9_dp,'k-component Laplacian')
   call check(minval(res%weights)>=-1e-12_dp,'k-component nonnegative weights')
   call check(relative_error(res%laplacian,ltrue)<0.2_dp,'k-component recovery')

   call symmetric_eigen_jacobi(ltrue,eig,vec,status)
   fixed=eig(2:size(eig))
   call learn_cospectral_graph(covariance,fixed,res,beta=100.0_dp,maxiter=3000,reltol=1e-7_dp)
   call check(relative_error(res%laplacian,ltrue)<0.2_dp,'cospectral recovery')

   call learn_bipartite_graph(covariance,res,use_qp=.true.,nu=1000.0_dp,maxiter=3000,reltol=1e-6_dp)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-8_dp,'bipartite Laplacian')
   call check_bipartite_spectrum(res%adjacency,'bipartite spectrum')

   call learn_bipartite_k_component_graph(covariance,res,use_qp=.true.,beta=100.0_dp,nu=100.0_dp, &
      maxiter=3000,reltol=1e-6_dp,record_objective=.true.)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-8_dp,'joint Laplacian')
   call check_bipartite_spectrum(res%adjacency,'joint bipartite spectrum')
   call check(allocated(res%objective),'joint objective')

   print '(a)', 'test_spectral_learning: PASS'
contains

   subroutine check_bipartite_spectrum(adjacency,name)
      real(dp), intent(in) :: adjacency(:,:)
      character(len=*), intent(in) :: name
      real(dp), allocatable :: values(:),vectors(:,:)
      real(dp) :: scale,symmetry_error
      integer :: eigen_status,n

      call symmetric_eigen_jacobi(adjacency,values,vectors,eigen_status)
      call check(eigen_status==sgt_ok .or. eigen_status==sgt_no_convergence, &
         trim(name)//' eigensolver status')
      n=size(values)
      call check(n>0,trim(name)//' nonempty')
      scale=max(1.0_dp,maxval(abs(values)))
      symmetry_error=maxval(abs(values+values(n:1:-1)))/scale
      call check(symmetry_error<1e-7_dp,name)
   end subroutine check_bipartite_spectrum

   subroutine check(condition,name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not.condition) then
         write(*,'(a)') 'FAIL: '//trim(name)
         error stop 1
      end if
   end subroutine check
end program test_spectral_learning
