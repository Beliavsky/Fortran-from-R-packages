! SPDX-License-Identifier: GPL-3.0-only
program test_updates_and_data
   use spectral_graph_topology, only : dp, graph_result, L, block_diag, learn_k_component_graph, &
      learn_bipartite_k_component_graph, build_initial_graph
   use sgt_linalg, only : symmetric_pseudoinverse
   use sgt_updates, only : laplacian_u_update, bipartite_v_update, laplacian_lambda_update, &
      bipartite_psi_update
   implicit none
   real(dp) :: w1(3),w2(3),data(6,3)
   real(dp), allocatable :: l1(:,:),l2(:,:),lt(:,:),cov(:,:),u(:,:),v(:,:),lambda(:),psi(:),aff(:,:)
   type(graph_result) :: r1,r2
   integer :: status

   w1=[1.0_dp,0.5_dp,0.8_dp]; w2=[0.7_dp,1.1_dp,0.6_dp]
   l1=L(w1); l2=L(w2); lt=block_diag(l1,l2)
   call laplacian_u_update(lt,2,u,status)
   call check(size(u,1)==6 .and. size(u,2)==4,'U dimensions')
   call check(maxval(abs(matmul(transpose(u),u)-identity4()))<1e-8_dp,'U orthonormal')
   call laplacian_lambda_update(1e-6_dp,1e4_dp,10.0_dp,u,lt,2,lambda,status)
   call check(size(lambda)==4 .and. all(lambda>0.0_dp),'lambda update')

   call bipartite_v_update(L([1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp]),0,v,status)
   call check(size(v,2)==4,'V dimensions')
   call bipartite_psi_update(v,diag_to_adjacency(v),psi,status)
   call check(size(psi)==4 .and. abs(psi(1)+psi(4))<1e-8_dp,'psi symmetry')

   call symmetric_pseudoinverse(lt,cov,status)
   call learn_k_component_graph(cov,r1,k=2,beta=50.0_dp,maxiter=3000,reltol=1e-7_dp)
   call learn_bipartite_k_component_graph(cov,r2,k=2,nu=0.0_dp,beta=50.0_dp,maxiter=3000,reltol=1e-7_dp)
   call check(maxval(abs(r1%weights-r2%weights))<1e-7_dp,'nu=0 consistency')

   data=reshape([0.0_dp,0.1_dp,0.2_dp,3.0_dp,3.1_dp,3.2_dp, &
                 0.1_dp,0.0_dp,0.2_dp,3.2_dp,3.0_dp,3.1_dp, &
                 0.0_dp,0.2_dp,0.1_dp,3.1_dp,3.2_dp,3.0_dp],shape(data))
   call build_initial_graph(data,2,aff,status)
   call check(size(aff,1)==6 .and. minval(aff)>=0.0_dp,'initial affinity')
   call learn_k_component_graph(data,r1,is_data_matrix=.true.,m=2,k=2,beta=10.0_dp,maxiter=500)
   call check(size(r1%laplacian,1)==6,'data-matrix path')
   call check(maxval(abs(sum(r1%laplacian,dim=1)))<1e-8_dp,'data-matrix Laplacian')

   print '(a)', 'test_updates_and_data: PASS'
contains
   pure function identity4() result(a)
      real(dp) :: a(4,4)
      integer :: i
      a=0.0_dp
      do i=1,4; a(i,i)=1.0_dp; end do
   end function identity4

   function diag_to_adjacency(q) result(a)
      real(dp), intent(in) :: q(:,:)
      real(dp) :: a(size(q,1),size(q,1))
      real(dp) :: values(size(q,2))
      integer :: i
      values=[-2.0_dp,-1.0_dp,1.0_dp,2.0_dp]
      a=0.0_dp
      do i=1,size(values)
         a=a+values(i)*spread(q(:,i),2,size(q,1))*spread(q(:,i),1,size(q,1))
      end do
   end function diag_to_adjacency

   subroutine check(condition,name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not.condition) then
         write(*,'(a)') 'FAIL: '//trim(name)
         error stop 1
      end if
   end subroutine check
end program test_updates_and_data
