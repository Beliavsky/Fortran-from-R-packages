! SPDX-License-Identifier: GPL-3.0-only
program test_other_learning
   use spectral_graph_topology, only : dp, graph_result, L, relative_error, &
      learn_smooth_approx_graph, cluster_k_component_graph, learn_smooth_graph, &
      learn_graph_sigrep, learn_laplacian_gle_mm, learn_laplacian_gle_admm, &
      learn_combinatorial_graph_laplacian
   use sgt_linalg, only : symmetric_pseudoinverse
   implicit none
   real(dp) :: data(6,3),signals(4,12),w(6),mask(4,4)
   real(dp), allocatable :: ltrue(:,:),covariance(:,:)
   type(graph_result) :: res
   integer :: i,status

   data=reshape([0.0_dp,0.1_dp,0.2_dp,3.0_dp,3.1_dp,3.2_dp, &
                 0.1_dp,0.0_dp,0.2_dp,3.2_dp,3.0_dp,3.1_dp, &
                 0.0_dp,0.2_dp,0.1_dp,3.1_dp,3.2_dp,3.0_dp],shape(data))
   call learn_smooth_approx_graph(data,2,res)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-12_dp,'smooth approximate row sums')
   call check(maxval(abs(res%laplacian-transpose(res%laplacian)))<1e-12_dp,'smooth approximate symmetry')

   call cluster_k_component_graph(data,res,k=2,m=2,maxiter=30)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-9_dp,'CLR row sums')
   call check(minval(res%adjacency)>=-1e-12_dp,'CLR nonnegative adjacency')

   do i=1,12
      signals(:,i)=[sin(0.2_dp*real(i,dp)),sin(0.2_dp*real(i,dp))+0.02_dp, &
                    cos(0.17_dp*real(i,dp)),cos(0.17_dp*real(i,dp))+0.02_dp]
   end do
   call learn_smooth_graph(signals,res,maxiter=2000,tolerance=1e-6_dp)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-8_dp,'Kalofolias row sums')
   call check(minval(res%weights)>=-1e-10_dp,'Kalofolias weights')

   call learn_graph_sigrep(signals,res,maxiter=50,ftol=1e-5_dp)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-8_dp,'signal representation row sums')
   call check(abs(sum([(res%laplacian(i,i),i=1,4)])-4.0_dp)<1e-5_dp,'signal representation trace')
   call check(allocated(res%smoothed_data),'signal representation data')

   w=[1.0_dp,0.8_dp,0.0_dp,0.7_dp,0.0_dp,0.9_dp]
   ltrue=L(w); call symmetric_pseudoinverse(ltrue,covariance,status)
   mask=merge(1.0_dp,0.0_dp,ltrue<0.0_dp)

   call learn_laplacian_gle_mm(covariance,res,a_mask=mask,maxiter=2000,reltol=1e-7_dp,record_objective=.true.)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-7_dp,'GLE-MM row sums')
   call check(relative_error(res%laplacian,ltrue)<0.35_dp,'GLE-MM recovery')

   call learn_laplacian_gle_admm(covariance,res,a_mask=mask,maxiter=3000,reltol=1e-7_dp,record_objective=.true.)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-6_dp,'GLE-ADMM row sums')
   call check(relative_error(res%laplacian,ltrue)<0.4_dp,'GLE-ADMM recovery')

   call learn_combinatorial_graph_laplacian(covariance,res,a_mask=mask,max_cycle=500,reltol=1e-7_dp, &
      record_objective=.true.)
   call check(maxval(abs(sum(res%laplacian,dim=1)))<1e-5_dp,'CGL row sums')
   call check(relative_error(res%laplacian,ltrue)<0.5_dp,'CGL recovery')

   print '(a)', 'test_other_learning: PASS'
contains
   subroutine check(condition,name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not.condition) then
         write(*,'(a)') 'FAIL: '//trim(name)
         error stop 1
      end if
   end subroutine check
end program test_other_learning
