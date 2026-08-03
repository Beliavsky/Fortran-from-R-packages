! SPDX-License-Identifier: GPL-3.0-only
program test_operators
   use spectral_graph_topology, only : dp, L, A, Linv, Ainv, Lstar, Astar, D, Dstar, &
      vecLmat, pairwise_matrix_rownorm2, upper_view_vec, metrics, graph_metrics, block_diag
   implicit none
   real(dp) :: w(6), expected_l(4,4), expected_a(4,4), y(3,3)
   real(dp), allocatable :: lw(:,:),aw(:,:),v(:,:),b(:,:)
   type(graph_metrics) :: m

   w=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   expected_l=reshape([6.0_dp,-1.0_dp,-2.0_dp,-3.0_dp, &
      -1.0_dp,10.0_dp,-4.0_dp,-5.0_dp, &
      -2.0_dp,-4.0_dp,12.0_dp,-6.0_dp, &
      -3.0_dp,-5.0_dp,-6.0_dp,14.0_dp],[4,4])
   expected_a=reshape([0.0_dp,1.0_dp,2.0_dp,3.0_dp, &
      1.0_dp,0.0_dp,4.0_dp,5.0_dp, &
      2.0_dp,4.0_dp,0.0_dp,6.0_dp, &
      3.0_dp,5.0_dp,6.0_dp,0.0_dp],[4,4])
   lw=L(w); aw=A(w)
   call check(maxval(abs(lw-expected_l))<1e-12_dp,'L operator')
   call check(maxval(abs(aw-expected_a))<1e-12_dp,'A operator')
   call check(maxval(abs(Linv(lw)-w))<1e-12_dp,'L inverse')
   call check(maxval(abs(Ainv(aw)-w))<1e-12_dp,'A inverse')
   call check(maxval(abs(sum(lw,dim=1)))<1e-12_dp,'L row sums')
   call check(maxval(abs(Lstar(lw)-matmul(transpose(vecLmat(4)),reshape(lw,[16]))))<1e-12_dp, &
      'Lstar consistency')
   call check(maxval(abs(Astar(aw)-2.0_dp*w))<1e-12_dp,'Astar consistency')
   call check(maxval(abs(Dstar(D(w))-Lstar(reshape([D(w),D(w),D(w),D(w)],[4,4]))))>=0.0_dp, &
      'degree operators execute')

   y=reshape([1.0_dp,3.0_dp,1.0_dp,2.0_dp,2.0_dp,1.0_dp,3.0_dp,0.0_dp,1.0_dp],[3,3])
   v=pairwise_matrix_rownorm2(y)
   call check(maxval(abs(v-reshape([0.0_dp,13.0_dp,5.0_dp,13.0_dp,0.0_dp,6.0_dp,5.0_dp,6.0_dp,0.0_dp],[3,3])))<1e-12_dp, &
      'pairwise row distances')
   call check(maxval(abs(upper_view_vec(lw)-[-1.0_dp,-2.0_dp,-3.0_dp,-4.0_dp,-5.0_dp,-6.0_dp]))<1e-12_dp, &
      'upper view')
   b=block_diag(L([1.0_dp]),L([2.0_dp]))
   call check(size(b,1)==4 .and. maxval(abs(sum(b,dim=1)))<1e-12_dp,'block diagonal')
   m=metrics(lw,lw,1e-8_dp)
   call check(abs(m%fscore-1.0_dp)<1e-12_dp .and. abs(m%accuracy-1.0_dp)<1e-12_dp,'metrics')
   print '(a)', 'test_operators: PASS'
contains
   subroutine check(condition,name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      if (.not.condition) then
         write(*,'(a)') 'FAIL: '//trim(name)
         error stop 1
      end if
   end subroutine check
end program test_operators
