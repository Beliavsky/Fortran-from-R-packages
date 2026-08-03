! SPDX-License-Identifier: Artistic-2.0
program test_linalg_structural
   use mts
   use test_support
   implicit none
   real(dp)::a(2,2),expected(2,2),x(12,2),yv(30),xv(30)
   real(dp),allocatable::ainv(:,:),root(:,:),v(:),restored(:,:),prod(:,:,:),pi(:,:,:),tab(:,:),ecc(:,:),pv(:,:),dy(:,:)
   real(dp)::pa(1,1,2),pb(1,1,2),phi(1,1,1),theta(1,1,1)
   integer::istat,i
   integer,allocatable::phi_id(:,:,:),theta_id(:,:,:)

   a=reshape([4.0_dp,2.0_dp,2.0_dp,3.0_dp],[2,2])
   call inverse_matrix(a,ainv,istat)
   expected=eye(2)
   call assert_true(istat==mts_success,'inverse status')
   call assert_matrix_close(matmul(a,ainv),expected,1.0e-11_dp,'matrix inverse')
   call matrix_sqrt_symmetric(a,root,istat)
   call assert_matrix_close(matmul(root,root),a,1.0e-9_dp,'symmetric square root')

   v=vech_matrix(a);call unvech_matrix(v,restored,istat)
   call assert_matrix_close(restored,a,1.0e-12_dp,'vech round trip')

   pa(1,1,:)=[1.0_dp,-0.5_dp];pb(1,1,:)=[1.0_dp,0.2_dp]
   call matrix_polynomial_product(pa,pb,prod)
   call assert_close(prod(1,1,1),1.0_dp,1.0e-12_dp,'polynomial c0')
   call assert_close(prod(1,1,2),-0.3_dp,1.0e-12_dp,'polynomial c1')
   call assert_close(prod(1,1,3),-0.1_dp,1.0e-12_dp,'polynomial c2')

   phi(1,1,1)=0.5_dp;theta(1,1,1)=0.2_dp
   call pi_weight_matrices(phi,theta,3,pi,istat)
   call assert_close(pi(1,1,0),1.0_dp,1.0e-12_dp,'pi zero')
   call assert_close(pi(1,1,1),-0.3_dp,1.0e-12_dp,'pi one')
   call assert_close(pi(1,1,2),-0.06_dp,1.0e-12_dp,'pi two')
   call kronecker_specification([2,1],phi_id,theta_id,istat)
   call assert_true(istat==mts_success.and.all(shape(phi_id)==[2,2,3]),'Kronecker specification')
   call assert_true(phi_id(1,2,0)==0.and.theta_id(1,1,0)==1,'Kronecker indicators')

   do i=1,12
      x(i,1)=real(i,dp);x(i,2)=2.0_dp*real(i,dp)
   end do
   dy=difference_matrix(x,d=1)
   call assert_true(all(shape(dy)==[11,2]),'difference shape')
   call assert_close(maxval(abs(dy(:,1)-1.0_dp)),0.0_dp,1.0e-12_dp,'first differences')

   xv(1)=sin(0.2_dp);yv(1)=0.8_dp*xv(1)
   do i=2,30
      xv(i)=sin(0.2_dp*real(i,dp));yv(i)=0.4_dp*yv(i-1)+0.8_dp*xv(i)
   end do
   call corner_table(yv,xv,3,3,tab,istat)
   call assert_true(istat==mts_success.and.all(shape(tab)==[3,3]),'corner table')
   x=0.0_dp
   do i=1,12
      x(i,1)=sin(real(i,dp));x(i,2)=cos(0.7_dp*real(i,dp))
   end do
   call extended_cross_correlation(x,1,2,ecc,pv,istat)
   call assert_true(istat==mts_success.and.all(shape(ecc)==[2,2]),'extended cross correlation')
   call assert_finite(reshape(ecc,[size(ecc)]),'finite extended cross correlations')

   print '(a)','test_linalg_structural: PASS'
end program test_linalg_structural
