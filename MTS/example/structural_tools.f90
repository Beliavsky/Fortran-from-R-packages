! SPDX-License-Identifier: Artistic-2.0
program structural_tools
   use mts
   implicit none
   real(dp)::phi(2,2,1),theta(2,2,1),a(2,2,2),b(2,2,2)
   real(dp),allocatable::pi(:,:,:),product(:,:,:)
   integer::status

   phi(:,:,1)=reshape([0.4_dp,0.0_dp,-0.1_dp,0.3_dp],[2,2])
   theta(:,:,1)=reshape([0.2_dp,0.0_dp,0.0_dp,0.1_dp],[2,2])
   call pi_weight_matrices(phi,theta,5,pi,status)
   a(:,:,1)=eye(2);a(:,:,2)=-phi(:,:,1)
   b(:,:,1)=eye(2);b(:,:,2)=-theta(:,:,1)
   call matrix_polynomial_product(a,b,product)

   print '(a)','# Pi-weight matrix at lag one'
   print '(2(1x,f9.5))',pi(:,:,1)
   print '(a)','# polynomial product coefficient at lag two'
   print '(2(1x,f9.5))',product(:,:,3)
end program structural_tools
