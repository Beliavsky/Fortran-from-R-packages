! SPDX-License-Identifier: GPL-3.0-only
program test_moments
   use highorderportfolios
   implicit none
   real(dp) :: x(8,3),w(3),m(4),y(8),meanx(3),expected(4),tensor_m3,tensor_m4
   type(sample_moments) :: s,sa
   integer :: i,j,k,l

   x(:,1)=[-0.08_dp,-0.03_dp,0.01_dp,0.02_dp,0.04_dp,0.05_dp,0.09_dp,0.14_dp]
   x(:,2)=[ 0.03_dp,-0.01_dp,0.02_dp,0.01_dp,0.00_dp,0.04_dp,0.02_dp,0.06_dp]
   x(:,3)=[-0.02_dp,0.01_dp,-0.01_dp,0.03_dp,0.02_dp,-0.02_dp,0.05_dp,0.04_dp]
   w=[0.2_dp,0.3_dp,0.5_dp]
   call estimate_sample_moments(x,s,store_tensors=.true.)
   call check(s%status==hop_success,'sample status')
   m=eval_portfolio_moments(w,s)
   meanx=sum(x,dim=1)/8.0_dp
   y=matmul(x-spread(meanx,1,8),w)
   expected(1)=dot_product(w,meanx)
   expected(2)=sum(y*y)/7.0_dp
   expected(3)=sum(y**3)/8.0_dp
   expected(4)=sum(y**4)/8.0_dp
   call check(maxval(abs(m-expected))<1.0e-12_dp,'direct moments')

   tensor_m3=0.0_dp
   tensor_m4=0.0_dp
   do i=1,3
      do j=1,3
         do k=1,3
            tensor_m3=tensor_m3+w(i)*w(j)*w(k)*s%coskewness(i,j,k)
            do l=1,3
               tensor_m4=tensor_m4+w(i)*w(j)*w(k)*w(l)*s%cokurtosis(i,j,k,l)
            end do
         end do
      end do
   end do
   call check(abs(tensor_m3-m(3))<1.0e-12_dp,'coskewness tensor')
   call check(abs(tensor_m4-m(4))<1.0e-12_dp,'cokurtosis tensor')

   call estimate_sample_moments(x,sa,adjust_magnitude=.true.)
   w=1.0_dp/3.0_dp
   m=eval_portfolio_moments(w,sa)
   call check(maxval(abs(abs(m)-1.0_dp))<1.0e-12_dp,'magnitude adjustment')
   print '(a)', 'test_moments: PASS'
contains
   subroutine check(ok,label)
      logical,intent(in)::ok
      character(len=*),intent(in)::label
      if(.not.ok) then
         print '(a)', 'FAIL: '//trim(label)
         error stop 1
      end if
   end subroutine check
end program test_moments
