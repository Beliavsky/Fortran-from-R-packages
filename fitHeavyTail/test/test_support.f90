! SPDX-License-Identifier: GPL-3.0-only
module test_support
   use fitheavytail, only: dp, random_mvt_identity
   implicit none
   private
   public :: check, make_data, max_abs, sample_cov_mle
contains
   subroutine check(condition,message)
      logical,intent(in)::condition
      character(len=*),intent(in)::message
      if(.not.condition) then
         write(*,'(a)') 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine check

   subroutine make_data(x)
      real(dp),intent(out)::x(:,:)
      real(dp),allocatable::z(:,:)
      integer::i,n,t
      t=size(x,1)
      n=size(x,2)
      allocate(z(t,n))
      call random_mvt_identity(t,n,7.0_dp,z,20260730)
      x=z
      if(n>=2) x(:,2)=0.35_dp*z(:,1)+0.9_dp*z(:,2)
      if(n>=3) x(:,3)=-0.2_dp*z(:,1)+0.25_dp*z(:,2)+0.8_dp*z(:,3)
      do i=1,t
         x(i,1)=x(i,1)+0.5_dp
         if(n>=2) x(i,2)=x(i,2)-0.3_dp
         if(n>=3) x(i,3)=x(i,3)+0.2_dp
      end do
   end subroutine make_data

   function max_abs(a) result(value)
      real(dp),intent(in)::a(:,:)
      real(dp)::value
      value=maxval(abs(a))
   end function max_abs

   function sample_cov_mle(x) result(cov)
      real(dp),intent(in)::x(:,:)
      real(dp)::cov(size(x,2),size(x,2)),mu(size(x,2)),row(size(x,2))
      integer::i
      mu=sum(x,dim=1)/real(size(x,1),dp)
      cov=0.0_dp
      do i=1,size(x,1)
         row=x(i,:)-mu
         cov=cov+spread(row,2,size(row))*spread(row,1,size(row))
      end do
      cov=cov/real(size(x,1),dp)
   end function sample_cov_mle
end module test_support
