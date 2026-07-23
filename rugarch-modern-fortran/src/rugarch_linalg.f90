! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_linalg
   use rugarch_kinds, only : dp
   implicit none
   private

   public :: invert_matrix, solve_linear_system, covariance_matrix
   public :: symmetric_trace_product

contains

   subroutine invert_matrix(a, ainv, info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(size(a,1),size(a,2))
      integer, intent(out) :: info
      real(dp), allocatable :: aug(:,:), rowtmp(:)
      real(dp) :: pivot
      integer :: n, i, j, k, ipiv

      n=size(a,1)
      info=0
      ainv=0.0_dp
      if(size(a,2)/=n)then
         info=1
         return
      end if
      allocate(aug(n,2*n),rowtmp(2*n))
      aug(:,1:n)=a
      aug(:,n+1:2*n)=0.0_dp
      do i=1,n
         aug(i,n+i)=1.0_dp
      end do

      do i=1,n
         ipiv=i
         do k=i+1,n
            if(abs(aug(k,i))>abs(aug(ipiv,i)))ipiv=k
         end do
         if(abs(aug(ipiv,i))<=100.0_dp*tiny(1.0_dp))then
            info=2
            return
         end if
         if(ipiv/=i)then
            rowtmp=aug(i,:)
            aug(i,:)=aug(ipiv,:)
            aug(ipiv,:)=rowtmp
         end if
         pivot=aug(i,i)
         aug(i,:)=aug(i,:)/pivot
         do j=1,n
            if(j==i)cycle
            pivot=aug(j,i)
            if(abs(pivot)>tiny(1.0_dp))aug(j,:)=aug(j,:)-pivot*aug(i,:)
         end do
      end do
      ainv=aug(:,n+1:2*n)
   end subroutine invert_matrix

   subroutine solve_linear_system(a,b,x,info)
      real(dp),intent(in)::a(:,:),b(:)
      real(dp),intent(out)::x(size(b))
      integer,intent(out)::info
      real(dp),allocatable::ainv(:,:)
      integer::n
      n=size(b)
      if(size(a,1)/=n .or. size(a,2)/=n)then
         info=1;x=0.0_dp;return
      end if
      allocate(ainv(n,n))
      call invert_matrix(a,ainv,info)
      if(info==0)then
         x=matmul(ainv,b)
      else
         x=0.0_dp
      end if
   end subroutine solve_linear_system

   function covariance_matrix(x,center) result(cov)
      real(dp),intent(in)::x(:,:)
      logical,intent(in),optional::center
      real(dp),allocatable::cov(:,:)
      real(dp),allocatable::xc(:,:)
      logical::do_center
      integer::n,j

      n=size(x,1)
      allocate(cov(size(x,2),size(x,2)),xc(size(x,1),size(x,2)))
      xc=x
      do_center=.true.;if(present(center))do_center=center
      if(do_center .and. n>0)then
         do j=1,size(x,2)
            xc(:,j)=xc(:,j)-sum(xc(:,j))/real(n,dp)
         end do
      end if
      if(n>1)then
         cov=matmul(transpose(xc),xc)/real(n-1,dp)
      else
         cov=0.0_dp
      end if
   end function covariance_matrix

   pure function symmetric_trace_product(a,b) result(value)
      real(dp),intent(in)::a(:,:),b(:,:)
      real(dp)::value
      integer::i,j,n
      n=min(size(a,1),min(size(a,2),min(size(b,1),size(b,2))))
      value=0.0_dp
      do i=1,n
         do j=1,n
            value=value+a(i,j)*b(j,i)
         end do
      end do
   end function symmetric_trace_product

end module rugarch_linalg
