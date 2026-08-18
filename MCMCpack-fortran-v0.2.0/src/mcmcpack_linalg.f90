! SPDX-License-Identifier: GPL-3.0-only
module mcmcpack_linalg
   use mcmcpack_kinds, only : dp
   implicit none
   private
   public :: chol_lower, solve_lower, solve_upper, solve_spd, inv_spd
   public :: outer_product, trace_mat, logdet_spd, symmetrize
   public :: least_squares_normal, covariance_matrix
contains
   subroutine chol_lower(a,l,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: l(size(a,1),size(a,2))
      integer, intent(out) :: info
      integer :: i,j,k,n
      real(dp) :: s
      n = size(a,1); l = 0.0_dp; info = 0
      if (size(a,2) /= n) then; info = -1; return; end if
      do i=1,n
         do j=1,i
            s = a(i,j)
            do k=1,j-1
               s = s-l(i,k)*l(j,k)
            end do
            if (i == j) then
               if (s <= 0.0_dp) then; info = i; return; end if
               l(i,j) = sqrt(s)
            else
               l(i,j) = s/l(j,j)
            end if
         end do
      end do
   end subroutine chol_lower

   subroutine solve_lower(l,b,x)
      real(dp), intent(in) :: l(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      integer :: i,j,n
      real(dp) :: s
      n=size(b)
      do i=1,n
         s=b(i)
         do j=1,i-1; s=s-l(i,j)*x(j); end do
         x(i)=s/l(i,i)
      end do
   end subroutine solve_lower

   subroutine solve_upper(u,b,x)
      real(dp), intent(in) :: u(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      integer :: i,j,n
      real(dp) :: s
      n=size(b)
      do i=n,1,-1
         s=b(i)
         do j=i+1,n; s=s-u(i,j)*x(j); end do
         x(i)=s/u(i,i)
      end do
   end subroutine solve_upper

   subroutine solve_spd(a,b,x,info)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      integer, intent(out) :: info
      real(dp), allocatable :: l(:,:),y(:)
      integer :: n
      n=size(b); allocate(l(n,n),y(n))
      call chol_lower(a,l,info); if (info /= 0) return
      call solve_lower(l,b,y)
      call solve_upper(transpose(l),y,x)
   end subroutine solve_spd

   subroutine inv_spd(a,ainv,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(size(a,1),size(a,2))
      integer, intent(out) :: info
      integer :: i,n
      real(dp), allocatable :: e(:),x(:)
      n=size(a,1); allocate(e(n),x(n)); ainv=0.0_dp
      do i=1,n
         e=0.0_dp; e(i)=1.0_dp
         call solve_spd(a,e,x,info); if (info /= 0) return
         ainv(:,i)=x
      end do
      call symmetrize(ainv)
   end subroutine inv_spd

   pure function outer_product(x,y) result(a)
      real(dp), intent(in) :: x(:),y(:)
      real(dp) :: a(size(x),size(y))
      integer :: i,j
      do j=1,size(y); do i=1,size(x); a(i,j)=x(i)*y(j); end do; end do
   end function outer_product

   pure real(dp) function trace_mat(a) result(t)
      real(dp), intent(in) :: a(:,:)
      integer :: i
      t=0.0_dp; do i=1,min(size(a,1),size(a,2)); t=t+a(i,i); end do
   end function trace_mat

   subroutine logdet_spd(a,value,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: value
      integer, intent(out) :: info
      real(dp), allocatable :: l(:,:)
      integer :: i,n
      n=size(a,1); allocate(l(n,n)); call chol_lower(a,l,info)
      if (info /= 0) then; value=-huge(1.0_dp); return; end if
      value=0.0_dp; do i=1,n; value=value+2.0_dp*log(l(i,i)); end do
   end subroutine logdet_spd

   subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:,:)
      integer :: i,j,n
      n=min(size(a,1),size(a,2))
      do j=1,n; do i=j+1,n
         a(i,j)=0.5_dp*(a(i,j)+a(j,i)); a(j,i)=a(i,j)
      end do; end do
   end subroutine symmetrize

   subroutine least_squares_normal(x,y,beta,sigma2,vcov,info)
      real(dp), intent(in) :: x(:,:),y(:)
      real(dp), intent(out) :: beta(size(x,2)),sigma2,vcov(size(x,2),size(x,2))
      integer, intent(out) :: info
      real(dp), allocatable :: xtx(:,:),xty(:),r(:),inv(:,:)
      integer :: n,k
      n=size(x,1); k=size(x,2)
      allocate(xtx(k,k),xty(k),r(n),inv(k,k))
      xtx=matmul(transpose(x),x); xty=matmul(transpose(x),y)
      call solve_spd(xtx,xty,beta,info); if (info /= 0) return
      r=y-matmul(x,beta); sigma2=dot_product(r,r)/real(max(1,n-k),dp)
      call inv_spd(xtx,inv,info); if (info /= 0) return
      vcov=sigma2*inv
   end subroutine least_squares_normal

   subroutine covariance_matrix(x,cov)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: cov(size(x,2),size(x,2))
      real(dp) :: mu(size(x,2)),d(size(x,2))
      integer :: i,n
      n=size(x,1); mu=sum(x,dim=1)/real(n,dp); cov=0.0_dp
      do i=1,n
         d=x(i,:)-mu; cov=cov+outer_product(d,d)
      end do
      cov=cov/real(max(1,n-1),dp)
   end subroutine covariance_matrix
end module mcmcpack_linalg
