! SPDX-License-Identifier: GPL-2.0-only
module ks_linalg
   use ks_kinds, only: dp
   use ks_lapack, only: dpotrf, dpotri, dpotrs, dsyev
   implicit none
   private
   public :: covariance_matrix, spd_inverse, spd_logdet, spd_cholesky
   public :: spd_solve, matrix_sqrt, matrix_power_sym, determinant_spd
   public :: is_symmetric, identity_matrix, trace_matrix, symmetric_eigen
contains
   function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i=1,n
         a(i,i)=1.0_dp
      end do
   end function identity_matrix

   function trace_matrix(a) result(t)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: t
      integer :: i
      t=0.0_dp
      do i=1,min(size(a,1),size(a,2))
         t=t+a(i,i)
      end do
   end function trace_matrix

   logical function is_symmetric(a, tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: tol
      real(dp) :: eps
      eps = 100.0_dp*epsilon(1.0_dp)
      if (present(tol)) eps=tol
      if (size(a,1)/=size(a,2)) then
         is_symmetric=.false.
      else
         is_symmetric=maxval(abs(a-transpose(a))) <= eps*max(1.0_dp,maxval(abs(a)))
      end if
   end function is_symmetric

   subroutine covariance_matrix(x, s, center)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: s(size(x,2),size(x,2))
      real(dp), intent(out), optional :: center(size(x,2))
      real(dp) :: mu(size(x,2)), z(size(x,1),size(x,2))
      integer :: n
      n=size(x,1)
      mu=sum(x,dim=1)/real(max(n,1),dp)
      if (present(center)) center=mu
      if (n<=1) then
         s=0.0_dp
         return
      end if
      z=x-spread(mu,1,n)
      s=matmul(transpose(z),z)/real(n-1,dp)
      s=0.5_dp*(s+transpose(s))
   end subroutine covariance_matrix

   subroutine spd_cholesky(a,l,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: l(size(a,1),size(a,2))
      integer, intent(out), optional :: info
      integer :: ierr,n,i,j
      n=size(a,1)
      l=a
      call dpotrf('L',n,l,n,ierr)
      if (ierr==0) then
         do j=1,n
            do i=1,j-1
               l(i,j)=0.0_dp
            end do
         end do
      end if
      if (present(info)) info=ierr
   end subroutine spd_cholesky

   subroutine spd_inverse(a,ainv,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(size(a,1),size(a,2))
      integer, intent(out), optional :: info
      integer :: ierr,n,i,j
      n=size(a,1)
      ainv=a
      call dpotrf('L',n,ainv,n,ierr)
      if (ierr==0) call dpotri('L',n,ainv,n,ierr)
      if (ierr==0) then
         do j=1,n
            do i=1,j-1
               ainv(i,j)=ainv(j,i)
            end do
         end do
      end if
      if (present(info)) info=ierr
   end subroutine spd_inverse

   function spd_logdet(a,info) result(v)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: info
      real(dp) :: v
      real(dp) :: l(size(a,1),size(a,2))
      integer :: ierr,i
      call spd_cholesky(a,l,ierr)
      if (ierr/=0) then
         v=huge(1.0_dp)
      else
         v=0.0_dp
         do i=1,size(a,1)
            v=v+2.0_dp*log(l(i,i))
         end do
      end if
      if (present(info)) info=ierr
   end function spd_logdet

   function determinant_spd(a,info) result(v)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: info
      real(dp) :: v
      integer :: ierr
      v=spd_logdet(a,ierr)
      if (ierr==0) then
         v=exp(v)
      else
         v=0.0_dp
      end if
      if (present(info)) info=ierr
   end function determinant_spd

   subroutine spd_solve(a,b,x,info)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), intent(out) :: x(size(b,1),size(b,2))
      integer, intent(out), optional :: info
      real(dp) :: l(size(a,1),size(a,2))
      integer :: ierr,n
      n=size(a,1)
      l=a
      x=b
      call dpotrf('L',n,l,n,ierr)
      if (ierr==0) call dpotrs('L',n,size(b,2),l,n,x,n,ierr)
      if (present(info)) info=ierr
   end subroutine spd_solve

   subroutine symmetric_eigen(a,eigenvalues,eigenvectors,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: eigenvalues(size(a,1)), eigenvectors(size(a,1),size(a,2))
      integer, intent(out), optional :: info
      real(dp), allocatable :: work(:)
      integer :: n,lwork,ierr
      n=size(a,1); eigenvectors=0.5_dp*(a+transpose(a)); lwork=max(1,3*n-1); allocate(work(lwork))
      call dsyev('V','U',n,eigenvectors,n,eigenvalues,work,lwork,ierr)
      if(present(info)) info=ierr
   end subroutine symmetric_eigen

   subroutine matrix_sqrt(a,root,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: root(size(a,1),size(a,2))
      integer, intent(out), optional :: info
      integer :: n,lwork,ierr,i
      real(dp), allocatable :: eig(:),work(:),v(:,:)
      n=size(a,1)
      allocate(eig(n),v(n,n))
      v=0.5_dp*(a+transpose(a))
      lwork=max(1,3*n-1)
      allocate(work(lwork))
      call dsyev('V','U',n,v,n,eig,work,lwork,ierr)
      if (ierr==0 .and. minval(eig)>=-100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(eig)))) then
         eig=max(eig,0.0_dp)
         do i=1,n
            v(:,i)=v(:,i)*sqrt(eig(i))
         end do
         root=matmul(v,transpose(v))
      else
         root=0.0_dp
         if (ierr==0) ierr=-1
      end if
      if (present(info)) info=ierr
   end subroutine matrix_sqrt

   subroutine matrix_power_sym(a,power,b,info)
      real(dp), intent(in) :: a(:,:), power
      real(dp), intent(out) :: b(size(a,1),size(a,2))
      integer, intent(out), optional :: info
      integer :: n,lwork,ierr,i
      real(dp), allocatable :: eig(:),work(:),v(:,:)
      n=size(a,1)
      allocate(eig(n),v(n,n),work(max(1,3*n-1)))
      v=0.5_dp*(a+transpose(a))
      lwork=size(work)
      call dsyev('V','U',n,v,n,eig,work,lwork,ierr)
      if (ierr==0) then
         if (power<0.0_dp .and. minval(eig)<=0.0_dp) then
            ierr=-1
            b=0.0_dp
         else
            do i=1,n
               if (eig(i)<0.0_dp .and. abs(power-nint(power))>100*epsilon(1.0_dp)) then
                  ierr=-2; b=0.0_dp; exit
               end if
               v(:,i)=v(:,i)*(eig(i)**power)
            end do
            if (ierr==0) b=matmul(v,transpose(v))
         end if
      else
         b=0.0_dp
      end if
      if (present(info)) info=ierr
   end subroutine matrix_power_sym
end module ks_linalg
