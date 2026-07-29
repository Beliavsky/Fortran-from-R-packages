! SPDX-License-Identifier: GPL-3.0-only
module linalg_mod
   use kind_mod, only: dp
   implicit none
   private
   public :: symmetric_pseudoinverse
contains
   pure subroutine symmetric_pseudoinverse(matrix, inverse, info)
      !! Jacobi eigensolver-based Moore-Penrose inverse for small symmetric matrices.
      real(dp), intent(in) :: matrix(:,:)
      real(dp), intent(out) :: inverse(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: a(:,:), v(:,:), eig(:)
      real(dp) :: app,aqq,apq,tau,t,c,s,maxeig,tol,temp
      integer :: n,p,q,i,j,sweep
      n=size(matrix,1); inverse=0.0_dp; info=0
      if(size(matrix,2)/=n .or. size(inverse,1)/=n .or. size(inverse,2)/=n)then
         info=1;return
      end if
      allocate(a(n,n),v(n,n),eig(n)); a=0.5_dp*(matrix+transpose(matrix));v=0.0_dp
      do i=1,n;v(i,i)=1.0_dp;end do
      do sweep=1,max(30,20*n*n)
         p=1;q=1;apq=0.0_dp
         do j=2,n;do i=1,j-1
            if(abs(a(i,j))>abs(apq))then;apq=a(i,j);p=i;q=j;end if
         end do;end do
         if(abs(apq)<=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a))))exit
         app=a(p,p);aqq=a(q,q);tau=(aqq-app)/(2.0_dp*apq)
         t=sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau));c=1.0_dp/sqrt(1.0_dp+t*t);s=t*c
         do j=1,n
            if(j/=p .and. j/=q)then
               temp=a(j,p);a(j,p)=c*temp-s*a(j,q);a(p,j)=a(j,p)
               a(j,q)=s*temp+c*a(j,q);a(q,j)=a(j,q)
            end if
         end do
         a(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         a(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
         a(p,q)=0.0_dp;a(q,p)=0.0_dp
         do j=1,n
            temp=v(j,p);v(j,p)=c*temp-s*v(j,q);v(j,q)=s*temp+c*v(j,q)
         end do
      end do
      do i=1,n;eig(i)=a(i,i);end do
      maxeig=maxval(abs(eig));tol=max(1.0e-12_dp,real(n,dp)*epsilon(1.0_dp))*max(1.0_dp,maxeig)
      do i=1,n
         if(abs(eig(i))>tol)inverse=inverse+outer(v(:,i),v(:,i))/eig(i)
      end do
      if(any(.not.(abs(inverse)<=huge(1.0_dp))))info=2
   contains
      pure function outer(x,y) result(z)
         real(dp),intent(in)::x(:),y(:)
         real(dp)::z(size(x),size(y))
         integer::ii,jj
         do jj=1,size(y);do ii=1,size(x);z(ii,jj)=x(ii)*y(jj);end do;end do
      end function outer
   end subroutine symmetric_pseudoinverse
end module linalg_mod
