! Matrix utilities corresponding to distr's PosSemDefSymmMatrix sqrt/solve methods.
! Copyright (C) 2005-2025 distr authors.
! SPDX-License-Identifier: LGPL-3.0-only
module distr_matrix
   use distr_kinds, only : dp, eps_dp
   implicit none
   private
   public :: symmetric_psd_sqrt, symmetric_pseudoinverse, solve_linear
contains
   subroutine jacobi_eigen(a,eval,evec)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: eval(size(a,1)),evec(size(a,1),size(a,1))
      real(dp), allocatable :: b(:,:)
      real(dp) :: app,aqq,apq,tau,t,c,s,tmp,bkp,bkq,maxoff
      integer :: n,p,q,k,it
      n=size(a,1); if (size(a,2)/=n) error stop 'jacobi_eigen requires square matrix'
      allocate(b(n,n)); b=a; evec=0.0_dp; do k=1,n; evec(k,k)=1.0_dp; end do
      do it=1,100*n*n
         maxoff=0.0_dp; p=1; q=min(2,n)
         do k=1,n-1
            do q=k+1,n
               if (abs(b(k,q))>maxoff) then; maxoff=abs(b(k,q)); p=k; end if
            end do
         end do
         if (maxoff<=64.0_dp*eps_dp*max(1.0_dp,maxval(abs(b)))) exit
         q=p+1; maxoff=0.0_dp
         do k=p+1,n; if (abs(b(p,k))>maxoff) then; maxoff=abs(b(p,k)); q=k; end if; end do
         if (maxoff==0.0_dp) cycle
         app=b(p,p); aqq=b(q,q); apq=b(p,q); tau=(aqq-app)/(2.0_dp*apq)
         if (tau>=0.0_dp) then; t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau)); else; t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau)); end if
         c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
         do k=1,n
            if (k/=p.and.k/=q) then
               bkp=b(k,p); bkq=b(k,q); b(k,p)=c*bkp-s*bkq; b(p,k)=b(k,p); b(k,q)=s*bkp+c*bkq; b(q,k)=b(k,q)
            end if
         end do
         b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq; b(p,q)=0.0_dp; b(q,p)=0.0_dp
         do k=1,n
            tmp=evec(k,p); evec(k,p)=c*tmp-s*evec(k,q); evec(k,q)=s*tmp+c*evec(k,q)
         end do
      end do
      do k=1,n; eval(k)=b(k,k); end do
   end subroutine jacobi_eigen

   subroutine symmetric_psd_sqrt(a,sqrta,tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: sqrta(size(a,1),size(a,2))
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: val(:),vec(:,:)
      real(dp) :: tv
      integer :: n,i,j,k
      n=size(a,1); allocate(val(n),vec(n,n)); call jacobi_eigen(a,val,vec)
      tv=64.0_dp*eps_dp*max(1.0_dp,maxval(abs(val))); if (present(tol)) tv=tol
      where(val<0.0_dp.and.abs(val)<=tv) val=0.0_dp
      if (any(val<0.0_dp)) error stop 'matrix is not positive semidefinite'
      sqrta=0.0_dp
      do k=1,n; do i=1,n; do j=1,n; sqrta(i,j)=sqrta(i,j)+vec(i,k)*sqrt(val(k))*vec(j,k); end do; end do; end do
   end subroutine symmetric_psd_sqrt

   subroutine symmetric_pseudoinverse(a,ainv,tol)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: ainv(size(a,1),size(a,2))
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: val(:),vec(:,:)
      real(dp) :: tv,w
      integer :: n,i,j,k
      n=size(a,1); allocate(val(n),vec(n,n)); call jacobi_eigen(a,val,vec)
      tv=sqrt(eps_dp)*max(1.0_dp,maxval(abs(val))); if (present(tol)) tv=tol
      ainv=0.0_dp
      do k=1,n
         if (abs(val(k))>tv) then
            w=1.0_dp/val(k); do i=1,n; do j=1,n; ainv(i,j)=ainv(i,j)+vec(i,k)*w*vec(j,k); end do; end do
         end if
      end do
   end subroutine symmetric_pseudoinverse

   subroutine solve_linear(a,b,x,tol)
      real(dp), intent(in) :: a(:,:),b(:)
      real(dp), intent(out) :: x(size(b))
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: aug(:,:)
      real(dp) :: tv,piv,tmp,f
      integer :: n,i,j,k,p
      n=size(a,1); if (size(a,2)/=n.or.size(b)/=n) error stop 'nonconformable solve'
      tv=sqrt(eps_dp); if (present(tol)) tv=tol
      allocate(aug(n,n+1)); aug(:,:n)=a; aug(:,n+1)=b
      do k=1,n
         p=k; do i=k+1,n; if (abs(aug(i,k))>abs(aug(p,k))) p=i; end do
         if (abs(aug(p,k))<=tv) error stop 'singular matrix in solve_linear'
         if (p/=k) then; do j=k,n+1; tmp=aug(k,j); aug(k,j)=aug(p,j); aug(p,j)=tmp; end do; end if
         piv=aug(k,k); aug(k,k:n+1)=aug(k,k:n+1)/piv
         do i=1,n
            if (i==k) cycle; f=aug(i,k); aug(i,k:n+1)=aug(i,k:n+1)-f*aug(k,k:n+1)
         end do
      end do
      x=aug(:,n+1)
   end subroutine solve_linear
end module distr_matrix
