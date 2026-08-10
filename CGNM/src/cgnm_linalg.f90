! SPDX-License-Identifier: MIT
module cgnm_linalg
   use cgnm_kinds, only : dp
   implicit none
   private
   public :: solve_sym_pinv, cgnr_reg, cgnr_reg_weight, weighted_local_linear
contains

   subroutine jacobi_sym(a, eval, evec, ierr)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: eval(size(a,1)), evec(size(a,1),size(a,1))
      integer, intent(out) :: ierr
      real(dp), allocatable :: d(:,:)
      real(dp) :: app, aqq, apq, tau, t, c, s, dkp, dkq, vkp, vkq
      real(dp) :: off, scale
      integer :: n, i, j, p, q, k, iter, maxit
      n = size(a,1); ierr = 0
      if (size(a,2) /= n) then
         ierr = 1; return
      end if
      d = 0.5_dp*(a+transpose(a))
      evec = 0.0_dp
      do i = 1, n
         evec(i,i) = 1.0_dp
      end do
      maxit = max(50,50*n*n)
      do iter = 1, maxit
         off = 0.0_dp; p = 1; q = min(2,n)
         do i = 1, n-1
            do j = i+1, n
               if (abs(d(i,j)) > off) then
                  off = abs(d(i,j)); p=i; q=j
               end if
            end do
         end do
         scale = max(1.0_dp,maxval(abs(d)))
         if (off <= 10.0_dp*epsilon(1.0_dp)*scale) exit
         apq=d(p,q); app=d(p,p); aqq=d(q,q)
         if (abs(apq) <= tiny(1.0_dp)) cycle
         tau=(aqq-app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t=1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
         else
            t=-1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
         end if
         c=1.0_dp/sqrt(1.0_dp+t*t); s=t*c
         do k=1,n
            if (k/=p .and. k/=q) then
               dkp=d(k,p); dkq=d(k,q)
               d(k,p)=c*dkp-s*dkq; d(p,k)=d(k,p)
               d(k,q)=s*dkp+c*dkq; d(q,k)=d(k,q)
            end if
         end do
         d(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         d(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
         d(p,q)=0.0_dp; d(q,p)=0.0_dp
         do k=1,n
            vkp=evec(k,p); vkq=evec(k,q)
            evec(k,p)=c*vkp-s*vkq
            evec(k,q)=s*vkp+c*vkq
         end do
      end do
      if (iter > maxit) ierr=2
      do i=1,n
         eval(i)=d(i,i)
      end do
   end subroutine jacobi_sym

   subroutine solve_sym_pinv(a,b,x,ierr,rtol)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      integer, intent(out) :: ierr
      real(dp), intent(in), optional :: rtol
      real(dp), allocatable :: eval(:), evec(:,:), tmp(:)
      real(dp) :: tol, mx
      integer :: n, i
      n=size(b); x=0.0_dp
      allocate(eval(n),evec(n,n),tmp(n))
      call jacobi_sym(a,eval,evec,ierr)
      if (ierr/=0) return
      mx=max(1.0_dp,maxval(abs(eval)))
      tol=1.0e-12_dp*mx
      if (present(rtol)) tol=rtol*mx
      tmp=matmul(transpose(evec),b)
      do i=1,n
         if (abs(eval(i))>tol) tmp(i)=tmp(i)/eval(i)
         if (abs(eval(i))<=tol) tmp(i)=0.0_dp
      end do
      x=matmul(evec,tmp)
   end subroutine solve_sym_pinv

   subroutine cgnr_reg(a,b,lambda,x)
      real(dp), intent(in) :: a(:,:), b(:), lambda
      real(dp), intent(out) :: x(size(a,2))
      real(dp), allocatable :: r(:), p(:), ap(:), gram_p(:)
      real(dp) :: rr_old, rr_new, den, alpha, beta
      integer :: n, i, maxit
      n=size(a,2); maxit=size(a,1)
      allocate(r(n),p(n),ap(size(a,1)),gram_p(n))
      x=0.0_dp
      r=matmul(transpose(a),b); p=r; rr_old=dot_product(r,r)
      do i=1,maxit
         ap=matmul(a,p)
         gram_p=matmul(transpose(a),ap)+lambda*p
         den=dot_product(p,gram_p)
         if (abs(den)<=tiny(1.0_dp)) exit
         if (rr_old<=tiny(1.0_dp)) exit
         alpha=rr_old/den
         x=x+alpha*p; r=r-alpha*gram_p
         rr_new=dot_product(r,r)
         if (rr_new<=tiny(1.0_dp)) exit
         beta=rr_new/rr_old; p=r+beta*p; rr_old=rr_new
      end do
   end subroutine cgnr_reg

   subroutine cgnr_reg_weight(a,b,lambda,weight,x)
      real(dp), intent(in) :: a(:,:), b(:), lambda, weight(:)
      real(dp), intent(out) :: x(size(a,2))
      real(dp), allocatable :: aw(:,:), bw(:)
      integer :: i
      allocate(aw(size(a,1),size(a,2)),bw(size(b)))
      aw=a; bw=b
      do i=1,min(size(weight),size(a,1))
         aw(i,:)=weight(i)*aw(i,:); bw(i)=weight(i)*bw(i)
      end do
      call cgnr_reg(aw,bw,lambda,x)
   end subroutine cgnr_reg_weight

   subroutine weighted_local_linear(dx,dy,weight,aout,ierr)
      ! Minimize || diag(weight) * (dx*a - dy) || for all dy columns.
      real(dp), intent(in) :: dx(:,:), dy(:,:), weight(:)
      real(dp), intent(out) :: aout(size(dy,2),size(dx,2))
      integer, intent(out) :: ierr
      real(dp), allocatable :: gram(:,:), rhs(:), sol(:)
      integer :: i,j,k,p,m
      p=size(dx,2); m=size(dx,1)
      allocate(gram(p,p),rhs(p),sol(p))
      gram=0.0_dp
      do i=1,m
         do j=1,p
            do k=1,p
               gram(j,k)=gram(j,k)+(weight(i)*dx(i,j))*(weight(i)*dx(i,k))
            end do
         end do
      end do
      ierr=0
      do k=1,size(dy,2)
         rhs=0.0_dp
         do i=1,m
            rhs=rhs+(weight(i)*weight(i)*dy(i,k))*dx(i,:)
         end do
         call solve_sym_pinv(gram,rhs,sol,ierr)
         if (ierr/=0) return
         aout(k,:)=sol
      end do
   end subroutine weighted_local_linear
end module cgnm_linalg
