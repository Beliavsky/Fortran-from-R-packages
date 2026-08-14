! SPDX-License-Identifier: GPL-2.0-or-later
module lsei_nnls
   use lsei_kinds, only : dp
   use lsei_types
   use lsei_linalg, only : pseudoinverse, matrix_rank
   implicit none
   private
   public :: nnls_solve, pnnls_solve, ldp_solve
contains
   subroutine nnls_solve(a,b,res,max_iter)
      real(dp), intent(in) :: a(:,:), b(:)
      type(ls_result), intent(out) :: res
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: x(:),z(:),w(:),r(:),ap(:,:),pinv(:,:),zp(:)
      logical, allocatable :: passive(:)
      integer, allocatable :: pidx(:)
      integer :: m,n,j,k,t,np,it,lim,rankp,st
      real(dp) :: best,alpha,eps
      m=size(a,1); n=size(a,2)
      allocate(res%x(max(0,n)),res%residuals(max(0,m)),res%dual(max(0,n)),res%index(max(0,n)))
      res%x=0.0_dp; res%residuals=0.0_dp; res%dual=0.0_dp; res%index=0
      if(m<=0 .or. n<=0 .or. size(b)/=m) then; res%mode=LSEI_BAD_DIMENSIONS; return; end if
      lim=3*n; if(present(max_iter)) lim=max_iter
      eps=100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)),maxval(abs(b)))
      allocate(x(n),z(n),w(n),r(m),passive(n)); x=0.0_dp; passive=.false.; it=0
      do
         r=b-matmul(a,x); w=matmul(transpose(a),r); t=0; best=eps
         do j=1,n
            if(.not.passive(j) .and. w(j)>best) then; best=w(j); t=j; end if
         end do
         if(t==0) exit
         passive(t)=.true.
         do
            np=count(passive); allocate(pidx(np),ap(m,np),pinv(np,m),zp(np)); k=0
            do j=1,n
               if(passive(j)) then; k=k+1; pidx(k)=j; ap(:,k)=a(:,j); end if
            end do
            call pseudoinverse(ap,pinv,rankp,1.0e-14_dp,st); zp=matmul(pinv,b); z=0.0_dp
            do k=1,np; z(pidx(k))=zp(k); end do
            deallocate(ap,pinv,zp,pidx)
            if(all(pack(z,passive)>eps)) then; x=z; exit; end if
            alpha=huge(1.0_dp)
            do j=1,n
               if(passive(j) .and. z(j)<=eps .and. x(j)>z(j)) alpha=min(alpha,x(j)/(x(j)-z(j)))
            end do
            if(alpha>=huge(1.0_dp)/2.0_dp) alpha=0.0_dp
            x=x+alpha*(z-x)
            do j=1,n
               if(passive(j) .and. x(j)<=eps) then; x(j)=0.0_dp; passive(j)=.false.; end if
            end do
            it=it+1
            if(it>lim) exit
         end do
         if(it>lim) exit
      end do
      res%x=x; res%residuals=b-matmul(a,x); res%rnorm=norm2(res%residuals)
      res%objective=res%rnorm**2; res%dual=matmul(transpose(a),res%residuals); res%iterations=it
      k=0
      do j=1,n; if(passive(j)) then; k=k+1; res%index(k)=j; end if; end do
      do j=1,n; if(.not.passive(j)) then; k=k+1; res%index(k)=j; end if; end do
      res%rank=count(passive); res%k=0
      if(it>lim) then; res%mode=LSEI_ITERATION_LIMIT; else; res%mode=LSEI_SUCCESS; end if
   end subroutine nnls_solve

   recursive subroutine pnnls_solve(a,b,kfree,res,sum_value,max_iter)
      real(dp), intent(in) :: a(:,:), b(:)
      integer, intent(in) :: kfree
      type(ls_result), intent(out) :: res
      real(dp), intent(in), optional :: sum_value
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: aw(:,:),bw(:),af(:,:),ac(:,:),afp(:,:),pperp(:,:),ared(:,:),bred(:),xf(:)
      type(ls_result) :: nr
      integer :: m,n,k0,nc,rf,st,j
      real(dp) :: s,t
      m=size(a,1); n=size(a,2); k0=max(0,min(kfree,n))
      if(present(sum_value)) then
         if(sum_value<=0.0_dp .or. k0>=n) then
            allocate(res%x(max(0,n)),res%residuals(max(0,m)),res%dual(max(0,n)),res%index(max(0,n)))
            res%x=0.0_dp; res%residuals=0.0_dp; res%dual=0.0_dp; res%index=0; res%mode=LSEI_BAD_DIMENSIONS; return
         end if
         s=sum_value; allocate(aw(m+1,n),bw(m+1)); aw(1:m,:)=a; bw(1:m)=b
         do j=k0+1,n; aw(1:m,j)=a(:,j)*s-b; end do
         aw(m+1,:)=0.0_dp; aw(m+1,k0+1:n)=1.0_dp; bw(m+1)=1.0_dp
         if(present(max_iter)) then; call pnnls_solve(aw,bw,k0,nr,max_iter=max_iter); else; call pnnls_solve(aw,bw,k0,nr); end if
         res=nr; t=sum(res%x(k0+1:n))
         if(abs(t)>tiny(1.0_dp)) then
            res%x=res%x/t; res%x(k0+1:n)=res%x(k0+1:n)*s
            res%rnorm=sqrt(max((res%rnorm/t)**2-(1.0_dp-1.0_dp/t)**2,0.0_dp))
            res%residuals=b-matmul(a,res%x); res%objective=res%rnorm**2
         end if
         res%k=k0; return
      end if
      nc=n-k0
      if(k0==0) then
         if(present(max_iter)) then; call nnls_solve(a,b,res,max_iter); else; call nnls_solve(a,b,res); end if
         return
      end if
      allocate(af(m,k0),afp(k0,m),pperp(m,m)); af=a(:,1:k0)
      call pseudoinverse(af,afp,rf,1.0e-14_dp,st); pperp=0.0_dp
      do j=1,m; pperp(j,j)=1.0_dp; end do
      pperp=pperp-matmul(af,afp)
      if(nc>0) then
         allocate(ac(m,nc),ared(m,nc),bred(m)); ac=a(:,k0+1:n); ared=matmul(pperp,ac); bred=matmul(pperp,b)
         if(present(max_iter)) then; call nnls_solve(ared,bred,nr,max_iter); else; call nnls_solve(ared,bred,nr); end if
         allocate(xf(k0)); xf=matmul(afp,b-matmul(ac,nr%x))
      else
         allocate(xf(k0)); xf=matmul(afp,b)
      end if
      allocate(res%x(n),res%residuals(m),res%dual(n),res%index(n)); res%x=0.0_dp; res%x(1:k0)=xf
      if(nc>0) res%x(k0+1:n)=nr%x
      res%residuals=b-matmul(a,res%x); res%rnorm=norm2(res%residuals); res%objective=res%rnorm**2
      res%dual=matmul(transpose(a),res%residuals); res%mode=merge(nr%mode,LSEI_SUCCESS,nc>0); res%k=rf; res%rank=rf
      res%index=[(j,j=1,n)]; if(nc>0) res%iterations=nr%iterations
   end subroutine pnnls_solve

   subroutine ldp_solve(e,f,res,tol)
      real(dp), intent(in) :: e(:,:),f(:)
      type(ls_result), intent(out) :: res
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: et(:,:),h(:),r(:)
      type(ls_result) :: nr
      real(dp) :: eps,den
      integer :: m,n
      m=size(e,1); n=size(e,2); eps=1.0e-15_dp; if(present(tol)) eps=tol
      allocate(et(n+1,m),h(n+1)); et(1:n,:)=transpose(e); et(n+1,:)=f; h=0.0_dp; h(n+1)=1.0_dp
      call nnls_solve(et,h,nr); allocate(r(n+1)); r=matmul(et,nr%x)-h; den=r(n+1)
      allocate(res%x(n),res%residuals(m),res%dual(n),res%index(n)); res%x=0.0_dp; res%residuals=0.0_dp
      res%dual=0.0_dp; res%index=[(m,m=1,n)]
      if(norm2(r)<=eps .or. abs(den)<=eps) then; res%mode=LSEI_INFEASIBLE; return; end if
      res%x=-r(1:n)/den; res%residuals=matmul(e,res%x)-f; res%rnorm=norm2(res%x); res%objective=res%rnorm**2
      if(minval(res%residuals)<-100.0_dp*eps) then; res%mode=LSEI_INFEASIBLE; else; res%mode=LSEI_SUCCESS; end if
   end subroutine ldp_solve
end module lsei_nnls
