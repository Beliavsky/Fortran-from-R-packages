! SPDX-License-Identifier: GPL-2.0-only
module nlsic_nnls
   use nlsic_kinds, only : dp
   use nlsic_linalg, only : pseudoinverse, vecnorm
   implicit none
   private
   public :: nnls_solve, ldp_solve
contains
   subroutine nnls_solve(a,b,x,rnorm,status,max_iter)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      real(dp), intent(out) :: rnorm
      integer, intent(out) :: status
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: z(:),w(:),r(:),ap(:,:),pinv(:,:),zp(:)
      logical, allocatable :: passive(:)
      integer, allocatable :: pidx(:)
      integer :: m,n,j,k,t,np,it,lim,rankp,st
      real(dp) :: best,alpha,eps,scale
      m=size(a,1); n=size(a,2); x=0.0_dp; rnorm=huge(1.0_dp); status=0
      if(m<=0 .or. n<=0 .or. size(b)/=m .or. size(x)/=n) then
         status=-1; return
      end if
      lim=3*n; if(present(max_iter)) lim=max_iter
      scale=max(1.0_dp,maxval(abs(a)),maxval(abs(b)))
      eps=100.0_dp*epsilon(1.0_dp)*scale
      allocate(z(n),w(n),r(m),passive(n)); passive=.false.; it=0
      do
         r=b-matmul(a,x); w=matmul(transpose(a),r); t=0; best=eps
         do j=1,n
            if(.not.passive(j) .and. w(j)>best) then
               best=w(j); t=j
            end if
         end do
         if(t==0) exit
         passive(t)=.true.
         do
            np=count(passive); allocate(pidx(np),ap(m,np),pinv(np,m),zp(np)); k=0
            do j=1,n
               if(passive(j)) then
                  k=k+1; pidx(k)=j; ap(:,k)=a(:,j)
               end if
            end do
            call pseudoinverse(ap,pinv,rankp,1.0e14_dp,st)
            zp=matmul(pinv,b); z=0.0_dp
            do k=1,np
               z(pidx(k))=zp(k)
            end do
            deallocate(ap,pinv,zp,pidx)
            if(all(pack(z,passive)>eps)) then
               x=z; exit
            end if
            alpha=huge(1.0_dp)
            do j=1,n
               if(passive(j) .and. z(j)<=eps .and. x(j)>z(j)) then
                  alpha=min(alpha,x(j)/(x(j)-z(j)))
               end if
            end do
            if(alpha>=huge(1.0_dp)/2.0_dp) alpha=0.0_dp
            x=x+alpha*(z-x)
            do j=1,n
               if(passive(j) .and. x(j)<=eps) then
                  x(j)=0.0_dp; passive(j)=.false.
               end if
            end do
            it=it+1
            if(it>lim) exit
         end do
         if(it>lim) exit
      end do
      r=b-matmul(a,x); rnorm=vecnorm(r)
      if(it>lim) status=1
   end subroutine nnls_solve

   subroutine ldp_solve(u,co,x,feasible,rcond)
      real(dp), intent(in) :: u(:,:),co(:)
      real(dp), intent(out) :: x(:)
      logical, intent(out) :: feasible
      real(dp), intent(in), optional :: rcond
      real(dp), allocatable :: e(:,:),f(:),z(:),r(:),uu(:,:),cc(:)
      logical, allocatable :: keep(:)
      real(dp) :: tol,rn,den,maxco,alpha
      integer :: m,n,i,nkeep,status
      m=size(u,1); n=size(u,2); x=0.0_dp; feasible=.false.
      if(size(co)/=m .or. size(x)/=n) return
      tol=1.0e-10_dp; if(present(rcond)) then; if(rcond>1.0_dp) tol=1.0_dp/rcond; end if
      if(m==0) then; feasible=.true.; return; end if
      allocate(keep(m)); keep=.true.
      do i=1,m
         if(maxval(abs(u(i,:)))<tol) then
            if(co(i)>100.0_dp*tol) return
            keep(i)=.false.
         end if
      end do
      nkeep=count(keep)
      if(nkeep==0) then; feasible=.true.; return; end if
      allocate(uu(nkeep,n),cc(nkeep)); uu=reshape(pack(u,spread(keep,2,n)),[nkeep,n]); cc=pack(co,keep)
      maxco=maxval(cc)
      if(maxco<=100.0_dp*tol) then; feasible=.true.; return; end if
      allocate(e(n+1,nkeep),f(n+1),z(nkeep),r(n+1))
      e(1:n,:)=transpose(uu); e(n+1,:)=cc; f=0.0_dp; f(n+1)=1.0_dp
      call nnls_solve(e,f,z,rn,status)
      r=matmul(e,z)-f; den=r(n+1)
      if(rn<=tol .or. abs(den)<=tol) return
      x=-r(1:n)/den
      if(any(matmul(uu,x)-cc < -100.0_dp*tol)) return
      if(all(matmul(uu,x)-cc>0.0_dp) .and. vecnorm(x)>tol) then
         alpha=huge(1.0_dp)
         do i=1,nkeep
            if(cc(i)>tol) then
               if(dot_product(uu(i,:),x)/cc(i)>=1.0_dp) then
                  alpha=min(alpha,dot_product(uu(i,:),x)/cc(i))
               end if
            end if
         end do
         if(alpha<huge(1.0_dp)/2.0_dp .and. alpha>1.0_dp) x=x/alpha
      end if
      feasible=.true.
   end subroutine ldp_solve
end module nlsic_nnls
