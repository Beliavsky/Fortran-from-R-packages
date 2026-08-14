! SPDX-License-Identifier: GPL-2.0-only
module ks_optimize
   use ks_kinds, only: dp
   implicit none
   private
   public :: golden_minimize, nelder_mead
   abstract interface
      function scalar_objective(x) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function scalar_objective
      function one_objective(x) result(f)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: f
      end function one_objective
   end interface
contains
   subroutine golden_minimize(f,a,b,xmin,fmin,tol,maxiter)
      procedure(one_objective) :: f
      real(dp),intent(in)::a,b
      real(dp),intent(out)::xmin,fmin
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      real(dp)::lo,hi,c,d,fc,fd,eps,gr
      integer::it,nit
      eps=sqrt(epsilon(1.0_dp));if(present(tol))eps=tol
      nit=500;if(present(maxiter))nit=maxiter
      gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp;lo=a;hi=b
      c=hi-gr*(hi-lo);d=lo+gr*(hi-lo);fc=f(c);fd=f(d)
      do it=1,nit
         if(abs(hi-lo)<=eps*max(1.0_dp,abs(0.5_dp*(lo+hi))))exit
         if(fc<fd)then
            hi=d;d=c;fd=fc;c=hi-gr*(hi-lo);fc=f(c)
         else
            lo=c;c=d;fc=fd;d=lo+gr*(hi-lo);fd=f(d)
         end if
      end do
      if(fc<fd)then;xmin=c;fmin=fc;else;xmin=d;fmin=fd;end if
   end subroutine golden_minimize

   subroutine sort_simplex(x,fv)
      real(dp),intent(inout)::x(:,:),fv(:)
      integer::i,j,k,n
      real(dp)::tf
      real(dp),allocatable::tx(:)
      n=size(fv);allocate(tx(size(x,1)))
      do i=2,n
         k=i;tf=fv(i);tx=x(:,i);j=i-1
         do while(j>=1)
            if(fv(j)<=tf)exit
            fv(j+1)=fv(j);x(:,j+1)=x(:,j);j=j-1
         end do
         fv(j+1)=tf;x(:,j+1)=tx
      end do
   end subroutine sort_simplex

   subroutine nelder_mead(f,x0,xbest,fbest,tol,maxiter,step)
      procedure(scalar_objective)::f
      real(dp),intent(in)::x0(:)
      real(dp),intent(out)::xbest(size(x0)),fbest
      real(dp),intent(in),optional::tol,step
      integer,intent(in),optional::maxiter
      integer::n,i,it,nit
      real(dp)::eps,st,fr,fe,fc,spreadf,spreadx
      real(dp),allocatable::x(:,:),fv(:),cent(:),xr(:),xe(:),xc(:)
      real(dp),parameter::alpha=1.0_dp,gamma=2.0_dp,rho=0.5_dp,sigma=0.5_dp
      n=size(x0);allocate(x(n,n+1),fv(n+1),cent(n),xr(n),xe(n),xc(n))
      eps=1e-8_dp;if(present(tol))eps=tol
      nit=2000;if(present(maxiter))nit=maxiter
      st=0.05_dp;if(present(step))st=step
      x(:,1)=x0
      do i=1,n
         x(:,i+1)=x0
         x(i,i+1)=x(i,i+1)+st*max(1.0_dp,abs(x0(i)))
      end do
      do i=1,n+1;fv(i)=f(x(:,i));end do
      do it=1,nit
         call sort_simplex(x,fv)
         spreadf=maxval(abs(fv-fv(1)))
         spreadx=maxval(abs(x-spread(x(:,1),2,n+1)))
         if(spreadf<=eps*max(1.0_dp,abs(fv(1))) .and. spreadx<=sqrt(eps)*max(1.0_dp,maxval(abs(x(:,1)))))exit
         cent=sum(x(:,1:n),dim=2)/real(n,dp)
         xr=cent+alpha*(cent-x(:,n+1));fr=f(xr)
         if(fr<fv(1))then
            xe=cent+gamma*(xr-cent);fe=f(xe)
            if(fe<fr)then;x(:,n+1)=xe;fv(n+1)=fe;else;x(:,n+1)=xr;fv(n+1)=fr;end if
         else if(fr<fv(n))then
            x(:,n+1)=xr;fv(n+1)=fr
         else
            if(fr<fv(n+1))then
               xc=cent+rho*(xr-cent)
            else
               xc=cent-rho*(cent-x(:,n+1))
            end if
            fc=f(xc)
            if(fc<min(fr,fv(n+1)))then
               x(:,n+1)=xc;fv(n+1)=fc
            else
               do i=2,n+1
                  x(:,i)=x(:,1)+sigma*(x(:,i)-x(:,1));fv(i)=f(x(:,i))
               end do
            end if
         end if
      end do
      call sort_simplex(x,fv);xbest=x(:,1);fbest=fv(1)
   end subroutine nelder_mead
end module ks_optimize
