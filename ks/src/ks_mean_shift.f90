! SPDX-License-Identifier: GPL-2.0-only
module ks_mean_shift
   use ks_kinds, only: dp
   use ks_kde, only: kde_model
   use ks_linalg, only: symmetric_eigen
   implicit none
   private
   public :: mean_shift_point, kms, density_ridge_point
contains
   subroutine mean_shift_point(model,start,mode,tol,maxiter,niter)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::start(:)
      real(dp),intent(out)::mode(size(start))
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      integer,intent(out),optional::niter
      real(dp)::eps,y(model%d),yn(model%d),dx(model%d),q,w,sumw
      integer::it,j,nit
      eps=1e-8_dp;if(present(tol))eps=tol;nit=500;if(present(maxiter))nit=maxiter;y=start
      do it=1,nit
         yn=0.0_dp;sumw=0.0_dp
         do j=1,model%n
            dx=y-model%x(j,:);q=dot_product(dx,matmul(model%Hinv,dx))
            if(q<1500.0_dp)then;w=model%w(j)*exp(-0.5_dp*q);yn=yn+w*model%x(j,:);sumw=sumw+w;end if
         end do
         if(sumw<=tiny(1.0_dp))exit
         yn=yn/sumw
         if(sqrt(sum((yn-y)**2))<=eps*max(1.0_dp,sqrt(sum(y*y))))then;y=yn;exit;end if
         y=yn
      end do
      mode=y;if(present(niter))niter=it
   end subroutine mean_shift_point

   subroutine kms(model,starts,modes,labels,tol_merge,tol_iter,maxiter)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::starts(:,:)
      real(dp),intent(out)::modes(size(starts,1),size(starts,2))
      integer,intent(out)::labels(size(starts,1))
      real(dp),intent(in),optional::tol_merge,tol_iter
      integer,intent(in),optional::maxiter
      real(dp)::tm,ti
      integer::i,j,nm
      tm=1e-4_dp;if(present(tol_merge))tm=tol_merge;ti=1e-8_dp;if(present(tol_iter))ti=tol_iter
      do i=1,size(starts,1);call mean_shift_point(model,starts(i,:),modes(i,:),ti,maxiter);end do
      nm=0;labels=0
      do i=1,size(starts,1)
         do j=1,i-1
            if(sqrt(sum((modes(i,:)-modes(j,:))**2))<=tm)then;labels(i)=labels(j);exit;end if
         end do
         if(labels(i)==0)then;nm=nm+1;labels(i)=nm;end if
      end do
   end subroutine kms

   subroutine kde_grad_hess(model,y,grad,hess,dens)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::y(:)
      real(dp),intent(out)::grad(model%d),hess(model%d,model%d),dens
      real(dp)::dx(model%d),u(model%d),q,phi,w
      integer::j,a,b
      grad=0.0_dp;hess=0.0_dp;dens=0.0_dp
      do j=1,model%n
         dx=y-model%x(j,:);u=matmul(model%Hinv,dx);q=dot_product(dx,u)
         if(q>=1500.0_dp)cycle
         phi=exp(model%log_norm-0.5_dp*q);w=model%w(j)/real(model%n,dp)
         dens=dens+w*phi;grad=grad-w*phi*u
         do b=1,model%d;do a=1,model%d
            hess(a,b)=hess(a,b)+w*phi*(u(a)*u(b)-model%Hinv(a,b))
         end do;end do
      end do
   end subroutine kde_grad_hess

   subroutine density_ridge_point(model,start,dim_ridge,point,tol,maxiter)
      type(kde_model),intent(in)::model
      real(dp),intent(in)::start(:)
      integer,intent(in)::dim_ridge
      real(dp),intent(out)::point(size(start))
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      real(dp)::y(model%d),yn(model%d),grad(model%d),hess(model%d,model%d),eig(model%d),vec(model%d,model%d)
      real(dp)::proj(model%d),ms(model%d),dx(model%d),q,w,sumw,dens,eps
      integer::it,nit,k,j,info
      eps=1e-7_dp;if(present(tol))eps=tol;nit=400;if(present(maxiter))nit=maxiter;y=start
      if(dim_ridge<0.or.dim_ridge>=model%d) error stop 'density_ridge_point: invalid ridge dimension'
      do it=1,nit
         call kde_grad_hess(model,y,grad,hess,dens);call symmetric_eigen(hess,eig,vec,info);if(info/=0)exit
         ms=0.0_dp;sumw=0.0_dp
         do j=1,model%n
            dx=y-model%x(j,:);q=dot_product(dx,matmul(model%Hinv,dx))
            if(q<1500.0_dp)then
               w=model%w(j)*exp(-0.5_dp*q);ms=ms+w*model%x(j,:);sumw=sumw+w
            end if
         end do
         if(sumw<=tiny(1.0_dp))exit
         ms=ms/sumw-y
         proj=0.0_dp
         do k=1,model%d-dim_ridge;proj=proj+dot_product(vec(:,k),ms)*vec(:,k);end do
         yn=y+proj
         if(sqrt(sum((yn-y)**2))<eps)then;y=yn;exit;end if
         y=yn
      end do
      point=y
   end subroutine density_ridge_point
end module ks_mean_shift
