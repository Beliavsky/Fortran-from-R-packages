module compoissonreg_regression
   use compoissonreg_kinds, only : dp
   use compoissonreg_types, only : cmp_control_t, cmp_init_t, cmp_fixed_t, cmp_offset_t
   use compoissonreg_types, only : cmp_fit_t, zicmp_fit_t, default_init, default_fixed, default_offset
   use compoissonreg_distributions, only : loglik_cmp, loglik_zicmp, cmp_stats
   use compoissonreg_numerics, only : logistic, invert_matrix
   implicit none
   private
   public :: fit_cmp_raw, fit_zicmp_raw, fitted_cmp, fitted_zicmp

   type :: reg_context_t
      integer, allocatable :: y(:)
      real(dp), allocatable :: xmat(:,:), smat(:,:), wmat(:,:)
      real(dp), allocatable :: offx(:), offs(:), offw(:)
      real(dp), allocatable :: base(:)
      logical, allocatable :: fixed(:)
      integer :: d1=0,d2=0,d3=0
      type(cmp_control_t) :: control
   end type reg_context_t

contains

   subroutine fitted_cmp(xmat,smat,beta,gamma,offx,offs,lambda,nu)
      real(dp),intent(in)::xmat(:,:),smat(:,:),beta(:),gamma(:),offx(:),offs(:)
      real(dp),intent(out)::lambda(:),nu(:)
      lambda=exp(matmul(xmat,beta)+offx)
      nu=exp(matmul(smat,gamma)+offs)
   end subroutine fitted_cmp

   subroutine fitted_zicmp(xmat,smat,wmat,beta,gamma,zeta,offx,offs,offw,lambda,nu,p)
      real(dp),intent(in)::xmat(:,:),smat(:,:),wmat(:,:),beta(:),gamma(:),zeta(:)
      real(dp),intent(in)::offx(:),offs(:),offw(:)
      real(dp),intent(out)::lambda(:),nu(:),p(:)
      integer::i
      lambda=exp(matmul(xmat,beta)+offx)
      nu=exp(matmul(smat,gamma)+offs)
      p=matmul(wmat,zeta)+offw
      do i=1,size(p);p(i)=logistic(p(i));end do
   end subroutine fitted_zicmp

   function reg_loglik(ctx, par) result(ll)
      type(reg_context_t),intent(in)::ctx
      real(dp),intent(in)::par(:)
      real(dp)::ll
      real(dp),allocatable::theta(:),beta(:),gamma(:),zeta(:),lambda(:),nu(:),p(:)
      integer::i,k,n
      n=size(ctx%y)
      allocate(theta(size(ctx%base)));theta=ctx%base;k=0
      do i=1,size(theta)
         if(.not.ctx%fixed(i))then;k=k+1;theta(i)=par(k);end if
      end do
      allocate(beta(ctx%d1),gamma(ctx%d2),zeta(ctx%d3),lambda(n),nu(n),p(n))
      if(ctx%d1>0)beta=theta(1:ctx%d1)
      if(ctx%d2>0)gamma=theta(ctx%d1+1:ctx%d1+ctx%d2)
      if(ctx%d3>0)zeta=theta(ctx%d1+ctx%d2+1:ctx%d1+ctx%d2+ctx%d3)
      if(ctx%d3==0)then
         call fitted_cmp(ctx%xmat,ctx%smat,beta,gamma,ctx%offx,ctx%offs,lambda,nu)
         ll=loglik_cmp(ctx%y,lambda,nu,ctx%control)
      else
         call fitted_zicmp(ctx%xmat,ctx%smat,ctx%wmat,beta,gamma,zeta,ctx%offx,ctx%offs,ctx%offw,lambda,nu,p)
         ll=loglik_zicmp(ctx%y,lambda,nu,p,ctx%control)
      end if
   end function reg_loglik

   subroutine reg_gradient(ctx,x,g)
      type(reg_context_t),intent(in)::ctx
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::g(:)
      real(dp),allocatable::theta(:),beta(:),gamma(:),zeta(:),lambda(:),nu(:),p(:),fullg(:)
      real(dp)::logz,ey,vy,elogf,lfy,score1,score2,score3,f0,den,tau
      integer::i,j,k,n
      n=size(ctx%y)
      allocate(theta(size(ctx%base)),fullg(size(ctx%base)))
      theta=ctx%base;k=0
      do i=1,size(theta)
         if(.not.ctx%fixed(i))then;k=k+1;theta(i)=x(k);end if
      end do
      allocate(beta(ctx%d1),gamma(ctx%d2),zeta(ctx%d3),lambda(n),nu(n),p(n))
      if(ctx%d1>0)beta=theta(1:ctx%d1)
      if(ctx%d2>0)gamma=theta(ctx%d1+1:ctx%d1+ctx%d2)
      if(ctx%d3>0)zeta=theta(ctx%d1+ctx%d2+1:)
      if(ctx%d3==0)then
         call fitted_cmp(ctx%xmat,ctx%smat,beta,gamma,ctx%offx,ctx%offs,lambda,nu)
      else
         call fitted_zicmp(ctx%xmat,ctx%smat,ctx%wmat,beta,gamma,zeta,ctx%offx,ctx%offs,ctx%offw,lambda,nu,p)
      end if
      fullg=0.0_dp
      do i=1,n
         call cmp_stats(lambda(i),nu(i),ctx%control,logz,ey,vy,elogf)
         lfy=log_gamma(real(ctx%y(i)+1,dp))
         if(ctx%d3==0)then
            score1=real(ctx%y(i),dp)-ey
            score2=nu(i)*(elogf-lfy)
            score3=0.0_dp
         else if(ctx%y(i)==0)then
            f0=exp(-logz)
            den=p(i)+(1.0_dp-p(i))*f0
            tau=(1.0_dp-p(i))*f0/max(den,tiny(1.0_dp))
            score1=-tau*ey
            score2=tau*nu(i)*elogf
            score3=p(i)*(1.0_dp-p(i))*(1.0_dp-f0)/max(den,tiny(1.0_dp))
         else
            score1=real(ctx%y(i),dp)-ey
            score2=nu(i)*(elogf-lfy)
            score3=-p(i)
         end if
         do j=1,ctx%d1
            fullg(j)=fullg(j)+ctx%xmat(i,j)*score1
         end do
         do j=1,ctx%d2
            fullg(ctx%d1+j)=fullg(ctx%d1+j)+ctx%smat(i,j)*score2
         end do
         do j=1,ctx%d3
            fullg(ctx%d1+ctx%d2+j)=fullg(ctx%d1+ctx%d2+j)+ctx%wmat(i,j)*score3
         end do
      end do
      k=0
      do i=1,size(fullg)
         if(.not.ctx%fixed(i))then;k=k+1;g(k)=fullg(i);end if
      end do
   end subroutine reg_gradient

   subroutine reg_hessian(ctx,x,hess)
      type(reg_context_t),intent(in)::ctx
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::hess(:,:)
      real(dp),allocatable::xpp(:),xpm(:),xmp(:),xmm(:),xp(:),xm(:)
      real(dp)::hi,hj,f0
      integer::i,j,n
      n=size(x);allocate(xpp(n),xpm(n),xmp(n),xmm(n),xp(n),xm(n));f0=reg_loglik(ctx,x)
      do i=1,n
         hi=sqrt(ctx%control%fd_step)*max(1.0_dp,abs(x(i)))
         xp=x;xm=x;xp(i)=xp(i)+hi;xm(i)=xm(i)-hi
         hess(i,i)=(reg_loglik(ctx,xp)-2.0_dp*f0+reg_loglik(ctx,xm))/(hi*hi)
         do j=i+1,n
            hj=sqrt(ctx%control%fd_step)*max(1.0_dp,abs(x(j)))
            xpp=x;xpm=x;xmp=x;xmm=x
            xpp(i)=xpp(i)+hi;xpp(j)=xpp(j)+hj
            xpm(i)=xpm(i)+hi;xpm(j)=xpm(j)-hj
            xmp(i)=xmp(i)-hi;xmp(j)=xmp(j)+hj
            xmm(i)=xmm(i)-hi;xmm(j)=xmm(j)-hj
            hess(i,j)=(reg_loglik(ctx,xpp)-reg_loglik(ctx,xpm)-reg_loglik(ctx,xmp)+reg_loglik(ctx,xmm)) &
               /(4.0_dp*hi*hj)
            hess(j,i)=hess(i,j)
         end do
      end do
   end subroutine reg_hessian

   subroutine bfgs_maximize(ctx,x,converged,iters)
      type(reg_context_t),intent(in)::ctx
      real(dp),intent(inout)::x(:)
      logical,intent(out)::converged
      integer,intent(out)::iters
      real(dp),allocatable::b(:,:),g(:),gnew(:),d(:),xnew(:),s(:),yv(:),iunit(:,:)
      real(dp)::f,fnew,alpha,gd,ys,rho,gnorm
      integer::n,i,iter
      n=size(x);converged=.false.;iters=0
      if(n==0)then;converged=.true.;return;end if
      allocate(b(n,n),g(n),gnew(n),d(n),xnew(n),s(n),yv(n),iunit(n,n))
      b=0.0_dp;iunit=0.0_dp
      do i=1,n;b(i,i)=1.0_dp;iunit(i,i)=1.0_dp;end do
      call reg_gradient(ctx,x,g);f=reg_loglik(ctx,x)
      do iter=1,ctx%control%max_iter
         iters=iter;gnorm=maxval(abs(g))
         if(gnorm<ctx%control%optim_tol)then;converged=.true.;exit;end if
         d=matmul(b,g)
         gd=dot_product(g,d)
         if(gd<=0.0_dp .or. .not.(gd<huge(1.0_dp)))then;d=g;b=iunit;gd=dot_product(g,d);end if
         if(maxval(abs(d))>2.0_dp) d=d*(2.0_dp/maxval(abs(d)))
         gd=dot_product(g,d)
         alpha=1.0_dp
         do
            xnew=x+alpha*d
            fnew=reg_loglik(ctx,xnew)
            if(fnew>=f+1.0e-4_dp*alpha*gd)exit
            alpha=0.5_dp*alpha
            if(alpha<1.0e-10_dp)exit
         end do
         if(alpha<1.0e-10_dp)then
            ! Reset the inverse Hessian and try a small steepest-ascent step.
            b=iunit;alpha=1.0e-3_dp/max(1.0_dp,sqrt(dot_product(g,g)))
            xnew=x+alpha*g;fnew=reg_loglik(ctx,xnew)
            if(fnew<=f)exit
         end if
         call reg_gradient(ctx,xnew,gnew)
         s=xnew-x;yv=gnew-g;ys=dot_product(yv,s)
         if(ys>1.0e-12_dp)then
            ! Maximization: inverse of negative Hessian uses y = g_old - g_new.
            yv=-yv;ys=dot_product(yv,s)
            if(ys>1.0e-12_dp)then
               rho=1.0_dp/ys
               b=matmul(iunit-rho*outer(s,yv),matmul(b,iunit-rho*outer(yv,s)))+rho*outer(s,s)
            else
               b=iunit
            end if
         else
            b=iunit
         end if
         x=xnew;g=gnew;f=fnew
         if(maxval(abs(s))<ctx%control%optim_tol*(1.0_dp+maxval(abs(x))))then
            converged=.true.;exit
         end if
      end do
   contains
      pure function outer(a,bv) result(c)
         real(dp),intent(in)::a(:),bv(:)
         real(dp)::c(size(a),size(bv))
         integer::ii,jj
         do jj=1,size(bv);do ii=1,size(a);c(ii,jj)=a(ii)*bv(jj);end do;end do
      end function outer
   end subroutine bfgs_maximize

   subroutine prepare_context(y,xmat,smat,wmat,init,fixed,offset,control,ctx,par)
      integer,intent(in)::y(:)
      real(dp),intent(in)::xmat(:,:),smat(:,:),wmat(:,:)
      type(cmp_init_t),intent(in)::init
      type(cmp_fixed_t),intent(in)::fixed
      type(cmp_offset_t),intent(in)::offset
      type(cmp_control_t),intent(in)::control
      type(reg_context_t),intent(out)::ctx
      real(dp),allocatable,intent(out)::par(:)
      integer::q,i,k
      ctx%d1=size(xmat,2);ctx%d2=size(smat,2);ctx%d3=size(wmat,2);ctx%control=control
      ctx%y=y;ctx%xmat=xmat;ctx%smat=smat;ctx%wmat=wmat
      ctx%offx=offset%x;ctx%offs=offset%s;ctx%offw=offset%w
      allocate(ctx%base(ctx%d1+ctx%d2+ctx%d3),ctx%fixed(ctx%d1+ctx%d2+ctx%d3))
      if(ctx%d1>0)then;ctx%base(1:ctx%d1)=init%beta;ctx%fixed(1:ctx%d1)=fixed%beta;end if
      if(ctx%d2>0)then
         ctx%base(ctx%d1+1:ctx%d1+ctx%d2)=init%gamma
         ctx%fixed(ctx%d1+1:ctx%d1+ctx%d2)=fixed%gamma
      end if
      if(ctx%d3>0)then
         ctx%base(ctx%d1+ctx%d2+1:)=init%zeta
         ctx%fixed(ctx%d1+ctx%d2+1:)=fixed%zeta
      end if
      q=count(.not.ctx%fixed);allocate(par(q));k=0
      do i=1,size(ctx%base);if(.not.ctx%fixed(i))then;k=k+1;par(k)=ctx%base(i);end if;end do
   end subroutine prepare_context

   subroutine unpack_solution(ctx,par,beta,gamma,zeta)
      type(reg_context_t),intent(in)::ctx
      real(dp),intent(in)::par(:)
      real(dp),intent(out)::beta(:),gamma(:),zeta(:)
      real(dp),allocatable::theta(:)
      integer::i,k
      theta=ctx%base;k=0
      do i=1,size(theta);if(.not.ctx%fixed(i))then;k=k+1;theta(i)=par(k);end if;end do
      if(ctx%d1>0)beta=theta(1:ctx%d1)
      if(ctx%d2>0)gamma=theta(ctx%d1+1:ctx%d1+ctx%d2)
      if(ctx%d3>0)zeta=theta(ctx%d1+ctx%d2+1:)
   end subroutine unpack_solution

   subroutine fit_cmp_raw(y,xmat,smat,fit,init,fixed,offset,control)
      integer,intent(in)::y(:)
      real(dp),intent(in)::xmat(:,:),smat(:,:)
      type(cmp_fit_t),intent(out)::fit
      type(cmp_init_t),intent(in),optional::init
      type(cmp_fixed_t),intent(in),optional::fixed
      type(cmp_offset_t),intent(in),optional::offset
      type(cmp_control_t),intent(in),optional::control
      type(cmp_init_t)::ini;type(cmp_fixed_t)::fix;type(cmp_offset_t)::off
      type(cmp_control_t)::ctrl;type(reg_context_t)::ctx
      real(dp),allocatable::wmat(:,:),par(:),zeta(:),neg_h(:,:)
      logical::ok
      integer::d1,d2,n,q
      n=size(y);d1=size(xmat,2);d2=size(smat,2)
      if(size(xmat,1)/=n.or.size(smat,1)/=n)error stop 'fit_cmp_raw: row mismatch'
      ini=default_init(d1,d2);if(present(init))ini=init
      fix=default_fixed(d1,d2);if(present(fixed))fix=fixed
      off=default_offset(n);if(present(offset))off=offset
      ctrl=cmp_control_t();if(present(control))ctrl=control
      allocate(wmat(n,0));call prepare_context(y,xmat,smat,wmat,ini,fix,off,ctrl,ctx,par)
      call bfgs_maximize(ctx,par,fit%converged,fit%iterations)
      allocate(fit%beta(d1),fit%gamma(d2),zeta(0))
      call unpack_solution(ctx,par,fit%beta,fit%gamma,zeta)
      fit%loglik=reg_loglik(ctx,par);q=size(par)
      allocate(fit%hessian(q,q),fit%covariance(q,q))
      if(q>0)then
         call reg_hessian(ctx,par,fit%hessian);neg_h=-fit%hessian
         call invert_matrix(neg_h,fit%covariance,ok)
         if(.not.ok)fit%covariance=huge(1.0_dp)
      end if
      fit%y=y;fit%xmat=xmat;fit%smat=smat;fit%offset=off;fit%fixed=fix;fit%control=ctrl
   end subroutine fit_cmp_raw

   subroutine fit_zicmp_raw(y,xmat,smat,wmat,fit,init,fixed,offset,control)
      integer,intent(in)::y(:)
      real(dp),intent(in)::xmat(:,:),smat(:,:),wmat(:,:)
      type(zicmp_fit_t),intent(out)::fit
      type(cmp_init_t),intent(in),optional::init
      type(cmp_fixed_t),intent(in),optional::fixed
      type(cmp_offset_t),intent(in),optional::offset
      type(cmp_control_t),intent(in),optional::control
      type(cmp_init_t)::ini;type(cmp_fixed_t)::fix;type(cmp_offset_t)::off
      type(cmp_control_t)::ctrl;type(reg_context_t)::ctx
      real(dp),allocatable::par(:),neg_h(:,:)
      logical::ok
      integer::d1,d2,d3,n,q
      n=size(y);d1=size(xmat,2);d2=size(smat,2);d3=size(wmat,2)
      if(size(xmat,1)/=n.or.size(smat,1)/=n.or.size(wmat,1)/=n) &
         error stop 'fit_zicmp_raw: row mismatch'
      ini=default_init(d1,d2,d3);if(present(init))ini=init
      fix=default_fixed(d1,d2,d3);if(present(fixed))fix=fixed
      off=default_offset(n);if(present(offset))off=offset
      ctrl=cmp_control_t();if(present(control))ctrl=control
      call prepare_context(y,xmat,smat,wmat,ini,fix,off,ctrl,ctx,par)
      call bfgs_maximize(ctx,par,fit%converged,fit%iterations)
      allocate(fit%beta(d1),fit%gamma(d2),fit%zeta(d3))
      call unpack_solution(ctx,par,fit%beta,fit%gamma,fit%zeta)
      fit%loglik=reg_loglik(ctx,par);q=size(par)
      allocate(fit%hessian(q,q),fit%covariance(q,q))
      if(q>0)then
         call reg_hessian(ctx,par,fit%hessian);neg_h=-fit%hessian
         call invert_matrix(neg_h,fit%covariance,ok)
         if(.not.ok)fit%covariance=huge(1.0_dp)
      end if
      fit%y=y;fit%xmat=xmat;fit%smat=smat;fit%wmat=wmat;fit%offset=off;fit%fixed=fix;fit%control=ctrl
   end subroutine fit_zicmp_raw

end module compoissonreg_regression
