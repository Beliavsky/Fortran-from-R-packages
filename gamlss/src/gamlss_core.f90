! Matrix-first GAMLSS fitting engine.
! Numerical implementation of the Rigby-Stasinopoulos (RS) and
! Cole-Green (CG) iteration structures used by upstream gamlss.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar, map_parameters, inverse_link, &
      family_logpdf, default_parameters
   use gamlss_linalg, only : solve_linear, invert_matrix, penalized_weighted_least_squares, matrix_rank
   use gamlss_types
   implicit none
   private
   public :: fit_gamlss_model, predict_gamlss_parameters

   type :: fit_block_t
      real(dp), allocatable :: x(:,:), penalty(:,:), offset(:)
      real(dp), allocatable :: beta(:), eta(:), covariance(:,:)
      real(dp), allocatable :: last_weights(:)
      real(dp) :: lambda = 0.0_dp
      real(dp) :: step = 1.0_dp
      logical :: fixed = .false.
      logical :: estimate_lambda = .false.
   end type fit_block_t

contains

   subroutine fit_gamlss_model(y,x_mu,family,result,method,x_sigma,x_nu,x_tau, &
      weights,offset_mu,offset_sigma,offset_nu,offset_tau, &
      penalty_mu,penalty_sigma,penalty_nu,penalty_tau, &
      lambda_mu,lambda_sigma,lambda_nu,lambda_tau,start, &
      fix_mu,fix_sigma,fix_nu,fix_tau,control)
      real(dp),intent(in)::y(:),x_mu(:,:)
      integer,intent(in)::family
      type(gamlss_result_t),intent(out)::result
      integer,intent(in),optional::method
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      real(dp),intent(in),optional::penalty_mu(:,:),penalty_sigma(:,:),penalty_nu(:,:),penalty_tau(:,:)
      real(dp),intent(in),optional::lambda_mu,lambda_sigma,lambda_nu,lambda_tau,start(:)
      logical,intent(in),optional::fix_mu,fix_sigma,fix_nu,fix_tau
      type(gamlss_control_t),intent(in),optional::control

      type(fit_block_t)::b(4)
      type(gamlss_control_t)::ctrl
      real(dp),allocatable::w(:)
      integer::n,np,meth,status,ncoef,meth_index

      result%family=family
      result%status=0
      n=size(y); np=family_npar(family)
      if(n<=0 .or. size(x_mu,1)/=n .or. np<1 .or. np>4)then
         result%status=10; return
      end if
      allocate(w(n)); w=1.0_dp
      if(present(weights))then
         if(size(weights)/=n .or. any(weights<0.0_dp))then
            result%status=11; return
         end if
         w=weights
      end if
      ctrl=gamlss_control_t(); if(present(control))ctrl=control
      meth=GAMLSS_METHOD_RS; if(present(method))meth=method
      result%method=meth

      call setup_block(b(1),n,x_mu,offset_mu,penalty_mu,lambda_mu,ctrl%mu_step,fix_mu, &
         ctrl%estimate_lambda_mu)
      if(np>=2)call setup_block(b(2),n,x_sigma,offset_sigma,penalty_sigma,lambda_sigma, &
         ctrl%sigma_step,fix_sigma,ctrl%estimate_lambda_sigma)
      if(np>=3)call setup_block(b(3),n,x_nu,offset_nu,penalty_nu,lambda_nu, &
         ctrl%nu_step,fix_nu,ctrl%estimate_lambda_nu)
      if(np>=4)call setup_block(b(4),n,x_tau,offset_tau,penalty_tau,lambda_tau, &
         ctrl%tau_step,fix_tau,ctrl%estimate_lambda_tau)
      if(any([(size(b(meth_index)%x,1)/=n,meth_index=1,np)]))then
         result%status=12; return
      end if
      ncoef=0
      block
         integer::j
         do j=1,np; ncoef=ncoef+size(b(j)%x,2); end do
      end block
      call initialize_blocks(b,np,y,family,start,status)
      if(status/=0)then; result%status=status; return; end if

      select case(meth)
      case(GAMLSS_METHOD_RS)
         call run_rs(b,np,y,w,family,ctrl,ctrl%n_cyc,result%iterations,result%converged,status)
      case(GAMLSS_METHOD_CG)
         ! Upstream CG is entered after parameter initialisation.  One blockwise
         ! Fisher pass provides the corresponding stable GLIM starting point.
         block
            integer :: init_it,cg_it
            logical :: init_cv,cg_cv
            cg_it=0;cg_cv=.false.
            call run_rs(b,np,y,w,family,ctrl,1,init_it,init_cv,status)
            if(status==0)call run_cg(b,np,y,w,family,ctrl,ctrl%n_cyc,cg_it,cg_cv,status)
            result%iterations=init_it+merge(cg_it,0,status==0)
            result%converged=(status==0 .and. cg_cv)
         end block
      case(GAMLSS_METHOD_MIXED)
         call run_rs(b,np,y,w,family,ctrl,max(1,ctrl%mixed_rs_cycles),result%iterations, &
            result%converged,status)
         if(status==0)then
            block
               integer::it2
               logical::cv2
               call run_cg(b,np,y,w,family,ctrl,max(1,ctrl%mixed_cg_cycles),it2,cv2,status)
               result%iterations=result%iterations+it2
               result%converged=cv2
            end block
         end if
      case default
         status=13
      end select
      result%status=status
      if(status/=0)return
      call finalize_result(b,np,y,w,family,result)
   end subroutine fit_gamlss_model

   subroutine setup_block(block,n,x,offset,penalty,lambda,step,fixed,estimate_lambda)
      type(fit_block_t),intent(out)::block
      integer,intent(in)::n
      real(dp),intent(in),optional::x(:,:),offset(:),penalty(:,:),lambda,step
      logical,intent(in),optional::fixed,estimate_lambda
      integer::p,i
      if(present(x))then
         if(size(x,1)/=n)then
            allocate(block%x(0,0)); return
         end if
         block%x=x
      else
         allocate(block%x(n,1)); block%x=1.0_dp
      end if
      p=size(block%x,2)
      allocate(block%offset(n)); block%offset=0.0_dp
      if(present(offset))then
         if(size(offset)==n)block%offset=offset
      end if
      allocate(block%penalty(p,p)); block%penalty=0.0_dp
      if(present(penalty))then
         if(size(penalty,1)==p .and. size(penalty,2)==p)block%penalty=penalty
      end if
      block%lambda=0.0_dp; if(present(lambda))block%lambda=max(0.0_dp,lambda)
      block%step=1.0_dp; if(present(step))block%step=step
      block%fixed=.false.; if(present(fixed))block%fixed=fixed
      block%estimate_lambda=.false.; if(present(estimate_lambda))block%estimate_lambda=estimate_lambda
      if(block%estimate_lambda .and. block%lambda<=0.0_dp)block%lambda=1.0_dp
      allocate(block%beta(p),block%eta(n),block%last_weights(n))
      block%beta=0.0_dp; block%eta=block%offset; block%last_weights=1.0_dp
      allocate(block%covariance(p,p)); block%covariance=0.0_dp
      do i=1,p
         if(.not.all(ieee_is_finite(block%x(:,i))))then
            block%x(:,i)=0.0_dp
         end if
      end do
   end subroutine setup_block

   subroutine initialize_blocks(b,np,y,family,start,status)
      type(fit_block_t),intent(inout)::b(4)
      integer,intent(in)::np,family
      real(dp),intent(in)::y(:)
      real(dp),intent(in),optional::start(:)
      integer,intent(out)::status
      real(dp)::m,s,a0,b0,c0,d0,pv(4),e0
      integer::j,p,pos,total
      status=0; total=0
      do j=1,np; total=total+size(b(j)%beta); end do
      if(present(start))then
         if(size(start)/=total)then; status=20; return; end if
         pos=0
         do j=1,np
            p=size(b(j)%beta); b(j)%beta=start(pos+1:pos+p); pos=pos+p
            b(j)%eta=matmul(b(j)%x,b(j)%beta)+b(j)%offset
         end do
         return
      end if
      m=sum(y)/real(size(y),dp)
      s=sqrt(max(sum((y-m)**2)/real(max(1,size(y)-1),dp),1.0e-8_dp))
      call default_parameters(family,m,s,a0,b0,c0,d0)
      pv=[a0,b0,c0,d0]
      do j=1,np
         b(j)%beta=0.0_dp
         e0=inverse_link(family,j,pv(j))
         if(size(b(j)%beta)>=1)then
            if(maxval(abs(b(j)%x(:,1)-1.0_dp))<1.0e-10_dp) &
               b(j)%beta(1)=e0-sum(b(j)%offset)/real(size(y),dp)
         end if
         b(j)%eta=matmul(b(j)%x,b(j)%beta)+b(j)%offset
      end do
   end subroutine initialize_blocks

   subroutine run_rs(b,np,y,w,family,ctrl,ncycles,iters,converged,status)
      type(fit_block_t),intent(inout)::b(4)
      integer,intent(in)::np,family,ncycles
      real(dp),intent(in)::y(:),w(:)
      type(gamlss_control_t),intent(in)::ctrl
      integer,intent(out)::iters,status
      logical,intent(out)::converged
      real(dp)::olddev,newdev
      integer::it,j
      status=0; converged=.false.; iters=0
      olddev=global_deviance(b,np,y,w,family)
      if(.not.ieee_is_finite(olddev))then; status=30; return; end if
      do it=1,ncycles
         do j=1,np
            if(.not.b(j)%fixed)then
               call update_block_rs(j,b,np,y,w,family,ctrl,status)
               if(status/=0)return
            end if
         end do
         newdev=global_deviance(b,np,y,w,family)
         iters=it
         if(abs(olddev-newdev)<ctrl%c_crit)then
            converged=.true.; exit
         end if
         olddev=newdev
      end do
   end subroutine run_rs

   subroutine update_block_rs(j,b,np,y,w,family,ctrl,status)
      integer,intent(in)::j,np,family
      type(fit_block_t),intent(inout)::b(4)
      real(dp),intent(in)::y(:),w(:)
      type(gamlss_control_t),intent(in)::ctrl
      integer,intent(out)::status
      real(dp),allocatable::z(:),ww(:),beta_fit(:),cov(:,:),oldbeta(:),oldeta(:),pen(:,:)
      real(dp)::olddev,newdev,frac,change,newlambda
      integer::inner,k,istat
      status=0
      allocate(z(size(y)),ww(size(y)),oldbeta(size(b(j)%beta)),oldeta(size(y)))
      allocate(pen(size(b(j)%beta),size(b(j)%beta)))
      pen=b(j)%lambda*b(j)%penalty
      olddev=global_deviance(b,np,y,w,family)
      do inner=1,ctrl%inner_cyc
         call working_values(j,b,np,y,w,family,z,ww,status)
         if(status/=0)return
         b(j)%last_weights=ww
         oldbeta=b(j)%beta; oldeta=b(j)%eta
         call penalized_weighted_least_squares(b(j)%x,z-b(j)%offset,ww,pen,beta_fit,cov,istat)
         if(istat/=0)then; status=31; return; end if
         frac=max(0.0_dp,min(1.0_dp,b(j)%step))
         b(j)%beta=oldbeta+frac*(beta_fit-oldbeta)
         b(j)%eta=matmul(b(j)%x,b(j)%beta)+b(j)%offset
         newdev=global_deviance(b,np,y,w,family)
         if(ctrl%autostep .and. newdev>olddev)then
            do k=1,8
               b(j)%beta=0.5_dp*(b(j)%beta+oldbeta)
               b(j)%eta=matmul(b(j)%x,b(j)%beta)+b(j)%offset
               newdev=global_deviance(b,np,y,w,family)
               if(newdev<=olddev)exit
            end do
         end if
         change=abs(newdev-olddev)
         b(j)%covariance=cov
         if(b(j)%estimate_lambda .and. maxval(abs(b(j)%penalty))>0.0_dp)then
            call estimate_block_lambda(b(j),z,ww,ctrl,newlambda,istat)
            if(istat==0)then
               if(abs(log(max(newlambda,ctrl%lambda_min))-log(max(b(j)%lambda,ctrl%lambda_min))) &
                  > ctrl%lambda_crit)then
                  ! Geometric damping is stable for variance-ratio updates spanning orders of magnitude.
                  b(j)%lambda=sqrt(max(ctrl%lambda_min,b(j)%lambda)*newlambda)
               else
                  b(j)%lambda=newlambda
               end if
               b(j)%lambda=min(ctrl%lambda_max,max(ctrl%lambda_min,b(j)%lambda))
               pen=b(j)%lambda*b(j)%penalty
            end if
         end if
         olddev=newdev
         if(change<ctrl%inner_crit)exit
      end do
   end subroutine update_block_rs


   subroutine estimate_block_lambda(block,z,ww,ctrl,lambda,status)
      type(fit_block_t),intent(in)::block
      real(dp),intent(in)::z(:),ww(:)
      type(gamlss_control_t),intent(in)::ctrl
      real(dp),intent(out)::lambda
      integer,intent(out)::status
      real(dp),allocatable::xtwx(:,:),a(:,:),ainv(:,:),hatcoef(:,:),res(:)
      real(dp)::edf,nullity,edf_pen,sig2,tau2,qform
      integer::i,r,c,p,n,rankp,istat
      p=size(block%beta); n=size(z); status=0
      allocate(xtwx(p,p)); xtwx=0.0_dp
      do i=1,n
         do r=1,p
            do c=1,p
               xtwx(r,c)=xtwx(r,c)+ww(i)*block%x(i,r)*block%x(i,c)
            end do
         end do
      end do
      a=xtwx+max(ctrl%lambda_min,block%lambda)*block%penalty
      call invert_matrix(a,ainv,istat)
      if(istat/=0)then; status=1; lambda=block%lambda; return; end if
      hatcoef=matmul(ainv,xtwx)
      edf=0.0_dp
      do i=1,p; edf=edf+hatcoef(i,i); end do
      rankp=matrix_rank(block%penalty)
      nullity=real(max(0,p-rankp),dp)
      edf_pen=max(1.0e-6_dp,edf-nullity)
      allocate(res(n)); res=z-block%offset-matmul(block%x,block%beta)
      sig2=sum(ww*res*res)/max(1.0_dp,sum(merge(1.0_dp,0.0_dp,ww>0.0_dp))-edf)
      qform=dot_product(block%beta,matmul(block%penalty,block%beta))
      tau2=max(1.0e-12_dp,qform/edf_pen)
      lambda=min(ctrl%lambda_max,max(ctrl%lambda_min,sig2/tau2))
   end subroutine estimate_block_lambda

   subroutine working_values(j,b,np,y,w,family,z,ww,status)
      integer,intent(in)::j,np,family
      type(fit_block_t),intent(in)::b(4)
      real(dp),intent(in)::y(:),w(:)
      real(dp),intent(out)::z(:),ww(:)
      integer,intent(out)::status
      real(dp)::score,hess
      integer::i
      status=0
      do i=1,size(y)
         call eta_derivative_1(i,j,b,np,y(i),family,score,hess,status)
         if(status/=0)return
         hess=min(hess,-1.0e-10_dp)
         ww(i)=w(i)*min(1.0e10_dp,max(1.0e-10_dp,-hess))
         z(i)=b(j)%eta(i)+score/max(1.0e-10_dp,-hess)
         if(.not.ieee_is_finite(z(i)))z(i)=b(j)%eta(i)
      end do
   end subroutine working_values

   subroutine eta_derivative_1(i,j,b,np,yi,family,score,hess,status)
      integer,intent(in)::i,j,np,family
      type(fit_block_t),intent(in)::b(4)
      real(dp),intent(in)::yi
      real(dp),intent(out)::score,hess
      integer,intent(out)::status
      real(dp)::e(4),h,l0,lp,lm
      integer::k
      e=0.0_dp
      do k=1,np;e(k)=b(k)%eta(i);end do
      h=1.0e-4_dp*(1.0_dp+abs(e(j)))
      l0=loglik_eta(yi,family,e)
      e(j)=e(j)+h; lp=loglik_eta(yi,family,e)
      e(j)=e(j)-2.0_dp*h; lm=loglik_eta(yi,family,e)
      if(.not.(ieee_is_finite(l0).and.ieee_is_finite(lp).and.ieee_is_finite(lm)))then
         score=0.0_dp; hess=-1.0e-8_dp; status=0; return
      end if
      score=(lp-lm)/(2.0_dp*h)
      hess=(lp-2.0_dp*l0+lm)/(h*h)
      status=0
   end subroutine eta_derivative_1

   subroutine eta_score_hessian(i,b,np,yi,family,score,hess)
      integer,intent(in)::i,np,family
      type(fit_block_t),intent(in)::b(4)
      real(dp),intent(in)::yi
      real(dp),intent(out)::score(4),hess(4,4)
      real(dp)::e(4),ep(4),h(4),l0,lp,lm,lpp,lpm,lmp,lmm
      integer::j,k,r
      e=0.0_dp; score=0.0_dp; hess=0.0_dp
      do r=1,np;e(r)=b(r)%eta(i);h(r)=1.0e-4_dp*(1.0_dp+abs(e(r)));end do
      l0=loglik_eta(yi,family,e)
      do j=1,np
         ep=e; ep(j)=e(j)+h(j); lp=loglik_eta(yi,family,ep)
         ep=e; ep(j)=e(j)-h(j); lm=loglik_eta(yi,family,ep)
         if(ieee_is_finite(lp).and.ieee_is_finite(lm).and.ieee_is_finite(l0))then
            score(j)=(lp-lm)/(2.0_dp*h(j))
            hess(j,j)=(lp-2.0_dp*l0+lm)/(h(j)*h(j))
         else
            hess(j,j)=-1.0e-8_dp
         end if
      end do
      do j=1,np-1
         do k=j+1,np
            ep=e; ep(j)=e(j)+h(j); ep(k)=e(k)+h(k); lpp=loglik_eta(yi,family,ep)
            ep=e; ep(j)=e(j)+h(j); ep(k)=e(k)-h(k); lpm=loglik_eta(yi,family,ep)
            ep=e; ep(j)=e(j)-h(j); ep(k)=e(k)+h(k); lmp=loglik_eta(yi,family,ep)
            ep=e; ep(j)=e(j)-h(j); ep(k)=e(k)-h(k); lmm=loglik_eta(yi,family,ep)
            if(ieee_is_finite(lpp).and.ieee_is_finite(lpm).and.ieee_is_finite(lmp).and.ieee_is_finite(lmm))then
               hess(j,k)=(lpp-lpm-lmp+lmm)/(4.0_dp*h(j)*h(k))
               hess(k,j)=hess(j,k)
            end if
         end do
      end do
   end subroutine eta_score_hessian

   subroutine run_cg(b,np,y,w,family,ctrl,ncycles,iters,converged,status)
      type(fit_block_t),intent(inout)::b(4)
      integer,intent(in)::np,family,ncycles
      real(dp),intent(in)::y(:),w(:)
      type(gamlss_control_t),intent(in)::ctrl
      integer,intent(out)::iters,status
      logical,intent(out)::converged
      real(dp),allocatable::grad(:),info(:,:),delta(:),theta(:),oldtheta(:)
      real(dp)::oldobj,newobj,frac
      integer::it,istat,k,total
      total=total_coefficients(b,np)
      allocate(grad(total),info(total,total),theta(total),oldtheta(total))
      call pack_beta(b,np,theta)
      oldobj=penalized_deviance(b,np,y,w,family)
      status=0; converged=.false.; iters=0
      do it=1,ncycles
         call cg_information(b,np,y,w,family,grad,info)
         call solve_linear(info,grad,delta,istat)
         if(istat/=0)then; status=40; return; end if
         oldtheta=theta; frac=1.0_dp
         theta=oldtheta+frac*delta
         call unpack_beta(b,np,theta)
         newobj=penalized_deviance(b,np,y,w,family)
         if(ctrl%autostep .and. newobj>oldobj)then
            do k=1,10
               frac=0.5_dp*frac
               theta=oldtheta+frac*delta
               call unpack_beta(b,np,theta)
               newobj=penalized_deviance(b,np,y,w,family)
               if(newobj<=oldobj)exit
            end do
         end if
         iters=it
         if(abs(oldobj-newobj)<ctrl%c_crit)then; converged=.true.; exit; end if
         oldobj=newobj
      end do
      call cg_covariance(b,np,y,w,family)
   end subroutine run_cg

   subroutine cg_information(b,np,y,w,family,grad,info)
      type(fit_block_t),intent(in)::b(4)
      integer,intent(in)::np,family
      real(dp),intent(in)::y(:),w(:)
      real(dp),intent(out)::grad(:),info(:,:)
      real(dp)::score(4),hess(4,4)
      integer::off(4),p(4),i,j,k,a,c,oj,ok
      call block_offsets(b,np,off,p)
      grad=0.0_dp; info=0.0_dp
      do i=1,size(y)
         call eta_score_hessian(i,b,np,y(i),family,score,hess)
         do j=1,np
            oj=off(j)
            do a=1,p(j)
               grad(oj+a)=grad(oj+a)+w(i)*score(j)*b(j)%x(i,a)
            end do
            do k=1,np
               ok=off(k)
               do a=1,p(j)
                  do c=1,p(k)
                     info(oj+a,ok+c)=info(oj+a,ok+c)-w(i)*hess(j,k)*b(j)%x(i,a)*b(k)%x(i,c)
                  end do
               end do
            end do
         end do
      end do
      do j=1,np
         oj=off(j)
         if(b(j)%lambda>0.0_dp)then
            grad(oj+1:oj+p(j))=grad(oj+1:oj+p(j))- &
               b(j)%lambda*matmul(b(j)%penalty,b(j)%beta)
            info(oj+1:oj+p(j),oj+1:oj+p(j))=info(oj+1:oj+p(j),oj+1:oj+p(j))+ &
               b(j)%lambda*b(j)%penalty
         end if
         if(b(j)%fixed)then
            do a=1,p(j)
               info(oj+a,:)=0.0_dp; info(:,oj+a)=0.0_dp
               info(oj+a,oj+a)=1.0_dp; grad(oj+a)=0.0_dp
            end do
         end if
      end do
      do a=1,size(grad); info(a,a)=info(a,a)+1.0e-9_dp; end do
   end subroutine cg_information

   subroutine cg_covariance(b,np,y,w,family)
      type(fit_block_t),intent(inout)::b(4)
      integer,intent(in)::np,family
      real(dp),intent(in)::y(:),w(:)
      real(dp),allocatable::g(:),inf(:,:),cov(:,:)
      integer::off(4),p(4),j,istat
      allocate(g(total_coefficients(b,np)),inf(total_coefficients(b,np),total_coefficients(b,np)))
      call cg_information(b,np,y,w,family,g,inf)
      call invert_matrix(inf,cov,istat)
      if(istat/=0)return
      call block_offsets(b,np,off,p)
      do j=1,np
         if(allocated(b(j)%covariance))deallocate(b(j)%covariance)
         allocate(b(j)%covariance(p(j),p(j)))
         b(j)%covariance=cov(off(j)+1:off(j)+p(j),off(j)+1:off(j)+p(j))
      end do
   end subroutine cg_covariance

   real(dp) function loglik_eta(y,family,e) result(lp)
      real(dp),intent(in)::y,e(4)
      integer,intent(in)::family
      real(dp)::a,b,c,d
      call map_parameters(family,e(1),e(2),e(3),e(4),a,b,c,d)
      lp=family_logpdf(family,y,a,b,c,d)
   end function loglik_eta

   real(dp) function global_deviance(b,np,y,w,family) result(dev)
      type(fit_block_t),intent(in)::b(4)
      integer,intent(in)::np,family
      real(dp),intent(in)::y(:),w(:)
      real(dp)::e(4),lp
      integer::i,j
      dev=0.0_dp
      do i=1,size(y)
         e=0.0_dp; do j=1,np;e(j)=b(j)%eta(i);end do
         lp=loglik_eta(y(i),family,e)
         if(.not.ieee_is_finite(lp))then; dev=huge(1.0_dp)/100.0_dp; return; end if
         dev=dev-2.0_dp*w(i)*lp
      end do
   end function global_deviance

   real(dp) function penalty_value(b,np) result(v)
      type(fit_block_t),intent(in)::b(4)
      integer,intent(in)::np
      integer::j
      v=0.0_dp
      do j=1,np
         if(b(j)%lambda>0.0_dp)v=v+b(j)%lambda*dot_product(b(j)%beta,matmul(b(j)%penalty,b(j)%beta))
      end do
   end function penalty_value

   real(dp) function penalized_deviance(b,np,y,w,family) result(v)
      type(fit_block_t),intent(in)::b(4)
      integer,intent(in)::np,family
      real(dp),intent(in)::y(:),w(:)
      v=global_deviance(b,np,y,w,family)+penalty_value(b,np)
   end function penalized_deviance

   subroutine finalize_result(b,np,y,w,family,result)
      type(fit_block_t),intent(inout)::b(4)
      integer,intent(in)::np,family
      real(dp),intent(in)::y(:),w(:)
      type(gamlss_result_t),intent(inout)::result
      real(dp)::e(4),a,c,d,s,lp,score,hess
      integer::i,j,istat
      result%global_deviance=global_deviance(b,np,y,w,family)
      result%penalized_deviance=result%global_deviance+penalty_value(b,np)
      call compute_edf(b,np,y,w,family)
      result%df_fit=0.0_dp
      do j=1,np; result%df_fit=result%df_fit+block_edf(b(j)); end do
      result%df_residual=sum(w)-result%df_fit
      result%aic=result%global_deviance+2.0_dp*result%df_fit
      result%sbc=result%global_deviance+log(max(1.0_dp,sum(w)))*result%df_fit
      allocate(result%case_deviance(size(y)),result%residuals(size(y)))
      do i=1,size(y)
         e=0.0_dp; do j=1,np;e(j)=b(j)%eta(i);end do
         lp=loglik_eta(y(i),family,e); result%case_deviance(i)=-2.0_dp*lp
         call eta_derivative_1(i,1,b,np,y(i),family,score,hess,istat)
         result%residuals(i)=score/sqrt(max(1.0e-12_dp,-hess))
      end do
      call copy_parameter(b(1),result%mu)
      if(np>=2)call copy_parameter(b(2),result%sigma)
      if(np>=3)call copy_parameter(b(3),result%nu)
      if(np>=4)call copy_parameter(b(4),result%tau)
      do i=1,size(y)
         e=0.0_dp; do j=1,np;e(j)=b(j)%eta(i);end do
         call map_parameters(family,e(1),e(2),e(3),e(4),a,s,c,d)
         result%mu%fitted(i)=a
         if(np>=2)result%sigma%fitted(i)=s
         if(np>=3)result%nu%fitted(i)=c
         if(np>=4)result%tau%fitted(i)=d
      end do
   end subroutine finalize_result

   subroutine compute_edf(b,np,y,w,family)
      type(fit_block_t),intent(inout)::b(4)
      integer,intent(in)::np,family
      real(dp),intent(in)::y(:),w(:)
      real(dp),allocatable::z(:),ww(:),xtwx(:,:),a(:,:),ainv(:,:)
      integer::j,r,c,i,istat
      do j=1,np
         allocate(z(size(y)),ww(size(y)))
         call working_values(j,b,np,y,w,family,z,ww,istat)
         if(istat/=0)then; deallocate(z,ww); cycle; end if
         allocate(xtwx(size(b(j)%beta),size(b(j)%beta))); xtwx=0.0_dp
         do i=1,size(y)
            do r=1,size(b(j)%beta)
               do c=1,size(b(j)%beta)
                  xtwx(r,c)=xtwx(r,c)+ww(i)*b(j)%x(i,r)*b(j)%x(i,c)
               end do
            end do
         end do
         a=xtwx+b(j)%lambda*b(j)%penalty
         call invert_matrix(a,ainv,istat)
         b(j)%last_weights=ww
         if(istat==0)then
            ! store EDF temporarily in covariance(1,1) is not acceptable; recomputed by block_edf.
            if(.not.allocated(b(j)%covariance))then
               allocate(b(j)%covariance(size(a,1),size(a,2))); b(j)%covariance=ainv
            end if
         end if
         deallocate(z,ww,xtwx,a)
         if(allocated(ainv))deallocate(ainv)
      end do
   end subroutine compute_edf

   real(dp) function block_edf(block) result(edf)
      type(fit_block_t),intent(in)::block
      real(dp),allocatable::xtwx(:,:),a(:,:),ainv(:,:)
      integer::i,r,c,istat,p
      p=size(block%beta); allocate(xtwx(p,p)); xtwx=0.0_dp
      do i=1,size(block%x,1)
         do r=1,p; do c=1,p
            xtwx(r,c)=xtwx(r,c)+block%last_weights(i)*block%x(i,r)*block%x(i,c)
         end do; end do
      end do
      a=xtwx+block%lambda*block%penalty
      call invert_matrix(a,ainv,istat)
      if(istat/=0)then; edf=real(p,dp); return; end if
      edf=0.0_dp
      a=matmul(ainv,xtwx)
      do i=1,p; edf=edf+a(i,i); end do
   end function block_edf

   subroutine copy_parameter(block,out)
      type(fit_block_t),intent(in)::block
      type(gamlss_parameter_result_t),intent(inout)::out
      out%coefficients=block%beta; out%eta=block%eta
      allocate(out%fitted(size(block%eta))); out%fitted=0.0_dp
      out%covariance=block%covariance
      out%lambda=block%lambda
      out%penalty=block%lambda*dot_product(block%beta,matmul(block%penalty,block%beta))
      out%edf=block_edf(block)
   end subroutine copy_parameter

   integer function total_coefficients(b,np) result(n)
      type(fit_block_t),intent(in)::b(4)
      integer,intent(in)::np
      integer::j
      n=0; do j=1,np;n=n+size(b(j)%beta);end do
   end function total_coefficients

   subroutine block_offsets(b,np,off,p)
      type(fit_block_t),intent(in)::b(4)
      integer,intent(in)::np
      integer,intent(out)::off(4),p(4)
      integer::j
      off=0;p=0
      p(1)=size(b(1)%beta); off(1)=0
      do j=2,np
         p(j)=size(b(j)%beta)
         off(j)=off(j-1)+p(j-1)
      end do
   end subroutine block_offsets

   subroutine pack_beta(b,np,theta)
      type(fit_block_t),intent(in)::b(4)
      integer,intent(in)::np
      real(dp),intent(out)::theta(:)
      integer::j,pos,p
      pos=0
      do j=1,np;p=size(b(j)%beta);theta(pos+1:pos+p)=b(j)%beta;pos=pos+p;end do
   end subroutine pack_beta

   subroutine unpack_beta(b,np,theta)
      type(fit_block_t),intent(inout)::b(4)
      integer,intent(in)::np
      real(dp),intent(in)::theta(:)
      integer::j,pos,p
      pos=0
      do j=1,np
         p=size(b(j)%beta);b(j)%beta=theta(pos+1:pos+p);pos=pos+p
         b(j)%eta=matmul(b(j)%x,b(j)%beta)+b(j)%offset
      end do
   end subroutine unpack_beta

   subroutine predict_gamlss_parameters(family,result,x_mu,mu,x_sigma,sigma,x_nu,nu,x_tau,tau, &
      offset_mu,offset_sigma,offset_nu,offset_tau,status)
      integer,intent(in)::family
      type(gamlss_result_t),intent(in)::result
      real(dp),intent(in)::x_mu(:,:)
      real(dp),allocatable,intent(out)::mu(:)
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:)
      real(dp),allocatable,intent(out),optional::sigma(:),nu(:),tau(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      integer,intent(out),optional::status
      real(dp),allocatable::e1(:),e2(:),e3(:),e4(:)
      real(dp)::a,b,c,d
      integer::n,np,i,istat
      istat=0; n=size(x_mu,1); np=family_npar(family)
      allocate(e1(n),e2(n),e3(n),e4(n)); e1=matmul(x_mu,result%mu%coefficients); e2=0;e3=0;e4=0
      if(present(offset_mu))then;if(size(offset_mu)==n)e1=e1+offset_mu;end if
      if(np>=2)then
         if(.not.present(x_sigma))then;istat=1;goto 900;end if
         e2=matmul(x_sigma,result%sigma%coefficients)
         if(present(offset_sigma))then;if(size(offset_sigma)==n)e2=e2+offset_sigma;end if
      end if
      if(np>=3)then
         if(.not.present(x_nu))then;istat=2;goto 900;end if
         e3=matmul(x_nu,result%nu%coefficients)
         if(present(offset_nu))then;if(size(offset_nu)==n)e3=e3+offset_nu;end if
      end if
      if(np>=4)then
         if(.not.present(x_tau))then;istat=3;goto 900;end if
         e4=matmul(x_tau,result%tau%coefficients)
         if(present(offset_tau))then;if(size(offset_tau)==n)e4=e4+offset_tau;end if
      end if
      allocate(mu(n)); if(present(sigma))allocate(sigma(n));if(present(nu))allocate(nu(n));if(present(tau))allocate(tau(n))
      do i=1,n
         call map_parameters(family,e1(i),e2(i),e3(i),e4(i),a,b,c,d);mu(i)=a
         if(present(sigma))sigma(i)=b;if(present(nu))nu(i)=c;if(present(tau))tau(i)=d
      end do
900   if(present(status))status=istat
   end subroutine predict_gamlss_parameters

end module gamlss_core
