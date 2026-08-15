! Fixed-denominator fitters for DBI/ZIBB/ZABB, which do not fit the generic four-predictor API.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_fit_v03
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gamlss_kinds, only : dp
   use gamlss_optim, only : bfgs_minimize, numerical_hessian
   use gamlss_linalg, only : invert_matrix
   use gamlss_fit, only : gamlss_fit_result_t
   use gamlss_discrete_v03, only : dDBI,dZIBB,dZABB
   implicit none
   private
   public :: fit_dbi,fit_zibb,fit_zabb

   type :: fixed_binom_context_t
      integer :: mode=0,n=0,p1=0,p2=0,p3=0,np=0,o2=0,o3=0
      integer,allocatable :: bd(:)
      real(dp),allocatable :: y(:),x_mu(:,:),x_sigma(:,:),x_nu(:,:)
   end type fixed_binom_context_t
   type(fixed_binom_context_t),save :: fixed_ctx
contains
   subroutine fit_dbi(y,bd,x_mu,x_sigma,result,start,max_iter,tol)
      real(dp),intent(in)::y(:),x_mu(:,:),x_sigma(:,:)
      integer,intent(in)::bd(:)
      type(gamlss_fit_result_t),intent(out)::result
      real(dp),intent(in),optional::start(:),tol
      integer,intent(in),optional::max_iter
      call fit_fixed_binom(1,y,bd,x_mu,x_sigma,result,start=start,max_iter=max_iter,tol=tol)
   end subroutine fit_dbi
   subroutine fit_zibb(y,bd,x_mu,x_sigma,x_nu,result,start,max_iter,tol)
      real(dp),intent(in)::y(:),x_mu(:,:),x_sigma(:,:),x_nu(:,:)
      integer,intent(in)::bd(:)
      type(gamlss_fit_result_t),intent(out)::result
      real(dp),intent(in),optional::start(:),tol
      integer,intent(in),optional::max_iter
      call fit_fixed_binom(2,y,bd,x_mu,x_sigma,result,x_nu,start,max_iter,tol)
   end subroutine fit_zibb
   subroutine fit_zabb(y,bd,x_mu,x_sigma,x_nu,result,start,max_iter,tol)
      real(dp),intent(in)::y(:),x_mu(:,:),x_sigma(:,:),x_nu(:,:)
      integer,intent(in)::bd(:)
      type(gamlss_fit_result_t),intent(out)::result
      real(dp),intent(in),optional::start(:),tol
      integer,intent(in),optional::max_iter
      call fit_fixed_binom(3,y,bd,x_mu,x_sigma,result,x_nu,start,max_iter,tol)
   end subroutine fit_zabb

   subroutine fit_fixed_binom(mode,y,bd,x_mu,x_sigma,result,x_nu,start,max_iter,tol)
      integer,intent(in)::mode,bd(:)
      real(dp),intent(in)::y(:),x_mu(:,:),x_sigma(:,:)
      type(gamlss_fit_result_t),intent(out)::result
      real(dp),intent(in),optional::x_nu(:,:),start(:),tol
      integer,intent(in),optional::max_iter
      integer::n,p1,p2,p3,np,o2,o3,i,istat
      real(dp),allocatable::th(:),h(:,:),cov(:,:)
      real(dp)::fv,tolerance,yy,bb,m0
      n=size(y); p1=size(x_mu,2); p2=size(x_sigma,2); p3=0
      if(mode>1)then
         if(.not.present(x_nu))then; result%status=21; return; end if
         p3=size(x_nu,2)
      end if
      if(size(bd)/=n.or.size(x_mu,1)/=n.or.size(x_sigma,1)/=n)then; result%status=22; return; end if
      if(mode>1)then; if(size(x_nu,1)/=n)then; result%status=23; return; end if; end if
      np=p1+p2+p3; o2=p1; o3=p1+p2; allocate(th(np)); th=0.0_dp
      if(present(start))then
         if(size(start)/=np)then; result%status=24; return; end if
         th=start
      else
         yy=0.0_dp; bb=0.0_dp
         do i=1,n; yy=yy+y(i); bb=bb+real(bd(i),dp); end do
         m0=min(0.95_dp,max(0.05_dp,yy/max(bb,1.0_dp)))
         if(p1>0)th(1)=log(m0/(1.0_dp-m0))
         if(p2>0)th(o2+1)=log(0.8_dp)
         if(p3>0)th(o3+1)=log(0.1_dp/0.9_dp)
      end if
      tolerance=1.0e-7_dp; if(present(tol))tolerance=tol
      call set_fixed_context(mode,y,bd,x_mu,x_sigma,p1,p2,p3,x_nu)
      call bfgs_minimize(fixed_binom_objective,th,fv,result%iterations_status,max_iter=max_iter,tol=tolerance)
      result%status=result%iterations_status; result%converged=(result%status==0)
      result%loglik=-fv; result%aic=2.0_dp*real(np,dp)-2.0_dp*result%loglik
      result%beta_mu=th(1:p1); result%beta_sigma=th(o2+1:o3)
      if(p3>0)result%beta_nu=th(o3+1:np)
      allocate(result%fitted_mu(n),result%fitted_sigma(n)); if(p3>0)allocate(result%fitted_nu(n))
      do i=1,n
         result%fitted_mu(i)=logistic(dot_product(x_mu(i,:),th(1:p1)))
         result%fitted_sigma(i)=exp(max(-40.0_dp,min(40.0_dp,dot_product(x_sigma(i,:),th(o2+1:o3)))))
         if(p3>0)result%fitted_nu(i)=logistic(dot_product(x_nu(i,:),th(o3+1:np)))
      end do
      allocate(h(np,np)); call numerical_hessian(fixed_binom_objective,th,h); call invert_matrix(h,cov,istat)
      if(istat==0)then; result%covariance=cov; else; allocate(result%covariance(0,0)); end if
      call clear_fixed_context()
   end subroutine fit_fixed_binom

   subroutine set_fixed_context(mode,y,bd,x_mu,x_sigma,p1,p2,p3,x_nu)
      integer,intent(in)::mode,bd(:),p1,p2,p3
      real(dp),intent(in)::y(:),x_mu(:,:),x_sigma(:,:)
      real(dp),intent(in),optional::x_nu(:,:)
      call clear_fixed_context()
      fixed_ctx%mode=mode;fixed_ctx%n=size(y);fixed_ctx%p1=p1;fixed_ctx%p2=p2;fixed_ctx%p3=p3
      fixed_ctx%o2=p1;fixed_ctx%o3=p1+p2;fixed_ctx%np=p1+p2+p3
      fixed_ctx%y=y;fixed_ctx%bd=bd;fixed_ctx%x_mu=x_mu;fixed_ctx%x_sigma=x_sigma
      if(present(x_nu))then;fixed_ctx%x_nu=x_nu;else;allocate(fixed_ctx%x_nu(size(y),0));end if
   end subroutine set_fixed_context

   subroutine clear_fixed_context()
      if(allocated(fixed_ctx%bd))deallocate(fixed_ctx%bd)
      if(allocated(fixed_ctx%y))deallocate(fixed_ctx%y)
      if(allocated(fixed_ctx%x_mu))deallocate(fixed_ctx%x_mu)
      if(allocated(fixed_ctx%x_sigma))deallocate(fixed_ctx%x_sigma)
      if(allocated(fixed_ctx%x_nu))deallocate(fixed_ctx%x_nu)
      fixed_ctx%mode=0;fixed_ctx%n=0
   end subroutine clear_fixed_context

   real(dp) function fixed_binom_objective(par) result(nll)
      real(dp),intent(in)::par(:)
      real(dp)::mu,sig,nu,lp
      integer::j
      nll=0.0_dp
      do j=1,fixed_ctx%n
         mu=logistic(dot_product(fixed_ctx%x_mu(j,:),par(1:fixed_ctx%p1)))
         sig=exp(max(-40.0_dp,min(40.0_dp,dot_product(fixed_ctx%x_sigma(j,:), &
            par(fixed_ctx%o2+1:fixed_ctx%o3)))))
         select case(fixed_ctx%mode)
         case(1)
            lp=dDBI(fixed_ctx%y(j),mu,sig,real(fixed_ctx%bd(j),dp),.true.)
         case(2)
            nu=logistic(dot_product(fixed_ctx%x_nu(j,:),par(fixed_ctx%o3+1:fixed_ctx%np)))
            lp=dZIBB(fixed_ctx%y(j),mu,sig,nu,fixed_ctx%bd(j),.true.)
         case default
            nu=logistic(dot_product(fixed_ctx%x_nu(j,:),par(fixed_ctx%o3+1:fixed_ctx%np)))
            lp=dZABB(fixed_ctx%y(j),mu,sig,nu,fixed_ctx%bd(j),.true.)
         end select
         if(ieee_is_finite(lp))then;nll=nll-lp;else;nll=nll+1.0e10_dp;end if
      end do
   end function fixed_binom_objective

   elemental real(dp) function logistic(x) result(v)
      real(dp),intent(in)::x
      if(x>=0.0_dp)then; v=1.0_dp/(1.0_dp+exp(-min(x,700.0_dp))); &
      else; v=exp(max(x,-700.0_dp))/(1.0_dp+exp(max(x,-700.0_dp))); end if
   end function logistic
end module gamlss_fit_v03
