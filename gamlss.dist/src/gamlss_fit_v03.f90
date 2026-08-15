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
      call bfgs_minimize(objective,th,fv,result%iterations_status,max_iter=max_iter,tol=tolerance)
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
      allocate(h(np,np)); call numerical_hessian(objective,th,h); call invert_matrix(h,cov,istat)
      if(istat==0)then; result%covariance=cov; else; allocate(result%covariance(0,0)); end if
   contains
      real(dp) function objective(par) result(nll)
         real(dp),intent(in)::par(:)
         real(dp)::mu,sig,nu,lp
         integer::j
         nll=0.0_dp
         do j=1,n
            mu=logistic(dot_product(x_mu(j,:),par(1:p1)))
            sig=exp(max(-40.0_dp,min(40.0_dp,dot_product(x_sigma(j,:),par(o2+1:o3)))))
            select case(mode)
            case(1)
               lp=dDBI(y(j),mu,sig,real(bd(j),dp),.true.)
            case(2)
               nu=logistic(dot_product(x_nu(j,:),par(o3+1:np)))
               lp=dZIBB(y(j),mu,sig,nu,bd(j),.true.)
            case default
               nu=logistic(dot_product(x_nu(j,:),par(o3+1:np)))
               lp=dZABB(y(j),mu,sig,nu,bd(j),.true.)
            end select
            if(ieee_is_finite(lp))then; nll=nll-lp; else; nll=nll+1.0e10_dp; end if
         end do
      end function objective
   end subroutine fit_fixed_binom

   elemental real(dp) function logistic(x) result(v)
      real(dp),intent(in)::x
      if(x>=0.0_dp)then; v=1.0_dp/(1.0_dp+exp(-min(x,700.0_dp))); &
      else; v=exp(max(x,-700.0_dp))/(1.0_dp+exp(max(x,-700.0_dp))); end if
   end function logistic
end module gamlss_fit_v03
