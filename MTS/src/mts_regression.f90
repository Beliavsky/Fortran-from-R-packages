! SPDX-License-Identifier: Artistic-2.0
module mts_regression
   use mts_kinds, only : dp
   use mts_types, only : varx_model, order_selection_result, mts_success, mts_invalid_input
   use mts_linalg, only : least_squares, inverse_matrix
   use mts_stats, only : determinant_ic
   use mts_var, only : var_psi_weights
   implicit none
   private

   public :: fit_varx, predict_varx, select_varx_order
   public :: multivariate_linear_model, recursive_least_squares
   public :: fit_regression_with_ar_errors

contains

   subroutine fit_varx(z,p,exog,m,model,include_mean,fixed,ridge)
      real(dp), intent(in) :: z(:,:),exog(:,:)
      integer, intent(in) :: p,m
      type(varx_model), intent(out) :: model
      logical, intent(in), optional :: include_mean
      logical, intent(in), optional :: fixed(:,:)
      real(dp), intent(in), optional :: ridge
      real(dp), allocatable :: design(:,:),y(:,:),beta(:,:),res(:,:),covb(:,:)
      real(dp), allocatable :: xm(:,:),yi(:,:),bi(:,:),ri(:,:),covi(:,:)
      logical, allocatable :: mask(:,:)
      logical :: mean_flag
      real(dp) :: s2
      integer :: n,k,kx,start,ne,d,offset,lag,j,i,ii,nfree,istat

      n=min(size(z,1),size(exog,1));k=size(z,2);kx=size(exog,2)
      mean_flag=.true.;if(present(include_mean))mean_flag=include_mean
      model%p=max(0,p);model%exog_lags=max(-1,m);model%n_series=k;model%n_exog=kx;model%include_mean=mean_flag
      start=max(p,m)+1
      if(p<0.or.m<0.or.n<start.or.k<1.or.kx<1)then
         model%status=mts_invalid_input;return
      end if
      ne=n-start+1;d=merge(1,0,mean_flag)+k*p+kx*(m+1)
      allocate(design(ne,d),y(ne,k),beta(d,k),res(ne,k),mask(d,k))
      y=z(start:n,:);offset=0
      if(mean_flag)then;design(:,1)=1.0_dp;offset=1;end if
      do lag=1,p
         design(:,offset+(lag-1)*k+1:offset+lag*k)=z(start-lag:n-lag,:)
      end do
      offset=offset+k*p
      do lag=0,m
         design(:,offset+lag*kx+1:offset+(lag+1)*kx)=exog(start-lag:n-lag,:)
      end do
      mask=.true.
      if(present(fixed))then
         if(all(shape(fixed)==shape(mask)))mask=fixed
      end if
      beta=0.0_dp;res=0.0_dp
      do j=1,k
         nfree=count(mask(:,j))
         if(nfree>0)then
            allocate(xm(ne,nfree),yi(ne,1));ii=0
            do i=1,d
               if(mask(i,j))then;ii=ii+1;xm(:,ii)=design(:,i);end if
            end do
            yi(:,1)=y(:,j)
            if(present(ridge))then
               call least_squares(xm,yi,bi,ri,covi,istat,ridge)
            else
               call least_squares(xm,yi,bi,ri,covi,istat)
            end if
            if(istat/=mts_success)then;model%status=istat;return;end if
            ii=0
            do i=1,d
               if(mask(i,j))then;ii=ii+1;beta(i,j)=bi(ii,1);end if
            end do
            res(:,j)=ri(:,1)
            deallocate(xm,yi,bi,ri,covi)
         else
            res(:,j)=y(:,j)
         end if
      end do
      model%residuals=res;model%sigma=matmul(transpose(res),res)/real(ne,dp)
      allocate(model%intercept(k),model%phi(k,k,p),model%beta(k,kx,0:m))
      model%intercept=0.0_dp;model%phi=0.0_dp;model%beta=0.0_dp;offset=0
      if(mean_flag)then;model%intercept=beta(1,:);offset=1;end if
      do lag=1,p
         model%phi(:,:,lag)=transpose(beta(offset+(lag-1)*k+1:offset+lag*k,:))
      end do
      offset=offset+k*p
      do lag=0,m
         model%beta(:,:,lag)=transpose(beta(offset+lag*kx+1:offset+(lag+1)*kx,:))
      end do
      call determinant_ic(model%sigma,ne,count(mask),model%aic,model%bic,s2,istat)
      model%status=istat
   end subroutine fit_varx

   subroutine predict_varx(model,history,exog_history,new_exog,h,forecast,forecast_covariance)
      type(varx_model), intent(in) :: model
      real(dp), intent(in) :: history(:,:),exog_history(:,:),new_exog(:,:)
      integer, intent(in) :: h
      real(dp), allocatable, intent(out) :: forecast(:,:)
      real(dp), allocatable, intent(out), optional :: forecast_covariance(:,:,:)
      real(dp), allocatable :: zwork(:,:),xwork(:,:),psi(:,:,:)
      integer :: n,nx,k,kx,t,j,idx
      n=size(history,1);nx=size(exog_history,1);k=size(history,2);kx=size(exog_history,2)
      allocate(forecast(max(0,h),k))
      if(h<=0.or.size(new_exog,1)<h.or.k/=model%n_series.or.kx/=model%n_exog)then
         if(h>0)forecast=0.0_dp;return
      end if
      allocate(zwork(n+h,k),xwork(nx+h,kx));zwork(1:n,:)=history;xwork(1:nx,:)=exog_history;xwork(nx+1:nx+h,:)=new_exog(1:h,:)
      do t=1,h
         zwork(n+t,:)=model%intercept
         do j=1,model%p
            idx=n+t-j;zwork(n+t,:)=zwork(n+t,:)+matmul(model%phi(:,:,j),zwork(idx,:))
         end do
         do j=0,model%exog_lags
            idx=nx+t-j
            if(idx>=1.and.idx<=nx+h)zwork(n+t,:)=zwork(n+t,:)+matmul(model%beta(:,:,j),xwork(idx,:))
         end do
      end do
      forecast=zwork(n+1:n+h,:)
      if(present(forecast_covariance))then
         call var_psi_weights(model%phi,h-1,psi)
         allocate(forecast_covariance(k,k,h));forecast_covariance=0.0_dp
         do t=1,h
            do j=0,t-1
               forecast_covariance(:,:,t)=forecast_covariance(:,:,t)+matmul(psi(:,:,j),matmul(model%sigma,transpose(psi(:,:,j))))
            end do
         end do
      end if
   end subroutine predict_varx

   subroutine select_varx_order(z,exog,max_p,max_m,aic,bic,hq,best_p,best_m,include_mean)
      real(dp),intent(in)::z(:,:),exog(:,:)
      integer,intent(in)::max_p,max_m
      real(dp),allocatable,intent(out)::aic(:,:),bic(:,:),hq(:,:)
      integer,intent(out)::best_p,best_m
      logical,intent(in),optional::include_mean
      type(varx_model)::model
      logical::mean_flag
      integer::p,m,loc(2)
      mean_flag=.true.;if(present(include_mean))mean_flag=include_mean
      allocate(aic(0:max_p,0:max_m),bic(0:max_p,0:max_m),hq(0:max_p,0:max_m))
      do m=0,max_m
         do p=0,max_p
            call fit_varx(z,p,exog,m,model,mean_flag)
            aic(p,m)=model%aic;bic(p,m)=model%bic
            hq(p,m)=model%aic+(2.0_dp*log(log(real(max(3,size(model%residuals,1)),dp)))-2.0_dp)* &
               real((merge(1,0,mean_flag)+size(z,2)*p+size(exog,2)*(m+1))*size(z,2),dp)/real(size(model%residuals,1),dp)
         end do
      end do
      loc=minloc(bic);best_p=loc(1)-1;best_m=loc(2)-1
   end subroutine select_varx_order

   subroutine multivariate_linear_model(y,z,coefficients,residuals,sigma,se_coefficients,include_constant,status)
      real(dp),intent(in)::y(:,:),z(:,:)
      real(dp),allocatable,intent(out)::coefficients(:,:),residuals(:,:),sigma(:,:)
      real(dp),allocatable,intent(out),optional::se_coefficients(:,:)
      logical,intent(in),optional::include_constant
      integer,intent(out),optional::status
      real(dp),allocatable::design(:,:),covb(:,:)
      logical::constant
      real(dp)::s2
      integer::i,j,istat
      constant=.true.;if(present(include_constant))constant=include_constant
      if(constant)then
         allocate(design(size(z,1),size(z,2)+1));design(:,1)=1.0_dp;design(:,2:)=z
      else;design=z;end if
      call least_squares(design,y,coefficients,residuals,covb,istat)
      sigma=matmul(transpose(residuals),residuals)/real(max(1,size(residuals,1)),dp)
      if(present(se_coefficients))then
         allocate(se_coefficients(size(coefficients,1),size(coefficients,2)))
         do j=1,size(coefficients,2)
            s2=sum(residuals(:,j)**2)/real(max(1,size(residuals,1)-size(design,2)),dp)
            do i=1,size(coefficients,1);se_coefficients(i,j)=sqrt(max(0.0_dp,covb(i,i)*s2));end do
         end do
      end if
      if(present(status))status=istat
   end subroutine multivariate_linear_model

   subroutine recursive_least_squares(y,x,initial_count,coefficients,coefficient_path,forgetting_factor,status)
      real(dp),intent(in)::y(:),x(:,:)
      integer,intent(in)::initial_count
      real(dp),allocatable,intent(out)::coefficients(:),coefficient_path(:,:)
      real(dp),intent(in),optional::forgetting_factor
      integer,intent(out),optional::status
      real(dp),allocatable::pmtx(:,:),beta0(:,:),res0(:,:),gain(:),xt(:)
      real(dp)::lambda,denom,error
      integer::n,k,t,istat
      n=size(x,1);k=size(x,2);lambda=1.0_dp;if(present(forgetting_factor))lambda=max(1.0e-6_dp,min(1.0_dp,forgetting_factor))
      if(size(y)/=n.or.initial_count<k.or.initial_count>n)then
         allocate(coefficients(k),coefficient_path(n,k));coefficients=0.0_dp;coefficient_path=0.0_dp
         if(present(status))status=mts_invalid_input;return
      end if
      call least_squares(x(1:initial_count,:),reshape(y(1:initial_count),[initial_count,1]),beta0,res0,pmtx,istat)
      allocate(coefficients(k),coefficient_path(n,k),gain(k),xt(k));coefficients=beta0(:,1);coefficient_path=0.0_dp
      coefficient_path(1:initial_count,:)=spread(coefficients,1,initial_count)
      do t=initial_count+1,n
         xt=x(t,:);denom=lambda+dot_product(xt,matmul(pmtx,xt));gain=matmul(pmtx,xt)/denom
         error=y(t)-dot_product(xt,coefficients);coefficients=coefficients+gain*error
         pmtx=(pmtx-reshape(gain,[k,1])*reshape(matmul(transpose(pmtx),xt),[1,k]))/lambda
         coefficient_path(t,:)=coefficients
      end do
      if(present(status))status=istat
   end subroutine recursive_least_squares

   subroutine fit_regression_with_ar_errors(y,x,p,coefficients,ar_model,residuals,include_mean,max_iterations,status)
      use mts_types, only : var_model
      use mts_var, only : fit_var
      real(dp),intent(in)::y(:,:),x(:,:)
      integer,intent(in)::p
      real(dp),allocatable,intent(out)::coefficients(:,:),residuals(:,:)
      type(var_model),intent(out)::ar_model
      logical,intent(in),optional::include_mean
      integer,intent(in),optional::max_iterations
      integer,intent(out),optional::status
      real(dp),allocatable::res(:,:),sigma(:,:),coef_new(:,:),se(:,:),xstar(:,:),ystar(:,:)
      integer::iter,maxit,t,j,istat,n,kx
      logical::mean_flag
      mean_flag=.true.;if(present(include_mean))mean_flag=include_mean
      maxit=20;if(present(max_iterations))maxit=max_iterations
      call multivariate_linear_model(y,x,coefficients,res,sigma,se,mean_flag,istat)
      n=size(y,1);kx=size(x,2)
      do iter=1,maxit
         call fit_var(res,p,ar_model,.false.)
         allocate(ystar(n-p,size(y,2)),xstar(n-p,kx));ystar=y(p+1:n,:);xstar=x(p+1:n,:)
         do j=1,p
            do t=1,n-p
               ystar(t,:)=ystar(t,:)-matmul(ar_model%phi(:,:,j),y(p+t-j,:))
               xstar(t,:)=xstar(t,:)-x(p+t-j,:)
            end do
         end do
         call multivariate_linear_model(ystar,xstar,coef_new,residuals,sigma,se,mean_flag,istat)
         if(maxval(abs(coef_new-coefficients))<1.0e-7_dp)then;coefficients=coef_new;exit;end if
         coefficients=coef_new
         res=y-matmul(merge_constant_design(x,mean_flag),coefficients)
         deallocate(ystar,xstar,coef_new,residuals,sigma,se)
      end do
      residuals=res
      if(present(status))status=istat
   end subroutine fit_regression_with_ar_errors

   function merge_constant_design(x,constant) result(design)
      real(dp),intent(in)::x(:,:)
      logical,intent(in)::constant
      real(dp),allocatable::design(:,:)
      if(constant)then;allocate(design(size(x,1),size(x,2)+1));design(:,1)=1.0_dp;design(:,2:)=x
      else;design=x;end if
   end function merge_constant_design

end module mts_regression
