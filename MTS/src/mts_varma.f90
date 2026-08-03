! SPDX-License-Identifier: Artistic-2.0
module mts_varma
   use mts_kinds, only : dp
   use mts_types, only : varma_model, mts_success, mts_invalid_input
   use mts_linalg, only : least_squares, log_determinant, eye, matrix_sqrt_symmetric
   use mts_stats, only : column_mean, covariance_matrix
   use mts_rng, only : random_multivariate_normal
   use mts_optimize, only : bfgs_minimize
   use mts_var, only : fit_var
   use mts_types, only : var_model
   implicit none
   private

   public :: varma_residuals, fit_varma, fit_vma
   public :: simulate_varma, varma_psi_weights, predict_varma
   public :: varma_impulse_response, varma_covariance
   public :: difference_matrix, seasonal_lag_matrices

contains

   subroutine varma_residuals(x,intercept,phi,theta,residuals,start_index)
      real(dp), intent(in) :: x(:,:),intercept(:),phi(:,:,:),theta(:,:,:)
      real(dp), intent(out) :: residuals(size(x,1),size(x,2))
      integer, intent(out), optional :: start_index
      integer :: n,k,p,q,t,j,start
      n=size(x,1);k=size(x,2);p=size(phi,3);q=size(theta,3)
      residuals=0.0_dp
      start=max(p,q)+1
      do t=1,n
         residuals(t,:)=x(t,:)-intercept
         do j=1,min(p,t-1)
            residuals(t,:)=residuals(t,:)-matmul(phi(:,:,j),x(t-j,:))
         end do
         do j=1,min(q,t-1)
            residuals(t,:)=residuals(t,:)+matmul(theta(:,:,j),residuals(t-j,:))
         end do
      end do
      if(present(start_index))start_index=start
   end subroutine varma_residuals

   subroutine fit_vma(x,q,model,include_mean,max_iterations,tolerance)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: q
      type(varma_model), intent(out) :: model
      logical, intent(in), optional :: include_mean
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      call fit_varma(x,0,q,model,include_mean,max_iterations,tolerance)
   end subroutine fit_vma

   subroutine fit_varma(x,p,q,model,include_mean,max_iterations,tolerance)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: p,q
      type(varma_model), intent(out) :: model
      logical, intent(in), optional :: include_mean
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: par(:),lo(:),hi(:),res(:,:),sigma(:,:)
      real(dp), allocatable :: intercept(:),phi(:,:,:),theta(:,:,:)
      real(dp) :: objective_value,tol
      integer :: n,k,np,iterations,istat,maxit,start
      logical :: mean_flag

      n=size(x,1);k=size(x,2)
      mean_flag=.true.;if(present(include_mean))mean_flag=include_mean
      model%p=max(0,p);model%q=max(0,q);model%n_series=k;model%include_mean=mean_flag
      if(p<0.or.q<0.or.n<=max(p,q)+2.or.k<1)then
         model%status=mts_invalid_input
         return
      end if
      np=merge(k,0,mean_flag)+k*k*(p+q)
      allocate(par(np),lo(np),hi(np))
      call initial_varma_parameters(x,p,q,mean_flag,par)
      lo=-2.5_dp;hi=2.5_dp
      if(mean_flag)then
         lo(1:k)=minval(x,dim=1)-5.0_dp*sqrt(max(1.0e-12_dp,diagonal(covariance_matrix(x))))
         hi(1:k)=maxval(x,dim=1)+5.0_dp*sqrt(max(1.0e-12_dp,diagonal(covariance_matrix(x))))
      end if
      maxit=300;if(present(max_iterations))maxit=max_iterations
      tol=1.0e-7_dp;if(present(tolerance))tol=tolerance
      call bfgs_minimize(objective,par,objective_value,istat,iterations,maxit,tol,lo,hi)
      call unpack_parameters(par,k,p,q,mean_flag,intercept,phi,theta)
      allocate(res(n,k))
      call varma_residuals(x,intercept,phi,theta,res,start)
      sigma=matmul(transpose(res(start:n,:)),res(start:n,:))/real(n-start+1,dp)
      model%intercept=intercept;model%phi=phi;model%theta=theta
      model%residuals=res(start:n,:);model%sigma=sigma
      model%log_likelihood=-objective_value
      model%aic=2.0_dp*objective_value+2.0_dp*real(np+k*(k+1)/2,dp)
      model%bic=2.0_dp*objective_value+log(real(n-start+1,dp))*real(np+k*(k+1)/2,dp)
      model%iterations=iterations;model%status=istat

   contains
      function objective(parameters) result(value)
         real(dp), intent(in) :: parameters(:)
         real(dp) :: value,logdet,quad
         real(dp), allocatable :: c(:),ph(:,:,:),th(:,:,:),e(:,:),s(:,:),sinv(:,:)
         integer :: t,st,sgn,status_local
         call unpack_parameters(parameters,k,p,q,mean_flag,c,ph,th)
         allocate(e(n,k));call varma_residuals(x,c,ph,th,e,st)
         s=matmul(transpose(e(st:n,:)),e(st:n,:))/real(n-st+1,dp)
         logdet=log_determinant(s,sgn,status_local)
         if(status_local/=mts_success.or.sgn<=0)then
            value=huge(1.0_dp)/100.0_dp;return
         end if
         call inverse_local(s,sinv,status_local)
         if(status_local/=mts_success)then
            value=huge(1.0_dp)/100.0_dp;return
         end if
         quad=0.0_dp
         do t=st,n
            quad=quad+dot_product(e(t,:),matmul(sinv,e(t,:)))
         end do
         value=0.5_dp*(real(n-st+1,dp)*(real(k,dp)*log(2.0_dp*acos(-1.0_dp))+logdet)+quad)
         value=value+stability_penalty(ph,th)
      end function objective
   end subroutine fit_varma

   subroutine inverse_local(a,ainv,status)
      use mts_linalg, only : inverse_matrix
      real(dp),intent(in)::a(:,:)
      real(dp),allocatable,intent(out)::ainv(:,:)
      integer,intent(out)::status
      call inverse_matrix(a,ainv,status)
   end subroutine inverse_local

   subroutine initial_varma_parameters(x,p,q,include_mean,par)
      real(dp),intent(in)::x(:,:)
      integer,intent(in)::p,q
      logical,intent(in)::include_mean
      real(dp),intent(out)::par(:)
      type(var_model)::prelim
      real(dp),allocatable::design(:,:),y(:,:),beta(:,:),resi(:,:),e0(:,:)
      integer::n,k,r,t,j,offset,idx,istat,preorder
      n=size(x,1);k=size(x,2);par=0.0_dp;offset=0
      if(include_mean)then
         par(1:k)=column_mean(x);offset=k
      end if
      preorder=min(max(p,q)+5,max(1,n/10))
      call fit_var(x,preorder,prelim,include_mean)
      if(prelim%status/=mts_success)then
         if(p>0)then
            call fit_var(x,p,prelim,include_mean)
            do j=1,p
               par(offset+(j-1)*k*k+1:offset+j*k*k)=reshape(prelim%phi(:,:,j),[k*k])
            end do
         end if
         return
      end if
      allocate(e0(n,k)); e0=0.0_dp
      e0(preorder+1:n,:)=prelim%residuals
      r=max(max(p,q),preorder)+1
      allocate(design(n-r+1,merge(1,0,include_mean)+k*(p+q)),y(n-r+1,k))
      idx=0
      if(include_mean)then;idx=1;design(:,1)=1.0_dp;end if
      do j=1,p
         design(:,idx+(j-1)*k+1:idx+j*k)=x(r-j:n-j,:)
      end do
      idx=idx+k*p
      do j=1,q
         design(:,idx+(j-1)*k+1:idx+j*k)=-e0(r-j:n-j,:)
      end do
      y=x(r:n,:)
      call least_squares(design,y,beta,resi,status=istat,ridge=1.0e-8_dp)
      if(istat==mts_success)then
         idx=0
         if(include_mean)then
            par(1:k)=beta(1,:);offset=k;idx=1
         end if
         do j=1,p
            par(offset+(j-1)*k*k+1:offset+j*k*k)=reshape(transpose(beta(idx+(j-1)*k+1:idx+j*k,:)),[k*k])
         end do
         offset=offset+k*k*p;idx=idx+k*p
         do j=1,q
            par(offset+(j-1)*k*k+1:offset+j*k*k)=reshape(transpose(beta(idx+(j-1)*k+1:idx+j*k,:)),[k*k])
         end do
      end if
   end subroutine initial_varma_parameters

   subroutine unpack_parameters(par,k,p,q,include_mean,intercept,phi,theta)
      real(dp),intent(in)::par(:)
      integer,intent(in)::k,p,q
      logical,intent(in)::include_mean
      real(dp),allocatable,intent(out)::intercept(:),phi(:,:,:),theta(:,:,:)
      integer::offset,j
      allocate(intercept(k),phi(k,k,p),theta(k,k,q));intercept=0.0_dp;phi=0.0_dp;theta=0.0_dp
      offset=0
      if(include_mean)then;intercept=par(1:k);offset=k;end if
      do j=1,p
         phi(:,:,j)=reshape(par(offset+(j-1)*k*k+1:offset+j*k*k),[k,k])
      end do
      offset=offset+k*k*p
      do j=1,q
         theta(:,:,j)=reshape(par(offset+(j-1)*k*k+1:offset+j*k*k),[k,k])
      end do
   end subroutine unpack_parameters

   function stability_penalty(phi,theta) result(penalty)
      real(dp),intent(in)::phi(:,:,:),theta(:,:,:)
      real(dp)::penalty,normp,normq
      integer::j
      normp=0.0_dp;normq=0.0_dp
      do j=1,size(phi,3);normp=normp+sqrt(sum(phi(:,:,j)**2));end do
      do j=1,size(theta,3);normq=normq+sqrt(sum(theta(:,:,j)**2));end do
      penalty=0.0_dp
      if(normp>2.5_dp)penalty=penalty+1.0e4_dp*(normp-2.5_dp)**2
      if(normq>2.5_dp)penalty=penalty+1.0e4_dp*(normq-2.5_dp)**2
   end function stability_penalty

   subroutine simulate_varma(n,intercept,phi,theta,sigma,x,innovations,burn_in)
      integer,intent(in)::n
      real(dp),intent(in)::intercept(:),phi(:,:,:),theta(:,:,:),sigma(:,:)
      real(dp),intent(out)::x(n,size(intercept))
      real(dp),intent(out),optional::innovations(n,size(intercept))
      integer,intent(in),optional::burn_in
      real(dp),allocatable::work(:,:),noise(:,:)
      real(dp)::z(size(intercept))
      integer::burn,nt,k,p,q,start,t,j
      burn=200;if(present(burn_in))burn=max(0,burn_in)
      k=size(intercept);p=size(phi,3);q=size(theta,3);start=max(p,q)+1;nt=n+burn
      allocate(work(nt+start,k),noise(nt+start,k));work=0.0_dp;noise=0.0_dp
      do t=start,nt+start-1
         call random_multivariate_normal(0.0_dp*intercept,sigma,z)
         noise(t,:)=z;work(t,:)=intercept+z
         do j=1,p;work(t,:)=work(t,:)+matmul(phi(:,:,j),work(t-j,:));end do
         do j=1,q;work(t,:)=work(t,:)-matmul(theta(:,:,j),noise(t-j,:));end do
      end do
      x=work(start+burn:start+burn+n-1,:)
      if(present(innovations))innovations=noise(start+burn:start+burn+n-1,:)
   end subroutine simulate_varma

   subroutine varma_psi_weights(phi,theta,max_lag,psi)
      real(dp),intent(in)::phi(:,:,:),theta(:,:,:)
      integer,intent(in)::max_lag
      real(dp),allocatable,intent(out)::psi(:,:,:)
      integer::k,p,q,h,j
      k=size(phi,1);p=size(phi,3);q=size(theta,3)
      if(k==0.and.size(theta,1)>0)k=size(theta,1)
      allocate(psi(k,k,0:max_lag));psi=0.0_dp;psi(:,:,0)=eye(k)
      do h=1,max_lag
         do j=1,min(p,h);psi(:,:,h)=psi(:,:,h)+matmul(phi(:,:,j),psi(:,:,h-j));end do
         if(h<=q)psi(:,:,h)=psi(:,:,h)-theta(:,:,h)
      end do
   end subroutine varma_psi_weights

   subroutine predict_varma(model,history,h,forecast,forecast_covariance)
      type(varma_model),intent(in)::model
      real(dp),intent(in)::history(:,:)
      integer,intent(in)::h
      real(dp),allocatable,intent(out)::forecast(:,:)
      real(dp),allocatable,intent(out),optional::forecast_covariance(:,:,:)
      real(dp),allocatable::work(:,:),res(:,:),psi(:,:,:)
      integer::n,k,t,j,idx
      n=size(history,1);k=size(history,2)
      allocate(forecast(max(0,h),k))
      if(h<=0.or.n<max(model%p,model%q))then;if(h>0)forecast=0.0_dp;return;end if
      allocate(work(n+h,k),res(n+h,k));work(1:n,:)=history
      call varma_residuals(history,model%intercept,model%phi,model%theta,res(1:n,:))
      res(n+1:n+h,:)=0.0_dp
      do t=1,h
         work(n+t,:)=model%intercept
         do j=1,model%p
            idx=n+t-j;work(n+t,:)=work(n+t,:)+matmul(model%phi(:,:,j),work(idx,:))
         end do
         do j=1,model%q
            idx=n+t-j;if(idx<=n)work(n+t,:)=work(n+t,:)-matmul(model%theta(:,:,j),res(idx,:))
         end do
      end do
      forecast=work(n+1:n+h,:)
      if(present(forecast_covariance))then
         call varma_psi_weights(model%phi,model%theta,h-1,psi)
         allocate(forecast_covariance(k,k,h));forecast_covariance=0.0_dp
         do t=1,h
            do j=0,t-1
               forecast_covariance(:,:,t)=forecast_covariance(:,:,t)+matmul(psi(:,:,j),matmul(model%sigma,transpose(psi(:,:,j))))
            end do
         end do
      end if
   end subroutine predict_varma

   subroutine varma_impulse_response(phi,theta,sigma,max_lag,response,orthogonal)
      real(dp),intent(in)::phi(:,:,:),theta(:,:,:),sigma(:,:)
      integer,intent(in)::max_lag
      real(dp),allocatable,intent(out)::response(:,:,:)
      logical,intent(in),optional::orthogonal
      real(dp),allocatable::psi(:,:,:),impact(:,:)
      logical::orth
      integer::h,istat,k
      orth=.true.;if(present(orthogonal))orth=orthogonal;k=size(sigma,1)
      call varma_psi_weights(phi,theta,max_lag,psi)
      if(orth)then
         call matrix_sqrt_symmetric(sigma,impact,istat);if(istat/=mts_success)impact=eye(k)
      else;impact=eye(k);end if
      allocate(response(k,k,0:max_lag))
      do h=0,max_lag;response(:,:,h)=matmul(psi(:,:,h),impact);end do
   end subroutine varma_impulse_response

   subroutine varma_covariance(phi,theta,sigma,max_lag,covariances,truncation)
      real(dp),intent(in)::phi(:,:,:),theta(:,:,:),sigma(:,:)
      integer,intent(in)::max_lag
      real(dp),allocatable,intent(out)::covariances(:,:,:)
      integer,intent(in),optional::truncation
      real(dp),allocatable::psi(:,:,:)
      integer::trun,h,j,k
      k=size(sigma,1);trun=max(100,max_lag+50);if(present(truncation))trun=max(max_lag,truncation)
      call varma_psi_weights(phi,theta,trun,psi)
      allocate(covariances(k,k,0:max_lag));covariances=0.0_dp
      do h=0,max_lag
         do j=0,trun-h
            covariances(:,:,h)=covariances(:,:,h)+matmul(psi(:,:,j+h),matmul(sigma,transpose(psi(:,:,j))))
         end do
      end do
   end subroutine varma_covariance

   function difference_matrix(x,d,seasonal_lag,seasonal_difference) result(y)
      real(dp),intent(in)::x(:,:)
      integer,intent(in),optional::d,seasonal_lag,seasonal_difference
      real(dp),allocatable::y(:,:),tmp(:,:)
      integer::nd,sd,s,i
      tmp=x;nd=0;if(present(d))nd=max(0,d)
      do i=1,nd
         if(size(tmp,1)<=1)then;allocate(y(0,size(x,2)));return;end if
         tmp=tmp(2:size(tmp,1),:)-tmp(1:size(tmp,1)-1,:)
      end do
      sd=0;if(present(seasonal_difference))sd=max(0,seasonal_difference)
      s=1;if(present(seasonal_lag))s=max(1,seasonal_lag)
      do i=1,sd
         if(size(tmp,1)<=s)then;allocate(y(0,size(x,2)));return;end if
         tmp=tmp(s+1:size(tmp,1),:)-tmp(1:size(tmp,1)-s,:)
      end do
      y=tmp
   end function difference_matrix

   subroutine seasonal_lag_matrices(nonseasonal,seasonal,season,combined)
      real(dp),intent(in)::nonseasonal(:,:,:),seasonal(:,:,:)
      integer,intent(in)::season
      real(dp),allocatable,intent(out)::combined(:,:,:)
      integer::k,n_nonseasonal,n_seasonal,i,j,maxlag
      k=size(nonseasonal,1);n_nonseasonal=size(nonseasonal,3);n_seasonal=size(seasonal,3)
      maxlag=n_nonseasonal+season*n_seasonal
      allocate(combined(k,k,maxlag));combined=0.0_dp
      do i=1,n_nonseasonal;combined(:,:,i)=combined(:,:,i)+nonseasonal(:,:,i);end do
      do j=1,n_seasonal;combined(:,:,season*j)=combined(:,:,season*j)+seasonal(:,:,j);end do
      do i=1,n_nonseasonal
         do j=1,n_seasonal
            combined(:,:,i+season*j)=combined(:,:,i+season*j)-matmul(nonseasonal(:,:,i),seasonal(:,:,j))
         end do
      end do
   end subroutine seasonal_lag_matrices

   pure function diagonal(a) result(d)
      real(dp),intent(in)::a(:,:)
      real(dp)::d(min(size(a,1),size(a,2)))
      integer::i
      do i=1,size(d);d(i)=a(i,i);end do
   end function diagonal

end module mts_varma
