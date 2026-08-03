! SPDX-License-Identifier: Artistic-2.0
module mts_var
   use mts_kinds, only : dp
   use mts_types, only : var_model, order_selection_result, diagnostic_result, &
      mts_success, mts_invalid_input
   use mts_linalg, only : least_squares, inverse_matrix, eye, matrix_sqrt_symmetric
   use mts_stats, only : covariance_matrix, determinant_ic, chi_square_survival
   use mts_rng, only : random_multivariate_normal
   implicit none
   private

   public :: fit_var, fit_sparse_var, refine_var
   public :: select_var_order, select_var_order_increasing
   public :: var_psi_weights, predict_var, simulate_var
   public :: var_impulse_response, forecast_error_variance_decomposition
   public :: granger_causality_test

contains

   subroutine fit_var(x,p,model,include_mean,fixed,ridge)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: p
      type(var_model), intent(out) :: model
      logical, intent(in), optional :: include_mean
      logical, intent(in), optional :: fixed(:,:)
      real(dp), intent(in), optional :: ridge
      real(dp), allocatable :: design(:,:), y(:,:), beta(:,:), residuals(:,:), xtx_inv(:,:)
      real(dp), allocatable :: xm(:,:), yi(:,:), bi(:,:), ri(:,:), covi(:,:)
      logical, allocatable :: mask(:,:)
      logical :: mean_flag
      real(dp) :: s2
      integer :: n,k,ne,d,offset,lag,i,j,ii,nfree,istat

      n=size(x,1); k=size(x,2)
      mean_flag=.true.
      if (present(include_mean)) mean_flag=include_mean
      model%order=max(0,p); model%n_series=k; model%include_mean=mean_flag
      if (p < 0 .or. n <= p .or. k < 1) then
         model%status=mts_invalid_input
         return
      end if
      ne=n-p
      d=k*p+merge(1,0,mean_flag)
      allocate(design(ne,d),y(ne,k),beta(d,k),residuals(ne,k),mask(d,k))
      y=x(p+1:n,:)
      offset=0
      if (mean_flag) then
         design(:,1)=1.0_dp
         offset=1
      end if
      do lag=1,p
         design(:,offset+(lag-1)*k+1:offset+lag*k)=x(p+1-lag:n-lag,:)
      end do
      mask=.true.
      if (present(fixed)) then
         if (size(fixed,1)==d .and. size(fixed,2)==k) mask=fixed
      end if
      beta=0.0_dp; residuals=0.0_dp; nfree=0
      allocate(model%se_coef(d,k)); model%se_coef=0.0_dp
      do j=1,k
         nfree=count(mask(:,j))
         if (nfree>0) then
            allocate(xm(ne,nfree),yi(ne,1))
            ii=0
            do i=1,d
               if (mask(i,j)) then
                  ii=ii+1
                  xm(:,ii)=design(:,i)
               end if
            end do
            yi(:,1)=y(:,j)
            if (present(ridge)) then
               call least_squares(xm,yi,bi,ri,covi,istat,ridge)
            else
               call least_squares(xm,yi,bi,ri,covi,istat)
            end if
            if (istat/=mts_success) then
               model%status=istat
               return
            end if
            ii=0
            do i=1,d
               if (mask(i,j)) then
                  ii=ii+1
                  beta(i,j)=bi(ii,1)
               end if
            end do
            residuals(:,j)=ri(:,1)
            s2=sum(ri(:,1)**2)/real(max(1,ne-nfree),dp)
            ii=0
            do i=1,d
               if (mask(i,j)) then
                  ii=ii+1
                  model%se_coef(i,j)=sqrt(max(0.0_dp,covi(ii,ii)*s2))
               end if
            end do
            deallocate(xm,yi,bi,ri,covi)
         else
            residuals(:,j)=y(:,j)
         end if
      end do
      model%coef=beta
      model%residuals=residuals
      model%sigma=matmul(transpose(residuals),residuals)/real(ne,dp)
      allocate(model%intercept(k),model%phi(k,k,p))
      model%intercept=0.0_dp
      model%phi=0.0_dp
      offset=0
      if (mean_flag) then
         model%intercept=beta(1,:)
         offset=1
      end if
      do lag=1,p
         model%phi(:,:,lag)=transpose(beta(offset+(lag-1)*k+1:offset+lag*k,:))
      end do
      model%n_parameters=count(mask)
      call determinant_ic(model%sigma,ne,model%n_parameters,model%aic,model%bic,model%hq,istat)
      model%status=istat
   end subroutine fit_var

   subroutine fit_sparse_var(x,lags,model,include_mean,fixed,ridge)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: lags(:)
      type(var_model), intent(out) :: model
      logical, intent(in), optional :: include_mean
      logical, intent(in), optional :: fixed(:,:)
      real(dp), intent(in), optional :: ridge
      logical, allocatable :: mask(:,:)
      integer :: p,k,d,offset,lag,j
      if (size(lags)<1) then
         if (present(ridge)) then
            call fit_var(x,0,model,include_mean,ridge=ridge)
         else if (present(include_mean)) then
            call fit_var(x,0,model,include_mean)
         else
            call fit_var(x,0,model)
         end if
         return
      end if
      p=maxval(lags);k=size(x,2)
      d=k*p+merge(1,0,present_and_true(include_mean))
      allocate(mask(d,k));mask=.false.
      offset=merge(1,0,present_and_true(include_mean))
      if (offset==1) mask(1,:)=.true.
      do j=1,size(lags)
         lag=lags(j)
         if (lag>=1 .and. lag<=p) mask(offset+(lag-1)*k+1:offset+lag*k,:)=.true.
      end do
      if (present(fixed)) then
         if (all(shape(fixed)==shape(mask))) mask=mask.and.fixed
      end if
      if (present(ridge)) then
         call fit_var(x,p,model,present_and_true(include_mean),mask,ridge)
      else
         call fit_var(x,p,model,present_and_true(include_mean),mask)
      end if
   end subroutine fit_sparse_var

   subroutine refine_var(x,initial_model,model,threshold)
      real(dp), intent(in) :: x(:,:)
      type(var_model), intent(in) :: initial_model
      type(var_model), intent(out) :: model
      real(dp), intent(in), optional :: threshold
      logical, allocatable :: mask(:,:)
      real(dp) :: thres
      integer :: i,j
      thres=1.0_dp
      if (present(threshold)) thres=max(0.0_dp,threshold)
      if (.not. allocated(initial_model%coef) .or. .not. allocated(initial_model%se_coef)) then
         model=initial_model
         model%status=mts_invalid_input
         return
      end if
      allocate(mask(size(initial_model%coef,1),size(initial_model%coef,2)))
      mask=.false.
      do j=1,size(mask,2)
         do i=1,size(mask,1)
            if (initial_model%se_coef(i,j)>0.0_dp) &
               mask(i,j)=abs(initial_model%coef(i,j)/initial_model%se_coef(i,j))>=thres
         end do
      end do
      call fit_var(x,initial_model%order,model,initial_model%include_mean,mask)
   end subroutine refine_var

   subroutine select_var_order(x,max_order,result,same_sample,include_mean)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: max_order
      type(order_selection_result), intent(out) :: result
      logical, intent(in), optional :: same_sample, include_mean
      logical :: common, mean_flag
      type(var_model) :: model
      real(dp), allocatable :: sample(:,:)
      real(dp) :: previous_logdet,current_logdet
      integer :: p,n,k,ne,df

      n=size(x,1);k=size(x,2)
      common=.true.;mean_flag=.true.
      if (present(same_sample)) common=same_sample
      if (present(include_mean)) mean_flag=include_mean
      result%max_order=max_order
      if (max_order<0 .or. n<=max_order+1) then
         result%status=mts_invalid_input
         return
      end if
      allocate(result%aic(0:max_order),result%bic(0:max_order),result%hq(0:max_order))
      allocate(result%m_stat(1:max_order),result%p_value(1:max_order))
      result%m_stat=0.0_dp;result%p_value=1.0_dp
      previous_logdet=0.0_dp
      do p=0,max_order
         if (common) then
            sample=x(max_order-p+1:n,:)
         else
            sample=x
         end if
         call fit_var(sample,p,model,mean_flag)
         result%aic(p)=model%aic
         result%bic(p)=model%bic
         result%hq(p)=model%hq
         current_logdet=model%aic-2.0_dp*real(model%n_parameters,dp)/real(size(model%residuals,1),dp)
         if (p>0) then
            ne=n-max_order
            df=k*k
            result%m_stat(p)=max(0.0_dp,(real(ne-k*p,dp)-1.5_dp)*(previous_logdet-current_logdet))
            result%p_value(p)=chi_square_survival(result%m_stat(p),df)
         end if
         previous_logdet=current_logdet
      end do
      result%aic_order=minloc(result%aic,dim=1)-1
      result%bic_order=minloc(result%bic,dim=1)-1
      result%hq_order=minloc(result%hq,dim=1)-1
      result%status=mts_success
   end subroutine select_var_order

   subroutine select_var_order_increasing(x,max_order,result,include_mean)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: max_order
      type(order_selection_result), intent(out) :: result
      logical, intent(in), optional :: include_mean
      call select_var_order(x,max_order,result,.false.,include_mean)
   end subroutine select_var_order_increasing

   subroutine var_psi_weights(phi,max_lag,psi)
      real(dp), intent(in) :: phi(:,:,:)
      integer, intent(in) :: max_lag
      real(dp), allocatable, intent(out) :: psi(:,:,:)
      integer :: k,p,h,j
      k=size(phi,1);p=size(phi,3)
      allocate(psi(k,k,0:max_lag));psi=0.0_dp
      psi(:,:,0)=eye(k)
      do h=1,max_lag
         do j=1,min(p,h)
            psi(:,:,h)=psi(:,:,h)+matmul(phi(:,:,j),psi(:,:,h-j))
         end do
      end do
   end subroutine var_psi_weights

   subroutine predict_var(model,history,h,forecast,forecast_covariance)
      type(var_model), intent(in) :: model
      real(dp), intent(in) :: history(:,:)
      integer, intent(in) :: h
      real(dp), allocatable, intent(out) :: forecast(:,:)
      real(dp), allocatable, intent(out), optional :: forecast_covariance(:,:,:)
      real(dp), allocatable :: work(:,:),psi(:,:,:)
      integer :: n,k,p,t,j,idx
      n=size(history,1);k=size(history,2);p=model%order
      allocate(forecast(max(0,h),k))
      if (h<=0 .or. k/=model%n_series .or. n<p) then
         if (h>0) forecast=0.0_dp
         if (present(forecast_covariance)) allocate(forecast_covariance(k,k,0))
         return
      end if
      allocate(work(n+h,k));work(1:n,:)=history
      do t=1,h
         work(n+t,:)=model%intercept
         do j=1,p
            idx=n+t-j
            work(n+t,:)=work(n+t,:)+matmul(model%phi(:,:,j),work(idx,:))
         end do
      end do
      forecast=work(n+1:n+h,:)
      if (present(forecast_covariance)) then
         allocate(forecast_covariance(k,k,h));forecast_covariance=0.0_dp
         call var_psi_weights(model%phi,h-1,psi)
         do t=1,h
            do j=0,t-1
               forecast_covariance(:,:,t)=forecast_covariance(:,:,t)+ &
                  matmul(psi(:,:,j),matmul(model%sigma,transpose(psi(:,:,j))))
            end do
         end do
      end if
   end subroutine predict_var

   subroutine simulate_var(n,intercept,phi,sigma,x,burn_in,initial)
      integer, intent(in) :: n
      real(dp), intent(in) :: intercept(:),phi(:,:,:),sigma(:,:)
      real(dp), intent(out) :: x(n,size(intercept))
      integer, intent(in), optional :: burn_in
      real(dp), intent(in), optional :: initial(:,:)
      real(dp), allocatable :: work(:,:)
      real(dp) :: innovation(size(intercept))
      integer :: burn,nt,k,p,t,j,start
      burn=200
      if (present(burn_in)) burn=max(0,burn_in)
      k=size(intercept);p=size(phi,3);nt=n+burn
      allocate(work(nt+p,k));work=0.0_dp
      if (present(initial)) then
         start=min(p,size(initial,1))
         if (start>0) work(p-start+1:p,:)=initial(size(initial,1)-start+1:size(initial,1),:)
      end if
      do t=p+1,nt+p
         work(t,:)=intercept
         do j=1,p
            work(t,:)=work(t,:)+matmul(phi(:,:,j),work(t-j,:))
         end do
         call random_multivariate_normal(0.0_dp*intercept,sigma,innovation)
         work(t,:)=work(t,:)+innovation
      end do
      x=work(p+burn+1:p+burn+n,:)
   end subroutine simulate_var

   subroutine var_impulse_response(phi,sigma,max_lag,response,orthogonal)
      real(dp), intent(in) :: phi(:,:,:),sigma(:,:)
      integer, intent(in) :: max_lag
      real(dp), allocatable, intent(out) :: response(:,:,:)
      logical, intent(in), optional :: orthogonal
      real(dp), allocatable :: psi(:,:,:),impact(:,:)
      logical :: orth
      integer :: h,istat,k
      orth=.true.
      if (present(orthogonal)) orth=orthogonal
      k=size(phi,1)
      call var_psi_weights(phi,max_lag,psi)
      if (orth) then
         call matrix_sqrt_symmetric(sigma,impact,istat)
         if (istat/=mts_success) impact=eye(k)
      else
         impact=eye(k)
      end if
      allocate(response(k,k,0:max_lag))
      do h=0,max_lag
         response(:,:,h)=matmul(psi(:,:,h),impact)
      end do
   end subroutine var_impulse_response

   subroutine forecast_error_variance_decomposition(phi,sigma,max_lag,decomposition,orthogonal)
      real(dp), intent(in) :: phi(:,:,:),sigma(:,:)
      integer, intent(in) :: max_lag
      real(dp), allocatable, intent(out) :: decomposition(:,:,:)
      logical, intent(in), optional :: orthogonal
      real(dp), allocatable :: irf(:,:,:)
      real(dp) :: denom
      integer :: h,i,j,s,k
      k=size(phi,1)
      call var_impulse_response(phi,sigma,max_lag,irf,orthogonal)
      allocate(decomposition(k,k,1:max_lag));decomposition=0.0_dp
      do h=1,max_lag
         do i=1,k
            denom=0.0_dp
            do s=0,h-1
               denom=denom+sum(irf(i,:,s)**2)
            end do
            do j=1,k
               do s=0,h-1
                  decomposition(i,j,h)=decomposition(i,j,h)+irf(i,j,s)**2
               end do
               if (denom>0.0_dp) decomposition(i,j,h)=decomposition(i,j,h)/denom
            end do
         end do
      end do
   end subroutine forecast_error_variance_decomposition

   subroutine granger_causality_test(x,p,causes,effects,result,include_mean)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: p,causes(:),effects(:)
      type(diagnostic_result), intent(out) :: result
      logical, intent(in), optional :: include_mean
      type(var_model) :: unrestricted,restricted
      logical, allocatable :: mask(:,:)
      real(dp) :: ssr_u,ssr_r
      integer :: k,d,offset,lag,c,e,df,ne
      logical :: mean_flag
      mean_flag=.true.;if(present(include_mean))mean_flag=include_mean
      call fit_var(x,p,unrestricted,mean_flag)
      k=size(x,2);d=k*p+merge(1,0,mean_flag);offset=merge(1,0,mean_flag)
      allocate(mask(d,k));mask=.true.
      do e=1,size(effects)
         if(effects(e)<1.or.effects(e)>k)cycle
         do lag=1,p
            do c=1,size(causes)
               if(causes(c)>=1.and.causes(c)<=k)mask(offset+(lag-1)*k+causes(c),effects(e))=.false.
            end do
         end do
      end do
      call fit_var(x,p,restricted,mean_flag,mask)
      ssr_u=sum(unrestricted%residuals(:,effects)**2)
      ssr_r=sum(restricted%residuals(:,effects)**2)
      df=p*size(causes)*size(effects)
      ne=size(unrestricted%residuals,1)
      result%degrees_freedom=df
      if(ssr_u>0.0_dp.and.df>0)then
         result%statistic=max(0.0_dp,real(ne,dp)*(ssr_r-ssr_u)/ssr_u)
         result%p_value=chi_square_survival(result%statistic,df)
         result%status=mts_success
      else
         result%status=mts_invalid_input
      end if
   end subroutine granger_causality_test

   pure logical function present_and_true(value)
      logical, intent(in), optional :: value
      present_and_true=.true.
      if(present(value))present_and_true=value
   end function present_and_true

end module mts_var
