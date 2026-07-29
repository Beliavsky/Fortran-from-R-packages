! SPDX-License-Identifier: GPL-3.0-only
module smoots_arma
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use smoots_kinds, only : dp
   use smoots_status, only : sm_ok, sm_invalid_input, sm_iteration_limit, sm_fit_failed
   use smoots_types, only : arma_model
   use smoots_linalg, only : solve_linear_system
   use smoots_stats, only : mean_value, random_normal
   implicit none
   private
   public :: fit_arma, arma_residuals, arma_point_forecast, simulate_arma
   public :: information_criterion_matrix, optimal_order, ma_infinity
   public :: estimate_cf0_ar, estimate_cf0_ma, estimate_cf0_arma
contains
   subroutine arma_residuals(x, ar, ma, mean_x, residuals, fitted)
      real(dp), intent(in) :: x(:), ar(:), ma(:), mean_x
      real(dp), intent(out) :: residuals(:), fitted(:)
      integer :: n,p,q,t,j
      real(dp) :: prediction
      n=size(x); p=size(ar); q=size(ma)
      residuals=0.0_dp; fitted=mean_x
      do t=1,n
         prediction=0.0_dp
         do j=1,min(p,t-1)
            prediction=prediction+ar(j)*(x(t-j)-mean_x)
         end do
         do j=1,min(q,t-1)
            prediction=prediction+ma(j)*residuals(t-j)
         end do
         fitted(t)=mean_x+prediction
         residuals(t)=x(t)-fitted(t)
      end do
   end subroutine arma_residuals

   subroutine residual_jacobian(x,ar,ma,mean_x,residuals,jacobian)
      real(dp),intent(in)::x(:),ar(:),ma(:),mean_x
      real(dp),intent(out)::residuals(:),jacobian(:,:)
      real(dp),allocatable::fitted(:)
      real(dp)::value
      integer::n,p,q,t,i,j
      n=size(x);p=size(ar);q=size(ma)
      allocate(fitted(n)); call arma_residuals(x,ar,ma,mean_x,residuals,fitted)
      jacobian=0.0_dp
      do i=1,p+q
         do t=1,n
            value=0.0_dp
            do j=1,min(q,t-1)
               value=value-ma(j)*jacobian(t-j,i)
            end do
            if (i<=p) then
               if (t>i) value=value-(x(t-i)-mean_x)
            else
               if (t>i-p) value=value-residuals(t-(i-p))
            end if
            jacobian(t,i)=value
         end do
      end do
   end subroutine residual_jacobian

   subroutine fit_arma(x,p,q,include_mean,model,status,tolerance,max_iterations)
      real(dp),intent(in)::x(:)
      integer,intent(in)::p,q
      logical,intent(in)::include_mean
      type(arma_model),intent(out)::model
      integer,intent(out)::status
      real(dp),intent(in),optional::tolerance
      integer,intent(in),optional::max_iterations
      real(dp),allocatable::res(:),trial_res(:),fit(:),jac(:,:),normal(:,:),damped(:,:),gradient(:),step(:)
      real(dp),allocatable::params(:),trial(:),trial_ar(:),trial_ma(:)
      real(dp)::tol,sse,trial_sse,lambda,relred,stepscale,mean_x
      integer::n,k,t0,iter,maxit,i,istat

      n=size(x); k=p+q; tol=1.0e-8_dp; if(present(tolerance))tol=tolerance
      maxit=150; if(present(max_iterations))maxit=max_iterations
      status=sm_ok
      model%p=p;model%q=q;model%include_mean=include_mean;model%status=sm_ok
      if(n<=max(p,q)+k+2 .or. p<0 .or. q<0)then
         status=sm_invalid_input;model%status=status;return
      end if
      allocate(model%ar(p),model%ma(q),model%residuals(n),model%fitted(n))
      model%ar=0.0_dp;model%ma=0.0_dp
      mean_x=merge(mean_value(x),0.0_dp,include_mean)
      model%mean_value=mean_x
      if(k==0)then
         call arma_residuals(x,model%ar,model%ma,mean_x,model%residuals,model%fitted)
         t0=1
         model%sigma2=sum(model%residuals(t0:n)**2)/real(n-t0+1,dp)
         call set_likelihood(model,n-t0+1)
         return
      end if
      allocate(res(n),trial_res(n),fit(n),jac(n,k),normal(k,k),damped(k,k),gradient(k),step(k))
      allocate(params(k),trial(k),trial_ar(p),trial_ma(q))
      params=0.0_dp
      call residual_jacobian(x,model%ar,model%ma,mean_x,res,jac)
      t0=max(p,q)+1
      sse=sum(res(t0:n)**2)
      normal=matmul(transpose(jac(t0:n,:)),jac(t0:n,:))
      gradient=matmul(transpose(jac(t0:n,:)),res(t0:n))
      lambda=1.0e-3_dp*max(1.0_dp,maxval([(abs(normal(i,i)),i=1,k)]))
      do iter=1,maxit
         if(maxval(abs(gradient))<=tol*(1.0_dp+sqrt(max(sse,0.0_dp))))exit
         damped=normal
         do i=1,k
            damped(i,i)=damped(i,i)+lambda*max(1.0_dp,normal(i,i))
         end do
         call solve_linear_system(damped,-gradient,step,istat)
         if(istat/=sm_ok)then
            lambda=lambda*10.0_dp
            if(lambda>1.0e18_dp)exit
            cycle
         end if
         trial=params+step
         if(maxval(abs(trial))>2.5_dp)then
            lambda=lambda*10.0_dp;cycle
         end if
         if(p>0)trial_ar=trial(1:p)
         if(q>0)trial_ma=trial(p+1:k)
         call arma_residuals(x,trial_ar,trial_ma,mean_x,trial_res,fit)
         trial_sse=sum(trial_res(t0:n)**2)
         if(ieee_is_finite(trial_sse).and.trial_sse<sse)then
            relred=(sse-trial_sse)/max(1.0_dp,sse)
            stepscale=sqrt(sum(step**2))/max(1.0_dp,sqrt(sum(params**2)))
            params=trial;sse=trial_sse
            if(p>0)model%ar=params(1:p)
            if(q>0)model%ma=params(p+1:k)
            call residual_jacobian(x,model%ar,model%ma,mean_x,res,jac)
            normal=matmul(transpose(jac(t0:n,:)),jac(t0:n,:))
            gradient=matmul(transpose(jac(t0:n,:)),res(t0:n))
            lambda=max(lambda/3.0_dp,1.0e-15_dp)
            if(relred<=tol.and.stepscale<=sqrt(tol))exit
         else
            lambda=lambda*10.0_dp
            if(lambda>1.0e18_dp)exit
         end if
      end do
      model%iterations=min(iter,maxit)
      call arma_residuals(x,model%ar,model%ma,mean_x,model%residuals,model%fitted)
      model%sigma2=sum(model%residuals(t0:n)**2)/real(n-t0+1,dp)
      if(.not.ieee_is_finite(model%sigma2).or.model%sigma2<=0.0_dp)then
         status=sm_fit_failed;model%status=status;return
      end if
      call set_likelihood(model,n-t0+1)
      if(iter>maxit)then
         status=sm_iteration_limit;model%status=status
      end if
   end subroutine fit_arma

   subroutine set_likelihood(model,neff)
      type(arma_model),intent(inout)::model
      integer,intent(in)::neff
      integer::k
      model%log_likelihood=-0.5_dp*real(neff,dp)*(log(2.0_dp*acos(-1.0_dp)*model%sigma2)+1.0_dp)
      k=model%p+model%q+merge(1,0,model%include_mean)+1
      model%aic=-2.0_dp*model%log_likelihood+2.0_dp*real(k,dp)
      model%bic=-2.0_dp*model%log_likelihood+log(real(neff,dp))*real(k,dp)
   end subroutine set_likelihood

   subroutine arma_point_forecast(x,residuals,ar,ma,mean_x,h,forecast,status)
      real(dp),intent(in)::x(:),residuals(:),ar(:),ma(:),mean_x
      integer,intent(in)::h
      real(dp),allocatable,intent(out)::forecast(:)
      integer,intent(out)::status
      real(dp),allocatable::xx(:),ee(:)
      integer::n,p,q,t,j
      n=size(x);p=size(ar);q=size(ma);allocate(forecast(h));forecast=0.0_dp
      if(h<1.or.size(residuals)/=n)then;status=sm_invalid_input;return;end if
      allocate(xx(n+h),ee(n+h));xx(1:n)=x-mean_x;ee=0.0_dp;ee(1:n)=residuals
      do t=n+1,n+h
         do j=1,p
            xx(t)=xx(t)+ar(j)*xx(t-j)
         end do
         do j=1,q
            xx(t)=xx(t)+ma(j)*ee(t-j)
         end do
      end do
      forecast=xx(n+1:n+h)+mean_x;status=sm_ok
   end subroutine arma_point_forecast

   subroutine simulate_arma(ar,ma,mean_x,n,burn,series,status,innovations,start_innovations)
      real(dp),intent(in)::ar(:),ma(:),mean_x
      integer,intent(in)::n,burn
      real(dp),allocatable,intent(out)::series(:)
      integer,intent(out)::status
      real(dp),intent(in),optional::innovations(:),start_innovations(:)
      real(dp),allocatable::z(:),e(:)
      integer::total,t,j,p,q
      p=size(ar);q=size(ma);total=n+burn;allocate(series(n))
      if(n<1.or.burn<0)then;status=sm_invalid_input;return;end if
      allocate(z(total),e(total));z=0.0_dp;e=0.0_dp
      if(present(start_innovations))then
         e(1:min(burn,size(start_innovations)))=start_innovations(1:min(burn,size(start_innovations)))
      else
         do t=1,burn;e(t)=random_normal();end do
      end if
      if(present(innovations))then
         if(size(innovations)<n)then;status=sm_invalid_input;return;end if
         e(burn+1:total)=innovations(1:n)
      else
         do t=burn+1,total;e(t)=random_normal();end do
      end if
      do t=1,total
         z(t)=e(t)
         do j=1,min(p,t-1);z(t)=z(t)+ar(j)*z(t-j);end do
         do j=1,min(q,t-1);z(t)=z(t)+ma(j)*e(t-j);end do
      end do
      series=z(burn+1:total)+mean_x;status=sm_ok
   end subroutine simulate_arma

   subroutine ma_infinity(ar,ma,m,coefficients,status)
      real(dp),intent(in)::ar(:),ma(:)
      integer,intent(in)::m
      real(dp),allocatable,intent(out)::coefficients(:)
      integer,intent(out)::status
      real(dp),allocatable::ma_pad(:),work(:)
      integer::p,q,i,j
      p=size(ar);q=size(ma);allocate(coefficients(m+1));coefficients=0.0_dp
      if(m<0)then;status=sm_invalid_input;return;end if
      allocate(ma_pad(max(m,q)));ma_pad=0.0_dp
      if(q>0)ma_pad(1:q)=ma
      allocate(work(max(1,p-1)+1+m));work=0.0_dp
      if(p==0)then
         coefficients(1)=1.0_dp
         if(m>0)coefficients(2:m+1)=ma_pad(1:m)
         status=sm_ok;return
      end if
      work(p)=1.0_dp
      do i=p+1,p+m
         do j=1,p
            work(i)=work(i)+ar(j)*work(i-j)
         end do
         work(i)=work(i)+ma_pad(i-p)
      end do
      coefficients=work(p:p+m);status=sm_ok
   end subroutine ma_infinity

   subroutine information_criterion_matrix(x,pmax,qmax,include_mean,use_bic,matrix,status)
      real(dp),intent(in)::x(:)
      integer,intent(in)::pmax,qmax
      logical,intent(in)::include_mean,use_bic
      real(dp),allocatable,intent(out)::matrix(:,:)
      integer,intent(out)::status
      type(arma_model)::model
      integer::p,q,istat
      allocate(matrix(pmax+1,qmax+1));matrix=huge(1.0_dp);status=sm_ok
      if(pmax<0.or.qmax<0)then;status=sm_invalid_input;return;end if
      do p=0,pmax
         do q=0,qmax
            call fit_arma(x,p,q,include_mean,model,istat,max_iterations=100)
            if(istat==sm_ok.or.istat==sm_iteration_limit)then
               if(use_bic)then
                  matrix(p+1,q+1)=-2.0_dp*model%log_likelihood+log(real(size(x),dp))*real(p+q,dp)
               else
                  matrix(p+1,q+1)=model%aic
               end if
            end if
         end do
      end do
   end subroutine information_criterion_matrix

   subroutine optimal_order(matrix,p,q,minimize,mask,status)
      real(dp),intent(in)::matrix(:,:)
      integer,intent(out)::p,q,status
      logical,intent(in),optional::minimize
      logical,intent(in),optional::mask(:,:)
      logical::do_min,allowed
      real(dp)::best
      integer::i,j
      do_min=.true.;if(present(minimize))do_min=minimize
      best=merge(huge(1.0_dp),-huge(1.0_dp),do_min);p=-1;q=-1
      do j=1,size(matrix,2);do i=1,size(matrix,1)
         allowed=.true.;if(present(mask))allowed=mask(i,j)
         if(allowed)then
            if((do_min.and.matrix(i,j)<best).or.((.not.do_min).and.matrix(i,j)>best))then
               best=matrix(i,j);p=i-1;q=j-1
            end if
         end if
      end do;end do
      status=merge(sm_ok,sm_invalid_input,p>=0)
   end subroutine optimal_order

   subroutine estimate_cf0_ar(x,cf0,model,status,pmax)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::cf0
      type(arma_model),intent(out)::model
      integer,intent(out)::status
      integer,intent(in),optional::pmax
      real(dp),allocatable::mat(:,:)
      integer::pm,p,q,istat
      pm=5;if(present(pmax))pm=pmax
      call information_criterion_matrix(x,pm,0,.true.,.true.,mat,istat)
      call optimal_order(mat,p,q,status=status)
      call fit_arma(x,p,0,.true.,model,status)
      cf0=model%sigma2/(1.0_dp-sum(model%ar))**2
   end subroutine estimate_cf0_ar

   subroutine estimate_cf0_ma(x,cf0,model,status,qmax)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::cf0
      type(arma_model),intent(out)::model
      integer,intent(out)::status
      integer,intent(in),optional::qmax
      real(dp),allocatable::mat(:,:)
      integer::qm,p,q,istat
      qm=5;if(present(qmax))qm=qmax
      call information_criterion_matrix(x,0,qm,.true.,.true.,mat,istat)
      call optimal_order(mat,p,q,status=status)
      call fit_arma(x,0,q,.true.,model,status)
      cf0=model%sigma2*(1.0_dp+sum(model%ma))**2
   end subroutine estimate_cf0_ma

   subroutine estimate_cf0_arma(x,cf0,model,status,pmax,qmax)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::cf0
      type(arma_model),intent(out)::model
      integer,intent(out)::status
      integer,intent(in),optional::pmax,qmax
      real(dp),allocatable::mat(:,:)
      integer::pm,qm,p,q,istat
      pm=5;qm=5;if(present(pmax))pm=pmax;if(present(qmax))qm=qmax
      call information_criterion_matrix(x,pm,qm,.true.,.true.,mat,istat)
      call optimal_order(mat,p,q,status=status)
      call fit_arma(x,p,q,.true.,model,status)
      cf0=model%sigma2*((1.0_dp+sum(model%ma))/(1.0_dp-sum(model%ar)))**2
   end subroutine estimate_cf0_arma
end module smoots_arma
