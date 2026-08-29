module forecast_modelar
   use forecast_kinds, only : dp
   use forecast_types, only : regression_model, forecast_result
   use forecast_regression, only : tslm_fit
   implicit none
   private
   public :: linear_ar_model, modelar_fit, modelar_forecast

   type :: linear_ar_model
      integer :: p = 0
      integer :: m = 1
      integer :: seasonal_p = 0
      integer, allocatable :: lags(:)
      type(regression_model) :: regression
      real(dp), allocatable :: fitted(:), residuals(:)
   end type linear_ar_model
contains
   function build_lags(p,m,seasonal_p) result(lags)
      integer,intent(in)::p,m,seasonal_p
      integer,allocatable::lags(:),tmp(:)
      integer::i,j,n,tmpval
      allocate(tmp(max(0,p)+max(0,seasonal_p)))
      n=0
      do i=1,p
         n=n+1
         tmp(n)=i
      end do
      do i=1,seasonal_p
         if (m*i <= 0) cycle
         if (.not.any(tmp(1:n)==m*i)) then
            n=n+1
            tmp(n)=m*i
         end if
      end do
      allocate(lags(n))
      if(n>0)lags=tmp(1:n)
      do i=1,n-1
         do j=i+1,n
            if(lags(j)<lags(i))then
               tmpval=lags(i)
               lags(i)=lags(j)
               lags(j)=tmpval
            end if
         end do
      end do
   end function build_lags

   function modelar_fit(y,p,m,seasonal_p) result(model)
      real(dp),intent(in)::y(:)
      integer,intent(in),optional::p,m,seasonal_p
      type(linear_ar_model)::model
      real(dp),allocatable::x(:,:),yt(:)
      integer::pp,mm,sp,maxlag,i,j,nobs
      real(dp)::best_aic,aic
      type(regression_model)::trial

      mm=1
      if(present(m))mm=max(1,m)
      sp=0
      if(present(seasonal_p))sp=max(0,seasonal_p)
      if(present(p))then
         pp=max(0,p)
      else
         best_aic=huge(1.0_dp)
         pp=0
         do i=0,min(10,max(0,size(y)/5))
            model%lags=build_lags(i,mm,sp)
            maxlag=0
            if(size(model%lags)>0)maxlag=maxval(model%lags)
            if(size(y)-maxlag<=max(3,size(model%lags)+1))cycle
            nobs=size(y)-maxlag
            allocate(x(nobs,size(model%lags)),yt(nobs))
            yt=y(maxlag+1:)
            do j=1,size(model%lags)
            x(:,j)=y(maxlag+1-model%lags(j):size(y)-model%lags(j))
            end do
            trial=tslm_fit(yt,x,.true.)
            aic=real(nobs,dp)*log(max(sum(trial%residuals**2)/real(nobs,dp),1.0e-300_dp))+2.0_dp*real(size(trial%coefficients), &
               & dp)
            if(aic<best_aic)then
            best_aic=aic
            pp=i
            end if
            deallocate(x,yt)
         end do
      end if
      model%p=pp
      model%m=mm
      model%seasonal_p=sp
      model%lags=build_lags(pp,mm,sp)
      maxlag=0
      if(size(model%lags)>0)maxlag=maxval(model%lags)
      nobs=size(y)-maxlag
      if(nobs<=size(model%lags)+1)error stop 'modelar_fit: insufficient observations'
      allocate(x(nobs,size(model%lags)),yt(nobs))
      yt=y(maxlag+1:)
      do j=1,size(model%lags)
      x(:,j)=y(maxlag+1-model%lags(j):size(y)-model%lags(j))
      end do
      model%regression=tslm_fit(yt,x,.true.)
      allocate(model%fitted(size(y)),model%residuals(size(y)))
      model%fitted=0.0_dp
      model%residuals=0.0_dp
      if(maxlag>0)model%fitted(1:maxlag)=y(1:maxlag)
      model%fitted(maxlag+1:)=model%regression%fitted
      model%residuals(maxlag+1:)=model%regression%residuals
   end function modelar_fit

   function modelar_forecast(model,y,h) result(fc)
      type(linear_ar_model),intent(in)::model
      real(dp),intent(in)::y(:)
      integer,intent(in)::h
      type(forecast_result)::fc
      real(dp),allocatable::history(:),row(:)
      integer::i,j,n0
      n0=size(y)
      allocate(history(n0+h))
      history(1:n0)=y
      allocate(row(size(model%lags)))
      allocate(fc%mean(h),fc%se(h))
      fc%se=sqrt(max(model%regression%sigma2,0.0_dp))
      do i=1,h
         do j=1,size(model%lags)
         row(j)=history(n0+i-model%lags(j))
         end do
         fc%mean(i)=model%regression%coefficients(1)+dot_product(model%regression%coefficients(2:),row)
         history(n0+i)=fc%mean(i)
      end do
   end function modelar_forecast
end module forecast_modelar
