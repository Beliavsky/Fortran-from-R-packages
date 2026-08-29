module forecast_diagnostics
   use forecast_kinds, only : dp
   use forecast_types, only : regression_model, arima_model, ets_model, bats_model
   use forecast_stats, only : mean_value
   use forecast_linalg, only : inverse_matrix
   implicit none
   private
   public :: cross_correlation, ljung_box, regression_cv_stats
   public :: model_df_arima, model_df_ets, model_df_bats
contains
   function cross_correlation(x,y,lag_max,covariance) result(ccf)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in)::lag_max
      logical,intent(in),optional::covariance
      real(dp),allocatable::ccf(:)
      real(dp)::mx,my,sx,sy,val
      integer::lag,i,n,count
      logical::covar
      if(size(x)/=size(y))error stop 'cross_correlation: size mismatch'
      n=size(x)
      mx=mean_value(x)
      my=mean_value(y)
      sx=sqrt(sum((x-mx)**2)/real(n,dp))
      sy=sqrt(sum((y-my)**2)/real(n,dp))
      covar=.false.
      if(present(covariance))covar=covariance
      allocate(ccf(2*lag_max+1))
      do lag=-lag_max,lag_max
         val=0.0_dp
         count=0
         do i=1,n
            if(i+lag>=1 .and. i+lag<=n)then
            val=val+(x(i)-mx)*(y(i+lag)-my)
            count=count+1
            end if
         end do
         if(covar)then
         ccf(lag+lag_max+1)=val/real(n,dp)
         else
         ccf(lag+lag_max+1)=val/(real(n,dp)*max(sx*sy,1.0e-300_dp))
         end if
      end do
   end function cross_correlation
   subroutine ljung_box(residuals,lag,df,statistic)
      real(dp),intent(in)::residuals(:)
      integer,intent(in)::lag
      integer,intent(in),optional::df
      real(dp),intent(out)::statistic
      real(dp)::mu,den,rho
      integer::n,k,dof
      n=size(residuals)
      mu=mean_value(residuals)
      den=sum((residuals-mu)**2)
      statistic=0.0_dp
      do k=1,min(lag,n-1)
         rho=sum((residuals(1:n-k)-mu)*(residuals(1+k:n)-mu))/max(den,1.0e-300_dp)
         statistic=statistic+rho*rho/real(n-k,dp)
      end do
      statistic=real(n*(n+2),dp)*statistic
      dof=0
      if(present(df))dof=df
      if(lag<=dof)statistic=0.0_dp
   end subroutine ljung_box
   function regression_cv_stats(model,x,y) result(stats)
      type(regression_model),intent(in)::model
      real(dp),intent(in)::x(:,:),y(:)
      real(dp)::stats(5)
      real(dp),allocatable::design(:,:),xtxinv(:,:),h(:)
      real(dp)::rss,tss,aic,aicc,bic,adjr2,cv
      integer::n,k,i,info,p
      n=size(y)
      p=size(model%coefficients)
      if(size(x,2)==p)then
      design=x
      else if(size(x,2)==p-1)then
      allocate(design(n,p))
      design(:,1)=1.0_dp
      design(:,2:)=x
      else
      error stop 'regression_cv_stats: dimensions'
      end if
      call inverse_matrix(matmul(transpose(design),design),xtxinv,info)
      if(info/=0)error stop 'regression_cv_stats: singular design'
      allocate(h(n))
      do i=1,n
      h(i)=dot_product(design(i,:),matmul(xtxinv,design(i,:)))
      end do
      rss=sum(model%residuals**2)
      k=p-1
      aic=real(n,dp)*log(max(rss/real(n,dp),1.0e-300_dp))+2.0_dp*real(k+2,dp)
      if(n>k+3)then
      aicc=aic+2.0_dp*real((k+2)*(k+3),dp)/real(n-k-3,dp)
      else
      aicc=huge(1.0_dp)
      end if
      bic=aic+real(k+2,dp)*(log(real(n,dp))-2.0_dp)
      cv=sum((model%residuals/max(1.0_dp-h,1.0e-12_dp))**2)/real(n,dp)
      tss=sum((y-mean_value(y))**2)
      adjr2=1.0_dp-(rss/real(max(1,n-p),dp))/(max(tss,1.0e-300_dp)/real(max(1,n-1),dp))
      stats=[cv,aic,aicc,bic,adjr2]
   end function regression_cv_stats
   integer function model_df_arima(model) result(df)
      type(arima_model),intent(in)::model
      df=size(model%ar)+size(model%ma)+size(model%sar)+size(model%sma)+merge(1,0,model%include_mean)+merge(1,0, &
         & model%include_drift)
      if(allocated(model%xreg_coef))df=df+size(model%xreg_coef)
   end function model_df_arima
   integer function model_df_ets(model) result(df)
      type(ets_model),intent(in)::model
      df=1+merge(1,0,model%trend_type/=0)+merge(1,0,model%season_type/=0)+merge(1,0,model%damped)
   end function model_df_ets
   integer function model_df_bats(model) result(df)
      type(bats_model),intent(in)::model
      df=1+merge(1,0,model%has_trend)+merge(1,0, &
         & model%phi<0.999999_dp)+size(model%gamma1)+size(model%gamma2)+model%p+model%q+merge(1,0,model%use_boxcox)
   end function model_df_bats
end module forecast_diagnostics
