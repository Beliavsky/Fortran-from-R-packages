module forecast_cv
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use forecast_kinds, only : dp
   use forecast_types, only : forecast_result, cvar_result, accuracy_result
   use forecast_accuracy, only : accuracy
   use forecast_diagnostics, only : ljung_box
   use urca_distributions, only : chi_square_sf
   implicit none
   private
   public :: ts_cv, cvar
   abstract interface
      function forecast_callback(y,h) result(fc)
         import dp,forecast_result
         real(dp),intent(in)::y(:)
         integer,intent(in)::h
         type(forecast_result)::fc
      end function forecast_callback
      function subset_fit_callback(y,train_mask) result(fitted)
         import dp
         real(dp),intent(in)::y(:)
         logical,intent(in)::train_mask(:)
         real(dp),allocatable::fitted(:)
      end function subset_fit_callback
   end interface
contains
   function ts_cv(y,fun,h,window,initial) result(errors)
      real(dp),intent(in)::y(:)
      procedure(forecast_callback)::fun
      integer,intent(in),optional::h,window,initial
      real(dp),allocatable::errors(:,:)
      type(forecast_result)::fc
      integer::hh,w,init,i,start,j,n
      n=size(y)
      hh=1
      if(present(h))hh=max(1,h)
      w=0
      if(present(window))w=max(1,window)
      init=0
      if(present(initial))init=max(0,initial)
      if(init>=n)error stop 'ts_cv: initial period too long'
      allocate(errors(n,hh))
      errors=ieee_value(0.0_dp,ieee_quiet_nan)
      do i=1+init,n-1
         if(w>0 .and. i<w+init)cycle
         start=1
         if(w>0)start=i-w+1
         fc=fun(y(start:i),hh)
         do j=1,min(hh,n-i)
            if(j<=size(fc%mean))errors(i,j)=y(i+j)-fc%mean(j)
         end do
      end do
   end function ts_cv

   function cvar(y,fun,k,blocked,lb_lags,fold_ids) result(out)
      real(dp),intent(in)::y(:)
      procedure(subset_fit_callback)::fun
      integer,intent(in),optional::k,lb_lags,fold_ids(:)
      logical,intent(in),optional::blocked
      type(cvar_result)::out
      integer,allocatable::fold(:),idx(:)
      logical,allocatable::train(:),test(:)
      real(dp),allocatable::fit(:)
      type(accuracy_result)::acc
      integer::n,kk,lag,i,j,nj,tmp
      real(dp)::u
      logical::blk
      n=size(y)
      kk=min(10,n)
      if(present(k))kk=min(max(2,k),n)
      if(n<2)error stop 'cvar: at least two observations are required'
      blk=.false.
      if(present(blocked))blk=blocked
      lag=min(24,max(1,n-1))
      if(present(lb_lags))lag=min(max(1,lb_lags),max(1,n-1))
      allocate(fold(n),idx(n),train(n),test(n),out%testfit(n),out%residuals(n))
      allocate(out%fold_accuracy(kk,7),out%cv_mean(7),out%cv_sd(7))
      out%testfit=ieee_value(0.0_dp,ieee_quiet_nan)
      if(present(fold_ids))then
         if(size(fold_ids)/=n)error stop 'cvar: fold_ids size mismatch'
         if(minval(fold_ids)<1 .or. maxval(fold_ids)>kk)error stop 'cvar: invalid fold_ids'
         fold=fold_ids
      else if(blk)then
         do i=1,n
            fold(i)=1+((i-1)*kk)/n
         end do
      else
         do i=1,n
            idx(i)=i
            fold(i)=1+mod(i-1,kk)
         end do
         do i=n,2,-1
            call random_number(u)
            j=1+int(u*real(i,dp))
            j=min(j,i)
            tmp=fold(i)
            fold(i)=fold(j)
            fold(j)=tmp
         end do
      end if
      do i=1,kk
         test=(fold==i)
         train=.not.test
         nj=count(test)
         if(nj==0)error stop 'cvar: empty fold'
         fit=fun(y,train)
         if(size(fit)/=n)error stop 'cvar: callback must return a fitted vector of length size(y)'
         out%testfit=merge(fit,out%testfit,test)
         acc=accuracy(pack(y,test),pack(fit,test),training=y)
         out%fold_accuracy(i,:)=[acc%me,acc%rmse,acc%mae,acc%mpe,acc%mape,acc%mase,acc%acf1]
      end do
      out%residuals=out%testfit-y
      do j=1,7
         out%cv_mean(j)=sum(out%fold_accuracy(:,j))/real(kk,dp)
         if(kk>1)then
            out%cv_sd(j)=sqrt(sum((out%fold_accuracy(:,j)-out%cv_mean(j))**2)/real(kk-1,dp))
         else
            out%cv_sd(j)=0.0_dp
         end if
      end do
      call ljung_box(out%residuals,lag,statistic=out%lb_statistic)
      out%lb_p_value=chi_square_sf(max(0.0_dp,out%lb_statistic),real(lag,dp))
      out%lb_lag=lag
      out%k=kk
   end function cvar
end module forecast_cv
