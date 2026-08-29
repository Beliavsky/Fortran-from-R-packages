module forecast_arfima
   use forecast_kinds, only : dp
   use forecast_types, only : forecast_result, arima_model
   use forecast_arima, only : arima_fit, auto_arima, arima_forecast
   use forecast_stats, only : normal_quantile
   use fracdiff, only : fracdiff_model, fracdiff_fit, diffseries, fractional_weights
   implicit none
   private
   public :: arfima_model, arfima_fit, arfima_forecast
   type :: arfima_model
      type(fracdiff_model) :: fractional
      real(dp) :: mean = 0.0_dp
      real(dp),allocatable :: x(:), fitted(:), residuals(:)
   end type
contains
   function arfima_fit(x,max_p,max_q) result(model)
      real(dp),intent(in)::x(:)
      integer,intent(in),optional::max_p,max_q
      type(arfima_model)::model
      type(fracdiff_model)::proxy,final
      type(arima_model)::arma
      real(dp),allocatable::centered(:),dx(:)
      integer::pmax,qmax,p,q,status
      model%mean=sum(x)/real(size(x),dp)
      centered=x-model%mean
      proxy=fracdiff_fit(centered,nar=min(2,max(0,size(x)/20)),nma=0)
      allocate(dx(size(x)))
      call diffseries(centered,proxy%d,dx,status)
      pmax=3
      if(present(max_p))pmax=max_p
      qmax=3
      if(present(max_q))qmax=max_q
      arma=auto_arima(dx,m=1,max_p=pmax,max_q=qmax,max_sp=0,max_sq=0,seasonal=.false.,stepwise=.true.)
      p=arma%p
      q=arma%q
      final=fracdiff_fit(centered,nar=p,nma=q)
      model%fractional=final
      model%x=x
      model%residuals=final%residuals
      model%fitted=x-model%residuals
   end function
   function arfima_forecast(model,h,levels) result(fc)
      type(arfima_model),intent(in)::model
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      type(arima_model)::arma
      real(dp),allocatable::centered(:),dx(:),fd(:),w(:),psi(:),lev(:)
      integer::status,n,k,j,p,q
      real(dp)::rhs,lhs,z
      centered=model%x-model%mean
      n=size(centered)
      allocate(dx(n))
      call diffseries(centered,model%fractional%d,dx,status)
      p=size(model%fractional%ar)
      q=size(model%fractional%ma)
      arma=arima_fit(dx,p,0,q,include_mean=.false.,optimize=.false.)
      if(p>0)arma%ar=model%fractional%ar
      if(q>0)arma%ma=-model%fractional%ma
      fc=arima_forecast(arma,dx,h,levels)
      fd=fc%mean
      if(allocated(fc%mean))deallocate(fc%mean)
      ! Fractional integration using generalized binomial weights.
      allocate(w(n+h+1))
      call fractional_weights(n+h+1,model%fractional%d,w,status)
      allocate(fc%mean(h))
      fc%mean=0.0_dp
      do k=1,h
         rhs=fd(k)
         do j=1,n
         rhs=rhs-w(j+1)*centered(n-j+1)
         end do
         lhs=0.0_dp
         if(k>1)then
            do j=1,k-1
            lhs=lhs+w(j+1)*fc%mean(k-j)
            end do
         end if
         fc%mean(k)=rhs-lhs
      end do
      fc%mean=fc%mean+model%mean
      if(present(levels))then
      lev=levels
      else
      lev=[80.0_dp,95.0_dp]
      end if
      if(allocated(fc%se))deallocate(fc%se)
      if(allocated(fc%level))deallocate(fc%level)
      if(allocated(fc%lower))deallocate(fc%lower)
      if(allocated(fc%upper))deallocate(fc%upper)
      allocate(fc%se(h),fc%level(size(lev)),fc%lower(h,size(lev)),fc%upper(h,size(lev)))
      fc%level=lev
      ! Fractional MA(infinity) approximation for forecast error variance.
      allocate(psi(h))
      psi=0.0_dp
      psi(1)=1.0_dp
      do k=2,h
      psi(k)=-w(k)
      end do
      do k=1,h
      fc%se(k)=sqrt(max(model%fractional%sigma**2*sum(psi(1:k)**2),0.0_dp))
      end do
      do j=1,size(lev)
      z=normal_quantile(0.5_dp+lev(j)/200.0_dp)
      fc%lower(:,j)=fc%mean-z*fc%se
      fc%upper(:,j)=fc%mean+z*fc%se
      end do
   end function
end module forecast_arfima
