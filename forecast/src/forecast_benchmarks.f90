module forecast_benchmarks
   use forecast_kinds, only : dp
   use forecast_types, only : forecast_result,croston_fit,theta_fit
   use forecast_stats, only : mean_value,variance_value,normal_quantile
   implicit none
   private
   public :: mean_forecast, random_walk_forecast, naive_forecast, seasonal_naive_forecast
   public :: fit_croston, forecast_croston, fit_theta, forecast_theta
contains
   function mean_forecast(y,h,levels) result(fc)
      real(dp),intent(in)::y(:)
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      real(dp),allocatable::lev(:)
      real(dp)::mu,sig,se,z
      integer::j
      if(present(levels))then
      lev=levels
      else
      lev=[80.0_dp,95.0_dp]
      end if
      allocate(fc%mean(h),fc%se(h),fc%level(size(lev)),fc%lower(h,size(lev)),fc%upper(h,size(lev)))
      mu=mean_value(y)
      sig=sqrt(variance_value(y,.true.))
      se=sig*sqrt(1.0_dp+1.0_dp/real(max(1,size(y)),dp))
      fc%mean=mu
      fc%se=se
      fc%level=lev
      do j=1,size(lev)
      z=normal_quantile(0.5_dp+lev(j)/200.0_dp)
      fc%lower(:,j)=mu-z*se
      fc%upper(:,j)=mu+z*se
      end do
   end function

   function random_walk_forecast(y,h,lag,drift,levels) result(fc)
      real(dp),intent(in)::y(:)
      integer,intent(in)::h
      integer,intent(in),optional::lag
      logical,intent(in),optional::drift
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      integer::l,i,j,n,step
      logical::dr
      real(dp)::b,s2,bse,z
      real(dp),allocatable::e(:),lev(:)
      l=1
      if(present(lag))l=max(1,lag)
      dr=.false.
      if(present(drift))dr=drift
      n=size(y)
      if(present(levels))then
      lev=levels
      else
      lev=[80.0_dp,95.0_dp]
      end if
      b=0.0_dp
      bse=0.0_dp
      if(dr .and. n>l)then
      b=(y(n)-y(1))/real(max(1,n-1),dp)
      end if
      if(n>l)then
      e=y(1+l:n)-y(1:n-l)-b*real(l,dp)
      s2=sum(e*e)/real(max(1,size(e)-merge(1,0,dr)),dp)
      else
      s2=0.0_dp
      end if
      if(dr .and. n>2) bse=sqrt(s2/real(n-1,dp))
      allocate(fc%mean(h),fc%se(h),fc%level(size(lev)),fc%lower(h,size(lev)),fc%upper(h,size(lev)))
      fc%level=lev
      do i=1,h
         step=(i-1)/l+1
         fc%mean(i)=y(n-mod(h-i,l))+real(step,dp)*b
         ! Match rwf structure: repeat last lag observations cyclically.
         fc%mean(i)=y(n-l+mod(i-1,l)+1)+real(step,dp)*b
         fc%se(i)=sqrt(max(0.0_dp,s2*real(step,dp)+(real(step,dp)*bse)**2))
      end do
      do j=1,size(lev)
      z=normal_quantile(0.5_dp+lev(j)/200.0_dp)
      fc%lower(:,j)=fc%mean-z*fc%se
      fc%upper(:,j)=fc%mean+z*fc%se
      end do
   end function
   function naive_forecast(y,h,levels) result(fc)
      real(dp),intent(in)::y(:)
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      fc=random_walk_forecast(y,h,1,.false.,levels)
   end function
   function seasonal_naive_forecast(y,h,m,levels) result(fc)
      real(dp),intent(in)::y(:)
      integer,intent(in)::h,m
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      fc=random_walk_forecast(y,h,m,.false.,levels)
   end function

   function fit_croston(y,alpha,variant) result(model)
      real(dp),intent(in)::y(:)
      real(dp),intent(in),optional::alpha
      integer,intent(in),optional::variant
      type(croston_fit)::model
      real(dp)::a,fd,fi,ratio,coeff
      integer::i,n,lastnz,interval,v
      a=0.1_dp
      if(present(alpha))a=alpha
      v=1
      if(present(variant))v=variant
      model%alpha=a
      model%variant=v
      n=size(y)
      allocate(model%fitted(n),model%residuals(n))
      model%fitted=0.0_dp
      lastnz=0
      fd=0.0_dp
      fi=1.0_dp
      ratio=0.0_dp
      select case(v)
      case(2)
      coeff=1.0_dp-a/2.0_dp
      case(3)
      coeff=1.0_dp-a/(2.0_dp-a)
      case default
      coeff=1.0_dp
      end select
      do i=1,n
         if(i>1)model%fitted(i)=ratio
         if(y(i)/=0.0_dp)then
            if(lastnz==0)then
            fd=y(i)
            fi=real(i,dp)
            else
            interval=i-lastnz
            fd=fd+a*(y(i)-fd)
            fi=fi+a*(real(interval,dp)-fi)
            end if
            lastnz=i
            ratio=coeff*fd/max(fi,tiny(1.0_dp))
         end if
      end do
      model%level=fd
      model%interval=fi
      model%forecast=ratio
      model%residuals=y-model%fitted
   end function
   function forecast_croston(model,h) result(fc)
      type(croston_fit),intent(in)::model
      integer,intent(in)::h
      type(forecast_result)::fc
      allocate(fc%mean(h))
      fc%mean=model%forecast
   end function

   function fit_theta(y,alpha) result(model)
      real(dp),intent(in)::y(:)
      real(dp),intent(in),optional::alpha
      type(theta_fit)::model
      real(dp)::a,l,b
      integer::i,n
      n=size(y)
      a=0.2_dp
      if(present(alpha))a=alpha
      ! Theta equals SES with drift
      ! estimate linear slope and choose SES alpha by one-dimensional grid if not supplied.
      if(n>1)then
         b=12.0_dp*sum(([(real(i,dp),i=1,n)]-0.5_dp*real(n+1,dp))*y)/real(n*(n*n-1),dp)
      else
      b=0.0_dp
      end if
      model%drift=0.5_dp*b
      model%alpha=a
      allocate(model%fitted(n),model%residuals(n))
      l=y(1)
      model%fitted(1)=y(1)
      do i=2,n
      model%fitted(i)=l+model%drift
      l=a*y(i)+(1.0_dp-a)*(l+model%drift)
      end do
      model%level=l
      model%residuals=y-model%fitted
      model%sigma2=sum(model%residuals(2:)**2)/real(max(1,n-1),dp)
   end function
   function forecast_theta(model,h,levels) result(fc)
      type(theta_fit),intent(in)::model
      integer,intent(in)::h
      real(dp),intent(in),optional::levels(:)
      type(forecast_result)::fc
      real(dp),allocatable::lev(:)
      integer::i,j
      real(dp)::z
      if(present(levels))then
      lev=levels
      else
      lev=[80.0_dp,95.0_dp]
      end if
      allocate(fc%mean(h),fc%se(h),fc%level(size(lev)),fc%lower(h,size(lev)),fc%upper(h,size(lev)))
      fc%level=lev
      do i=1,h
      fc%mean(i)=model%level+real(i,dp)*model%drift
      fc%se(i)=sqrt(model%sigma2*real(i,dp))
      end do
      do j=1,size(lev)
      z=normal_quantile(0.5_dp+lev(j)/200.0_dp)
      fc%lower(:,j)=fc%mean-z*fc%se
      fc%upper(:,j)=fc%mean+z*fc%se
      end do
   end function
end module forecast_benchmarks
