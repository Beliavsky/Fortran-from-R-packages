module forecast_accuracy
   use forecast_kinds, only : dp
   use forecast_types, only : accuracy_result
   use forecast_stats, only : mean_value, autocorrelation, normal_cdf
   implicit none
   private
   public :: accuracy, dm_test
contains
   function accuracy(actual,predicted,training,m,nonseasonal_d,seasonal_d) result(out)
      real(dp),intent(in) :: actual(:),predicted(:)
      real(dp),intent(in),optional :: training(:)
      integer,intent(in),optional :: m,nonseasonal_d,seasonal_d
      type(accuracy_result) :: out
      real(dp),allocatable :: e(:),pe(:),z(:)
      real(dp)::scale
      integer::n,mm,dd,DDs,k
      n=min(size(actual),size(predicted))
      if(n==0)return
      e=actual(1:n)-predicted(1:n)
      allocate(pe(n))
      pe=0.0_dp
      where(abs(actual(1:n))>tiny(1.0_dp)) pe=100.0_dp*e/actual(1:n)
      out%me=mean_value(e)
      out%rmse=sqrt(sum(e*e)/real(n,dp))
      out%mae=sum(abs(e))/real(n,dp)
      out%mpe=mean_value(pe)
      out%mape=mean_value(abs(pe))
      if(n>1)out%acf1=autocorrelation(e,1)
      if(present(training)) then
         mm=1
         dd=0
         DDs=0
         if(present(m))mm=max(1,m)
         if(present(nonseasonal_d))dd=nonseasonal_d
         if(present(seasonal_d))DDs=seasonal_d
         z=training
         do k=1,DDs
            if(size(z)<=mm)exit
            z=z(1+mm:)-z(:size(z)-mm)
         end do
         do k=1,dd
            if(size(z)<=1)exit
            z=z(2:)-z(:size(z)-1)
         end do
         if(size(z)>0)then
         scale=mean_value(abs(z))
         else
         scale=0.0_dp
         end if
         if(scale>tiny(1.0_dp))out%mase=out%mae/scale
      end if
   end function

   subroutine dm_test(e1,e2,h,power,bartlett,statistic,p_value)
      real(dp),intent(in)::e1(:),e2(:)
      integer,intent(in),optional::h,power
      logical,intent(in),optional::bartlett
      real(dp),intent(out)::statistic,p_value
      real(dp),allocatable::d(:)
      real(dp)::mu,cov,dvar,w,kfac
      integer::hh,pw,n,lag
      logical::bw
      n=min(size(e1),size(e2))
      hh=1
      pw=2
      bw=.false.
      if(present(h))hh=max(1,min(h,n))
      if(present(power))pw=power
      if(present(bartlett))bw=bartlett
      d=abs(e1(1:n))**pw-abs(e2(1:n))**pw
      mu=mean_value(d)
      cov=sum((d-mu)**2)/real(n,dp)
      dvar=cov
      do lag=1,hh-1
         cov=sum((d(1:n-lag)-mu)*(d(1+lag:n)-mu))/real(n,dp)
         w=1.0_dp
         if(bw)w=1.0_dp-real(lag,dp)/real(hh,dp)
         dvar=dvar+2.0_dp*w*cov
      end do
      dvar=dvar/real(n,dp)
      if(dvar<=0.0_dp)then
      statistic=0.0_dp
      p_value=1.0_dp
      return
      end if
      statistic=mu/sqrt(dvar)
      kfac=sqrt(real(n+1-2*hh,dp)/real(n,dp)+real(hh*(hh-1),dp)/real(n*n,dp))
      statistic=statistic*kfac
      ! Normal approximation
      ! excellent for moderate samples and dependency-free.
      p_value=2.0_dp*(1.0_dp-normal_cdf(abs(statistic)))
   end subroutine
end module forecast_accuracy
