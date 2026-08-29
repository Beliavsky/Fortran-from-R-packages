module forecast_unitroot
   use forecast_kinds, only : dp, pi
   use urca, only : ur_test_result, kpss_test, adf_test, pp_test, KPSS_MU, KPSS_TAU, UR_DRIFT, UR_TREND, &
      PP_ZTAU, PP_CONSTANT, PP_TREND, LAG_FIXED
   use forecast_features, only : difference_series, seasonal_difference
   use forecast_stats, only : autocorrelation, variance_value
   use forecast_linalg, only : least_squares, inverse_matrix
   use forecast_decompose, only : seasonal_strength
   implicit none
   private
   public :: ndiffs, nsdiffs, ocsb_test, ocsb_critical_value
contains
   pure real(dp) function critical_at_alpha(values,probs,alpha) result(cv)
      real(dp),intent(in)::values(:),probs(:),alpha
      integer::i
      real(dp)::t
      if(size(values)/=size(probs))error stop 'critical_at_alpha: size mismatch'
      if(alpha<=minval(probs))then
         cv=values(minloc(probs,dim=1))
         return
      else if(alpha>=maxval(probs))then
         cv=values(maxloc(probs,dim=1))
         return
      end if
      do i=1,size(probs)-1
         if((alpha-probs(i))*(alpha-probs(i+1))<=0.0_dp)then
            t=(alpha-probs(i))/(probs(i+1)-probs(i))
            cv=values(i)+t*(values(i+1)-values(i))
            return
         end if
      end do
      cv=values(1)
   end function critical_at_alpha

   integer function ndiffs(y,max_d,test,type,alpha,max_lags) result(d)
      real(dp),intent(in)::y(:)
      integer,intent(in),optional::max_d,max_lags
      character(len=*),intent(in),optional::test,type
      real(dp),intent(in),optional::alpha
      real(dp),allocatable::z(:)
      type(ur_test_result)::ut
      integer::md,ml,model
      real(dp)::a,cv
      character(len=8)::which,det
      logical::dodiff
      md=2
      if(present(max_d))md=max(0,max_d)
      ml=1
      if(present(max_lags))ml=max(0,max_lags)
      which='kpss'
      if(present(test))which=adjustl(test)
      det='level'
      if(present(type))det=adjustl(type)
      a=0.05_dp
      if(present(alpha))a=max(0.01_dp,min(0.10_dp,alpha))
      z=y
      d=0
      do while(d<md .and. size(z)>=10)
         select case(trim(which))
         case('adf','ADF')
            model=merge(UR_TREND,UR_DRIFT,trim(det)=='trend' .or. trim(det)=='TREND')
            ut=adf_test(z,model,ml,LAG_FIXED)
            if(ut%info/=0)exit
            cv=critical_at_alpha(ut%critical_values(1,:),[0.01_dp,0.05_dp,0.10_dp],a)
            dodiff=ut%statistic(1)>cv
         case('pp','PP')
            model=merge(PP_TREND,PP_CONSTANT,trim(det)=='trend' .or. trim(det)=='TREND')
            ut=pp_test(z,PP_ZTAU,model)
            if(ut%info/=0)exit
            cv=critical_at_alpha(ut%critical_values(1,:),[0.01_dp,0.05_dp,0.10_dp],a)
            dodiff=ut%statistic(1)>cv
         case default
            model=merge(KPSS_TAU,KPSS_MU,trim(det)=='trend' .or. trim(det)=='TREND')
            ut=kpss_test(z,model,use_lag=int(3.0_dp*sqrt(real(size(z),dp))/13.0_dp))
            if(ut%info/=0)exit
            cv=critical_at_alpha(ut%critical_values(1,:),[0.10_dp,0.05_dp,0.025_dp,0.01_dp],a)
            dodiff=ut%statistic(1)>cv
         end select
         if(.not.dodiff)exit
         z=difference_series(z,1)
         d=d+1
      end do
   end function ndiffs

   pure real(dp) function ocsb_critical_value(period) result(cv)
      integer,intent(in)::period
      real(dp)::logm
      if(period<=1)then
         cv=huge(1.0_dp)
         return
      end if
      logm=log(real(period,dp))
      cv=-0.2937411_dp*exp(-0.2850853_dp*(logm-0.7656451_dp) &
         -0.05983644_dp*(logm-0.7656451_dp)**2)-1.652202_dp
   end function ocsb_critical_value

   subroutine lag_regression_series(y,lag,maxlag,response,xreg)
      ! Build the common-sample AR regression used by forecast::ocsb.test:
      ! y_t on y_(t-1),...,y_(t-lag), with the first maxlag values dropped.
      real(dp),intent(in)::y(:)
      integer,intent(in)::lag,maxlag
      real(dp),allocatable,intent(out)::response(:),xreg(:,:)
      integer::start,nobs,t,j
      start=maxlag+1
      if(start>size(y))then
         allocate(response(0),xreg(0,lag))
         return
      end if
      nobs=size(y)-start+1
      allocate(response(nobs),xreg(nobs,lag))
      response=y(start:)
      if(lag>0)then
         do t=1,nobs
            do j=1,lag
               xreg(t,j)=y(start+t-1-j)
            end do
         end do
      end if
   end subroutine lag_regression_series

   subroutine fit_ocsb_lag(x,period,lag,maxlag,stat,rss,nobs,ncoef,info)
      real(dp),intent(in)::x(:)
      integer,intent(in)::period,lag,maxlag
      real(dp),intent(out)::stat,rss
      integer,intent(out)::nobs,ncoef,info
      real(dp),allocatable::ds(:),d1(:),ydd(:),resp(:),xlags(:,:),arcoef(:),arres(:)
      real(dp),allocatable::z4(:),z5(:),xf(:,:),yf(:),coef(:),resid(:),xtx(:,:),inv(:,:)
      logical,allocatable::valid4(:),valid5(:)
      real(dp)::sigma2,se
      integer::n,t,j,row,rank,k,tt,idx

      n=size(x)
      stat=huge(1.0_dp)
      rss=huge(1.0_dp)
      nobs=0
      ncoef=0
      info=1
      if(period<=1 .or. n<=period+maxlag+4)return

      ds=seasonal_difference(x,period,1)
      d1=difference_series(x,1)
      ydd=difference_series(ds,1)
      call lag_regression_series(ydd,lag,maxlag,resp,xlags)
      if(size(resp)<=lag+2)return

      if(lag>0)then
         call least_squares(xlags,resp,arcoef,arres,rank,info)
         if(info/=0)return
      else
         allocate(arcoef(0),arres(size(resp)))
         arres=resp
      end if

      allocate(z4(n),z5(n),valid4(n),valid5(n))
      z4=0.0_dp
      z5=0.0_dp
      valid4=.false.
      valid5=.false.

      ! ds index i corresponds to original time t=period+i.
      do idx=lag+1,size(ds)
         tt=period+idx
         z4(tt)=ds(idx)
         do j=1,lag
            z4(tt)=z4(tt)-arcoef(j)*ds(idx-j)
         end do
         valid4(tt)=.true.
      end do
      ! d1 index i corresponds to original time t=1+i.
      do idx=lag+1,size(d1)
         tt=1+idx
         z5(tt)=d1(idx)
         do j=1,lag
            z5(tt)=z5(tt)-arcoef(j)*d1(idx-j)
         end do
         valid5(tt)=.true.
      end do

      ! ydd index i corresponds to original time t=period+1+i.
      row=0
      do idx=maxlag+1,size(ydd)
         tt=period+1+idx
         if(tt-1<1 .or. tt-period<1)cycle
         if(.not.valid4(tt-1) .or. .not.valid5(tt-period))cycle
         row=row+1
      end do
      nobs=row
      ncoef=lag+2
      if(nobs<=ncoef)return
      allocate(yf(nobs),xf(nobs,ncoef))
      row=0
      do idx=maxlag+1,size(ydd)
         tt=period+1+idx
         if(tt-1<1 .or. tt-period<1)cycle
         if(.not.valid4(tt-1) .or. .not.valid5(tt-period))cycle
         row=row+1
         yf(row)=ydd(idx)
         do j=1,lag
            xf(row,j)=ydd(idx-j)
         end do
         xf(row,lag+1)=z4(tt-1)
         xf(row,lag+2)=z5(tt-period)
      end do
      call least_squares(xf,yf,coef,resid,rank,info)
      if(info/=0)return
      rss=sum(resid**2)
      sigma2=rss/real(max(1,nobs-ncoef),dp)
      xtx=matmul(transpose(xf),xf)
      call inverse_matrix(xtx,inv,info)
      if(info/=0)return
      se=sqrt(max(sigma2*inv(ncoef,ncoef),tiny(1.0_dp)))
      stat=coef(ncoef)/se
      info=0
   end subroutine fit_ocsb_lag

   subroutine ocsb_test(y,period,maxlag,lag_method,statistic,critical,lag_order,info)
      real(dp),intent(in)::y(:)
      integer,intent(in)::period
      integer,intent(in),optional::maxlag
      character(len=*),intent(in),optional::lag_method
      real(dp),intent(out)::statistic,critical
      integer,intent(out)::lag_order,info
      integer::ml,lag,nobs,ncoef,istat,best_lag
      real(dp)::stat,rss,ic,bestic,ll
      character(len=8)::method

      ml=0
      if(present(maxlag))ml=max(0,maxlag)
      method='fixed'
      if(present(lag_method))method=adjustl(lag_method)
      critical=ocsb_critical_value(period)
      best_lag=ml
      if(ml>0 .and. trim(method)/='fixed' .and. trim(method)/='FIXED')then
         bestic=huge(1.0_dp)
         best_lag=0
         do lag=1,ml
            call fit_ocsb_lag(y,period,lag,ml,stat,rss,nobs,ncoef,istat)
            if(istat/=0 .or. rss<=0.0_dp)cycle
            ll=-0.5_dp*real(nobs,dp)*(log(2.0_dp*pi)+1.0_dp+log(rss/real(nobs,dp)))
            select case(trim(method))
            case('AIC','aic')
               ic=-2.0_dp*ll+2.0_dp*real(ncoef+1,dp)
            case('BIC','bic')
               ic=-2.0_dp*ll+log(real(nobs,dp))*real(ncoef+1,dp)
            case('AICc','aicc','AICC')
               ic=-2.0_dp*ll+2.0_dp*real(ncoef+1,dp)
               if(nobs>ncoef+2)then
                  ic=ic+2.0_dp*real((ncoef+1)*(ncoef+2),dp)/real(nobs-ncoef-2,dp)
               else
                  ic=huge(1.0_dp)
               end if
            case default
               ic=-2.0_dp*ll+2.0_dp*real(ncoef+1,dp)
            end select
            if(ic< bestic)then
               bestic=ic
               ! forecast::ocsb.test uses which.min(fits)-1 after fitting lags 1:maxlag.
               best_lag=lag-1
            end if
         end do
      end if
      lag_order=best_lag
      call fit_ocsb_lag(y,period,best_lag,best_lag,statistic,rss,nobs,ncoef,info)
   end subroutine ocsb_test

   integer function nsdiffs(y,m,max_D,test,maxlag) result(Ds)
      real(dp),intent(in)::y(:)
      integer,intent(in)::m
      integer,intent(in),optional::max_D,maxlag
      character(len=*),intent(in),optional::test
      real(dp),allocatable::z(:),sdiff(:)
      real(dp)::r,strength,v0,v1,stat,crit
      integer::md,ml,lagord,info
      character(len=8)::method
      logical::dodiff

      md=1
      if(present(max_D))md=max_D
      ml=3
      if(present(maxlag))ml=maxlag
      method='seas'
      if(present(test))method=adjustl(test)
      Ds=0
      if(m<=1 .or. size(y)<2*m)return
      z=y
      do while(Ds<md .and. size(z)>2*m)
         if(trim(method)=='ocsb' .or. trim(method)=='OCSB')then
            call ocsb_test(z,m,ml,'AIC',stat,crit,lagord,info)
            dodiff=(info==0 .and. stat>crit)
         else
            strength=seasonal_strength(z,m)
            dodiff=(strength>0.64_dp)
         end if
         if(.not.dodiff)exit
         z=seasonal_difference(z,m,1)
         Ds=Ds+1
      end do
   end function nsdiffs
end module forecast_unitroot
